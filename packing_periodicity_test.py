#!/usr/bin/python3

import argparse
import math
import sys
import re

## 目前这个脚本验证周期性是可以适用于任何周期性情况，但是没有专门区分，虽然有 periodic_axes 命令行输入 但是没有用到具体的逻辑里
## 需要同时满足三个检测方法的精度都保持在 1e-4


def read_beads(file_path):
    beads = {}
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            x, y, z, d = map(float, line.split())
            beads[i] = {'x': x, 'y': y, 'z': z, 'radius': d/2}
    return beads


def calculate(bead_info, axis, value):
    if axis == 'corner':
        # 如果 axis 是 'corner'，则 value 应为一个元组表示的角落坐标 (x, y, z)
        corner_x, corner_y, corner_z = value
        # 计算球心到角落的距离
        distance = math.sqrt((bead_info['x'] - corner_x)**2 + (bead_info['y'] - corner_y)**2) ## 只考虑了x和y
        if distance > bead_info['radius']:
            return None
        # 计算投影半径
        ## projected_radius = math.sqrt(bead_info['radius']**2 - distance**2)
        ## 这个计算没有什么意义 可以直接用 radius - distance表示
        projected_radius = bead_info['radius'] - distance ## stick out length difference
        
        inter_data = bead_info.copy()
        inter_data['radius'] = projected_radius
        inter_data['cutplane_dis'] = {
            'x': corner_x, 'y': corner_y, 'z': corner_z,  # 包括z轴信息
            'distance': distance,
            'projected_radius': projected_radius
        }
        
        return {'cutplane_dis': inter_data}
    else:
        # 处理普通轴的情况
        distance = abs(bead_info[axis] - value)
        if distance > bead_info['radius']:
            return None
        
        projected_radius = math.sqrt(bead_info['radius']**2 - distance**2)
        
        inter_data = bead_info.copy()
        inter_data['radius'] = projected_radius
        
        return {'cutplane_dis': inter_data}





def find_intersections(beads, cut_bounds):
    """
    Find out which cut plane each bead intersects, including diagonal intersections on the XY plane.
    Record the full 3D coordinates for corners to provide complete spatial information.
    cut_bounds: positions of the cut planes [x_min, y_min, z_min, x_max, y_max, z_max]
    """
    for bead_id, bead_info in beads.items():
        intersections = []

        # Only consider x and y for periodic boundaries (assuming z is not periodic)
        for i, (axis, min_val, max_val) in enumerate(zip('xy', cut_bounds[:2], cut_bounds[3:5])):
            if bead_info[axis] - bead_info['radius'] < min_val:
                inter = calculate(bead_info, axis, min_val)
                if inter:
                    intersections.append({'cutplane': f"{axis}={min_val}", **inter})
                    
            if bead_info[axis] + bead_info['radius'] > max_val:
                inter = calculate(bead_info, axis, max_val)
                if inter:
                    intersections.append({'cutplane': f"{axis}={max_val}", **inter})

        # Check diagonal intersections only in the XY plane, but include the Z coordinate in the corner information
        corners = [
            (cut_bounds[0], cut_bounds[1], bead_info['z']),  # min_x, min_y, current bead z
            (cut_bounds[3], cut_bounds[4], bead_info['z']),  # max_x, max_y, current bead z
            (cut_bounds[0], cut_bounds[4], bead_info['z']),
            (cut_bounds[3], cut_bounds[1], bead_info['z'])

        ]
        for corner in corners:
            dist_to_corner = ((bead_info['x'] - corner[0]) ** 2 + (bead_info['y'] - corner[1]) ** 2) ** 0.5
            if dist_to_corner < bead_info['radius']:
                inter = calculate(bead_info, 'corner', corner)
                if inter:
                    intersections.append({'cutplane': f"corner={corner}", **inter})

        if intersections:
            bead_info['intersections'] = intersections
    return beads


def match_counterparts(beads, cut_bounds, periodic_axes):
    for bead_id, bead_info in beads.items():
        bead_info['counterbeads'] = []

        for other_bead_id, other_bead_info in beads.items():
            if bead_id == other_bead_id:
                continue

            # Check if the positions and intersections match across PBC
            match_found = True

            if match_found and 'intersections' in bead_info and 'intersections' in other_bead_info:
                matched_dirs = set()
                for intersection in bead_info['intersections']:
                    matching_intersection = False
                    
                    # Check if this is a corner intersection on the XY plane
                    if 'corner' in intersection['cutplane']:
                        axis = 'corner'
                        corners = [
                            (cut_bounds[0], cut_bounds[1], bead_info['z']),  # min_x, min_y, current bead z
                            (cut_bounds[3], cut_bounds[4], bead_info['z']),  # max_x, max_y, current bead z
                            (cut_bounds[0], cut_bounds[4], bead_info['z']),
                            (cut_bounds[3], cut_bounds[1], bead_info['z'])
                        ]
                        corner = tuple(map(float, intersection['cutplane'].split('=')[1].strip('()').split(',')))


                        
                        if corner == corners[0]:  
                            opposite_corner = corners[1]  
                        elif corner == corners[1]: 
                            opposite_corner = corners[0]  
                        elif corner == corners[2]: 
                            opposite_corner = corners[3] 
                        elif corner == corners[3]: 
                            opposite_corner = corners[2]
                        else:
                            opposite_corner = None 

                        if opposite_corner:
                            opposite_cutplane = f"corner={opposite_corner}"
                        else:
                            print("No valid corner match found.")
                            opposite_cutplane = None

                    else:
                        axis = intersection['cutplane'][0]
                        cutplane = intersection['cutplane']
                        if '={}'.format(cut_bounds[0]) in cutplane:
                            opposite_cutplane = cutplane.replace('={}'.format(cut_bounds[0]), '={}'.format(cut_bounds[3]))
                        else:
                            opposite_cutplane = cutplane.replace('={}'.format(cut_bounds[3]), '={}'.format(cut_bounds[0]))
                            ## 这里不考虑x,y是因为正方形长度一样 不然应该分别考虑xy

                    # Check against intersections of the other bead
                    for other_dir, other_intersection in enumerate(other_bead_info['intersections']):
                        if other_dir in matched_dirs:
                            continue
                        if opposite_cutplane == other_intersection['cutplane']:
                            if axis == 'corner':
                                # Special handling for corner cases, assume match found directly
                                if (intersection['cutplane_dis']['radius'] - other_intersection['cutplane_dis']['radius'] < 1e-4):
                                    matching_intersection = True
                                    matched_dirs.add(other_dir)
                                    bead_info['counterbeads'].append((other_bead_id, f"corner match at {opposite_cutplane}"))
                            else:
                                all_axes = ['x', 'y', 'z']
                                all_axes.remove(axis)
                                if all(abs(intersection['cutplane_dis'][a] - other_intersection['cutplane_dis'][a]) < 1e-4 for a in all_axes if a in intersection['cutplane_dis']) and \
                                    abs(intersection['cutplane_dis']['radius'] - other_intersection['cutplane_dis']['radius']) < 1e-4:
                                    matching_intersection = True
                                    matched_dirs.add(other_dir)
                                    bead_info['counterbeads'].append((other_bead_id, opposite_cutplane))
                            break

                    if not matching_intersection:
                        match_found = False

    return beads




def sort_counterparts(beads):
    for bead_id, bead_info in beads.items():
        if 'counterbeads' in bead_info:
            # sort the list of counterbeads, based on the first character (x, y, z) of the cutplane information
            bead_info['counterbeads'] = sorted(bead_info['counterbeads'], key=lambda x: x[1][0])
    return beads

def print_counterparts(beads, output_file):
    with open(output_file, 'w') as f:
        for bead_id, bead_info in beads.items():
            direction_to_id = {'x': None, 'y': None, 'z': None}
        
            if 'counterbeads' in bead_info:
                for counterpart in bead_info['counterbeads']:
                    counterpart_id, cutplane = counterpart
                    direction = cutplane[0]  # get direction (x, y, z)
                    direction_to_id[direction] = (counterpart_id, cutplane)
        
            counterparts_str = ', '.join([f"{dir}:{id_and_plane[0]} ({id_and_plane[1]})" for dir, id_and_plane in direction_to_id.items() if id_and_plane is not None])
            
            intersections_str = ', '.join([f"{intersection['cutplane']}" for intersection in bead_info.get('intersections', [])])
        
            if counterparts_str or intersections_str:
                f.write(f"Bead {bead_id} counterparts: {counterparts_str if counterparts_str else 'None'}, intersections: {intersections_str if intersections_str else 'None'}\n")
            else:
                f.write(f"Bead {bead_id} has no counterparts or intersections.\n")



def parse_details(detail_str):
    detail_dict = {}
    details_order = []
    items = detail_str.split(', ')
    for item in items:
        if ':' in item:
            key, value = item.split(':')
            key = key.strip().split(' ')[-1]
            value = value.strip()
            if '(' in value:
                # 提取括号中的数据并转换为字典
                detail_match = re.findall(r'(\w+)=(\d+\.?\d*)', value)
                value_dict = {k: float(v) if '.' in v else int(v) for k, v in detail_match}
                detail_dict[key] = value_dict
            elif 'corner match' in value:
                # 处理corner数据
                corner_data = re.search(r'corner=\((.*?)\)', value).group(1)
                corner_values = tuple(map(float, corner_data.split(', ')))
                detail_dict[key] = corner_values
            else:
                detail_dict[key] = int(value)
            details_order.append(int(re.search(r'\d+', value.split(' ')[0]).group()))
    return detail_dict, details_order




'''
# def check_beads(file_path, output_file):
#     input_start = input("请输入起始珠子编号（留空检查所有珠子）: ")
#     input_end = input("请输入结束珠子编号（留空检查所有珠子）: ")
    
#     start_bead = int(input_start) if input_start.strip() else None
#     end_bead = int(input_end) if input_end.strip() else None

#     errors = []
#     first_bead_detail = None
#     last_bead_detail = None
#     previous_ids = None

#     with open(file_path, 'r') as file:
#         lines = file.readlines()

#     for line in lines:
#         bead_number = int(line.split()[1])
#         if (start_bead is None or bead_number >= start_bead) and (end_bead is None or bead_number <= end_bead):
#             if 'counterparts:' in line and 'intersections:' in line:
#                 parts = line.split('counterparts:')
#                 second_part = parts[1].split('intersections:')
#                 counterparts, ids_order = parse_details(second_part[0].strip())
#                 intersections, _ = parse_details(second_part[1].strip())

#                 if first_bead_detail is None:
#                     first_bead_detail = ids_order
#                 last_bead_detail = ids_order

#                 if previous_ids and previous_ids[-1] + 1 != ids_order[0]:
#                     errors.append(f"Bead {bead_number}: Counterpart IDs are not consecutive. Expected start {previous_ids[-1] + 1}, found {ids_order[0]}.")

#                 for key, counterpart_details in counterparts.items():
#                     if isinstance(counterpart_details, dict):  # 处理括号内的值
#                         for axis, counterpart_value in counterpart_details.items():
#                             if axis in intersections:
#                                 intersection_value = intersections[axis]
#                                 expected_value = 10 if counterpart_value == 0 else 0
#                                 if intersection_value != expected_value:
#                                     errors.append(f"Bead {bead_number}: Mismatch for {axis} in {key}, expected {expected_value}, found {intersection_value}.")
#                     else:
#                         if key in intersections:
#                             intersection_value = intersections[key]
#                             expected_value = 10 if counterpart_details == 0 else 0
#                             if intersection_value != expected_value:
#                                 errors.append(f"Bead {bead_number}: Mismatch for {key}, expected {expected_value}, found {intersection_value}.")
#                         else:
#                             errors.append(f"Bead {bead_number}: No intersection data for key '{key}'.")

#                 # 检查corner数据
#                 if 'corner' in intersections and 'corner' in counterparts:
#                     corner_intersections = intersections['corner']
#                     corner_counterparts = counterparts['corner']
#                     if isinstance(corner_intersections, tuple) and isinstance(corner_counterparts, tuple):
#                         # 检查x和y的互为0和1，z值相同
#                         if corner_intersections[0] != 10 - corner_counterparts[0]:
#                             errors.append(f"Bead {bead_number}: Mismatch for x in corner, expected {10 - corner_counterparts[0]}, found {corner_intersections[0]}.")
#                         if corner_intersections[1] != 10 - corner_counterparts[1]:
#                             errors.append(f"Bead {bead_number}: Mismatch for y in corner, expected {10 - corner_counterparts[1]}, found {corner_intersections[1]}.")
#                         if corner_intersections[2] != corner_counterparts[2]:
#                             errors.append(f"Bead {bead_number}: Mismatch for z in corner, expected {corner_counterparts[2]}, found {corner_intersections[2]}.")
#                     else:
#                         errors.append(f"Bead {bead_number}: Corner data format mismatch.")

#                 previous_ids = ids_order

#             else:
#                 if 'counterparts:' in line or 'intersections:' in line:
#                     errors.append(f"Bead {bead_number}: Incomplete data for counterparts or intersections.")

#     with open(output_file, 'a') as file:
#         file.write("\n\n--- New Check Results ---\n")
#         if errors:
#             for error in errors:
#                 file.write(error + '\n')
#         else:
#             file.write("All beads checks passed successfully.\n")
        
#         if first_bead_detail:
#             file.write(f"First bead counterparts start at: {first_bead_detail[0]}\n")
#         if last_bead_detail:
#             file.write(f"Last bead counterparts end at: {last_bead_detail[-1]}\n")

#     print(f"Results appended to {output_file}")
'''




def main():
    parser = argparse.ArgumentParser(description="Find and print bead counterparts.")
    parser.add_argument("-f","--file_path", type=str, help='Path to the beads_used.txt.')
    parser.add_argument("-b", "--cutbounds", type=float, nargs=6, metavar=('XMIN', 'YMIN', 'ZMIN', 'XMAX', 'YMAX', 'ZMAX'), help='Cutting plane bounds for the beads.')
    parser.add_argument("-p", "--periodic_axes", type=str, default="", help="Axes with periodic boundary conditions (e.g., 'x', 'xy', 'xyz')")

    ## default=[0, 0, 0, 6, 6, 6]

    args = parser.parse_args()

    int_cutbounds = [int(x) for x in args.cutbounds]
    periodic_axes = [axis for axis in args.periodic_axes]

    errors = [] 

    if not args.file_path:
        errors.append("The beads_used.txt file path is required. Use -f or --file_path to specify the path to the beads_used_xyzd file.")
    if not args.cutbounds or len(args.cutbounds) != 6:
        errors.append("Exactly six values are required for cutbounds. Use -b or --cutbounds to specify the bounds as six floats.")

    if errors:
        parser.print_usage()
        print(f"{parser.prog}: error: " + "\n".join(errors), file=sys.stderr)
        sys.exit(1)

    beads = read_beads(args.file_path)
    beads = find_intersections(beads, int_cutbounds)
    beads = match_counterparts(beads, int_cutbounds, periodic_axes)
    beads = sort_counterparts(beads)
    print_counterparts(beads, 'counterparts_results_v3.txt')

    # check_beads("counterparts_results_v3.txt", "counterparts_results_v3.txt")

    ## beads_216 = '/home/jovyan/Msc_thesis/data/beads_used.txt' ## file path arg like this 
    ## beads_216 = '/home/jovyan/Msc_thesis/data/packing.txt'
    ## cut_bounds = [0, 0, 0, 6, 6, 6]
    

    '''
    beads = read_beads(beads_216)
    beads = find_intersections(beads, cut_bounds)
    beads = match_counterparts(beads, cut_bounds)
    beads = sort_counterparts(beads)
    print_counterparts(beads)
    '''

if __name__ == "__main__":
    main()