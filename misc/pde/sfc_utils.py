import math

def to_integer_coords(data_points):
    levels = 19

    # regular cube box method

    bounding_box = [math.inf,-math.inf]         # -limit to +limit
    for d in data_points:
        bounding_box[0] = min(min(d),bounding_box[0])
        bounding_box[1] = max(max(d),bounding_box[1])
    
    leaf_node_length = (bounding_box[1] - bounding_box[0])/(2**levels)

    integer_data_points = []

    for d in data_points:
        d_new = []
        for dim_i in range(3):
            d_new.append(int((d[dim_i]-bounding_box[0])/leaf_node_length))
        integer_data_points.append(d_new)
    return integer_data_points


def morton_compare(p1_data, p2_data):
    p1 = p1_data[0]
    p2 = p2_data[0]
    if p1[0]==p2[0] and p1[1]==p2[1] and p1[2]==p2[2]:      # same point
        return 0
    temp_x = p1[0] ^ p2[0]
    temp_y = p1[1] ^ p2[1]
    temp_z = p1[2] ^ p2[2]

    maxC = temp_z
    yOrx = temp_y

    if (yOrx < temp_x):
        if ((temp_x ^ yOrx) >= yOrx):
            yOrx = temp_x

    if (maxC < yOrx):
        if ((maxC ^ yOrx) >= maxC):
            maxC = yOrx
        
    if (maxC == temp_z):
        return -1 if p1[2] < p2[2] else 1
    elif (maxC == temp_y):
        return -1 if p1[1] < p2[1] else 1
    else:
        return -1 if p1[0] < p2[0] else 1