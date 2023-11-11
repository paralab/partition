import random
import pprint
from operator import itemgetter
import networkx as nx
import matplotlib.pyplot as plt

pp = pprint.PrettyPrinter(indent=2)
random.seed(42)
def display(x):
    print("")
    pp.pprint(x)
    print("")
# display = pp.pprint
def xyz_to_morton_index (coord):
    z,y,x = coord
    answer = 0
    for i in range(0,32//3):
        answer = answer | ((x & (1 << i)) << (2*i)) | ((y & (1 << i)) << (2*i+1)) | ((z & (1 << i)) << (2*i+2))

    return answer


def xy_to_morton_index (coord):
    x, y = coord
    answer = 0
    for i in range(0,32//2):
        answer = answer | ((x & (1 << i)) << (2*i)) | ((y & (1 << i)) << (2*i+1))

    return answer

N=5
data_points = []
# for i in range(N):
#     for j in range(N):
#         for k in range(N):
#             data_points.append([i,j,k])


for i in range(N):
    for j in range(N):
        x = random.uniform(-1,1)
        y = random.uniform(-1,1)
        data_points.append([x,y])


display(data_points)


all_coord_values = []
for idx, d in enumerate(data_points):
    # keep pair [value,location]
    all_coord_values.append([
        d[0],[idx,0]
    ])
    all_coord_values.append([
        d[1],[idx,1]
    ])
display(all_coord_values)
# for i in range(N*N*N):
#     print(xyz_to_morton_index(data_points[i]))


all_coord_values_sorted = sorted(all_coord_values, key=itemgetter(0))
display(all_coord_values_sorted)
data_points_adjusted = [[-1,-1] for _ in range(len(data_points))]
for i, x in enumerate(all_coord_values_sorted):
    data_points_adjusted[x[1][0]][x[1][1]] = i

display(data_points_adjusted)

morton_indices_temp = [[xy_to_morton_index(data),i] for i,data in enumerate(data_points_adjusted)]

# for i in range(N*N):
#     print(xy_to_morton_index(data_points_adjusted[i]))

display(morton_indices_temp)
morton_indices_temp_sorted = sorted(morton_indices_temp, key=itemgetter(0))
display(morton_indices_temp_sorted)

morton_indices_final = [-1 for _ in range((len(data_points)))]

for i,x in enumerate(morton_indices_temp_sorted):
    morton_indices_final[x[1]] = i

display(morton_indices_final)


G = nx.Graph()

for i,_ in enumerate(data_points):
    G.add_node(i)

plt.figure(figsize=(20,20))
nx.draw_networkx(G, data_points)