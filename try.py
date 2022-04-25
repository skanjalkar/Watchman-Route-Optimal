'''This file is for testing don't delete '''
from shapely.geometry import Polygon
from shapely.geometry import Point
# Poly = [(4,4),(8,4),(8,0),(14,-5),(20,0),(20,6),(15,6),(15,10),\
#     (20,10),(20,14),(16,14),(16,16),(10,16),(10,14),(6,14),(6,16),(2,16)\
#        ,(2,14),(0,14),(-5,7),(0,0),(2,-2),(4,0)]
# Poly = [(0,0),(10,0),(10,1),(10,5),(9,5),(8,1),(8,5),(7,5)\
#         ,(6,1),(6,5),(5,5),(4,1),(4,5),(3,5),(2,1),(2,5),(1,5),(0,1)]
# Poly = [(250,190),(236,192),(207,221),(227,241),(193,275)\
#     ,(173,255),(155,272),(155,302),(164,302),(164,315)\
#     ,(206,315),(206,337),(233,337),(233,311),(257,311)\
#     ,(257,414),(167,414),(167,394),(100,394),(100,414)\
#     ,(43,414),(43,394),(13,394),(13,313),(35,313)\
#     ,(35,343),(62,343),(62,290),(47,290),(47,265)\
#     ,(20,265),(20,286),(0,286),(0,211),(119,211),(124,216)\
#     ,(168,172),(146,149),(164,131),(144,112),(126,130)\
#     ,(90,94),(114,70),(114,48),(195,48),(195,92),(262,92)\
#     ,(262,70),(320,70),(320,97),(365,97),(365,150)\
#     ,(343,150),(343,128),(314,128),(314,192),(343,192)\
#     ,(343,170),(374,170),(374,234),(343,234),(343,260)\
#     ,(283,260),(283,242),(249,242)]


# SP = Shrink_Polygon_AGP
#
# row,col = zip(*P)
# row = list(row)
# col = list(col)
# row.append(row[0])
# col.append(row[0])
# ax.plot(row, col)
# P.reverse()  # already clockwise
# P.append(P[0])
# points = SP.shrink(P)
# points.append(points[0])
# Pb = SP.create_point_pair(P)
# Yx = SP.non_intersecting_diag(points, P, Pb)
# Yn = SP.mini_chk_pts(Pb, points, P, Yx)
# Final_Diagonals = SP.clean_up_final(Yn)
# guards = SP.Guards(Final_Diagonals)
# Guards = tuple(tuple(map(int, tup)) for tup in guards)





''' 
This code is used for random polygon testing
    # P = draw_polygon(ax, 8, lim_x,  lim_y)
    # vor = Voronoi(P)
    # for i in vor.vertices:
    #     if i[0] < 0 or i[1] < 0 or i[0] > lim_x or i[1] > lim_y:
    #         continue
    #     x = int(i[0])
    #     y = int(i[1])
    #     vor_int.append((x, y))
'''
'''
lim_x, lim_y = 150, 150
    xr = random.randint(5,10)
    P = draw_polygon(ax, xr, lim_x,  lim_y)
    vor = Voronoi(P)
    final_guards = []
    for i in vor.vertices:
        if i[0] < 0 or i[1] < 0 or i[0] > lim_x or i[1] > lim_y:
            continue
        x = int(i[0])
        y = int(i[1])
        final_guards.append((x, y))
        '''

# print(track)
    # print()
    # print(back_track)
    # print(time.time()-s)
    # point_list = [i for i in combinations(watchman_route_pts, 2)]
    # print(point_list)

# answer = []
#
#
# def TSP(graph, v, current_pos, num_nodes, count, cost):   # brute force with recursion
#     if count == num_nodes and graph[current_pos][0]:
#         answer.append(cost + graph[current_pos][0])
#         return
#
#     for i in range(num_nodes):
#         if v[i] == False and graph[current_pos][i]:
#             v[i] = True
#             TSP(graph, v, i, num_nodes, count+1, cost + graph[current_pos][i])
#             v[i] = False