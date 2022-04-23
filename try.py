'''This file is for testing don't delete '''

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

import math

from node import Node
from environment import *
import heapq


def astar(grid, start, goal):
    path, node = [], []
    for i, row in enumerate(grid):
        col = []
        for j, val in enumerate(row):
            col.append(Node(i, j, val, 0))
        node.append(col)

    for i in range(len(node)):
        for j in range(len(node[i])):
            node[i][j].h = math.hypot(goal[0]-node[i][j].x, goal[1]-node[i][j].y)  # h value for each node

    seen = [[False for _ in range(len(node))]
            for _ in range(len(grid))]

    # initialize the start
    node[start[0]][start[1]].g = 0
    node[start[0]][start[1]].cost = node[start[0]][start[1]].g + node[start[0]][start[1]].h
    current_node = not_seen_min_cost(node, seen)
    seen[current_node.x][current_node.y] = True
    path_length = 0

    while current_node:
        # goal condition
        if current_node.x == goal[0] and current_node.y == goal[1]:
            path.append(goal)
            print(f'Current goal is {goal}')
            while path[-1] != start:
                print(f'Im here {current_node.x, current_node.y}')
                path_length += distance(current_node, current_node.parent)
                path.append((current_node.parent.x, current_node.parent.y))
                current_node = current_node.parent
            path.reverse()
            break

        neighbor_locations = travel(current_node, grid, seen)
        for neighbor in neighbor_locations:
            # check f value
            if node[neighbor[0]][neighbor[1]].g is None:
                node[neighbor[0]][neighbor[1]].g = math.inf
            temp_g = current_node.g + distance(current_node, node[neighbor[0]][neighbor[1]])
            if temp_g < node[neighbor[0]][neighbor[1]].g:
                node[neighbor[0]][neighbor[1]].parent = current_node
                node[neighbor[0]][neighbor[1]].g = temp_g
                node[neighbor[0]][neighbor[1]].cost = node[neighbor[0]][neighbor[1]].g + node[neighbor[0]][
                    neighbor[1]].h
        current_node = not_seen_min_cost(node, seen)
        # if no path exists
        if current_node is None:
            break
        seen[current_node.x][current_node.y] = True  # won't visit this node again

    return path, path_length
