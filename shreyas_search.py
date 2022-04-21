import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import time
import math
from queue import Queue
from matplotlib.animation import FuncAnimation
from itertools import combinations
from itertools import pairwise
from itertools import permutations
import sys
global answer
answer = []


class Node:
    def __init__(self, x, y, inside, h):
        self.x = x
        self.y = y
        self.inside = inside
        self.h = h
        self.g = 0
        self.cost = math.inf
        self.parent = None


def grid(limx, limy, polygon):
    grid = []
    for i in range(0, limx+1):
        row = []
        for j in range(0, limy+1):
            if polygon.contains(Point(i, j)):
                row.append(1)
            else:
                row.append(0)
        grid.append(row)
    return grid


def mark_visited(neighbours, visited):
    for (x, y) in neighbours:
        visited[x][y] = True


def valid(x, y, grid, visited):
    if (x >= 0 and y >= 0) and (x < len(grid) and y < len(grid[0])):
        if grid[x][y]==1 and visited[x][y] is False:
            return True
    return False

def travel(node, grid, visited):
    possible_travel = []
    if valid(node.x, node.y+1, grid, visited):  # move right
        possible_travel.append([node.x, node.y+1])

    if valid(node.x+1, node.y, grid, visited):  # move down
        possible_travel.append([node.x+1, node.y])

    if valid(node.x, node.y-1, grid, visited):  # move left
        possible_travel.append([node.x, node.y-1])

    if valid(node.x-1, node.y, grid, visited):  # move up
        possible_travel.append([node.x-1, node.y])

    if valid(node.x-1, node.y-1, grid, visited):
        possible_travel.append([node.x-1, node.y-1])

    if valid(node.x+1, node.y-1, grid, visited):  # move up
        possible_travel.append([node.x+1, node.y-1])

    if valid(node.x-1, node.y+1, grid, visited):  # move up
        possible_travel.append([node.x-1, node.y+1])

    if valid(node.x+1, node.y+1, grid, visited):  # move up
        possible_travel.append([node.x+1, node.y+1])
    return possible_travel  # returns all possible node locations that you can travel to


def not_seen_min_cost(node, seen):
    min_node = None  # initialize min node
    for i in range(len(node)):
        for j in range(len(node[i])):
            if seen[i][j] or node[i][j].cost is None:  # if true
                continue
            if min_node is None:  # updating node
                min_node = node[i][j]
            elif min_node.cost > node[i][j].cost:  # cost
                min_node = node[i][j]
            elif min_node.cost == node[i][j].cost:
                if min_node.h > node[i][j].h:  # comparing h when cost is same
                    min_node = node[i][j]

    return min_node

def distance(start, end):
    distance = math.hypot(abs(end.x-start.x), abs(end.y - start.y))
    return distance


def search(grid, start, goal, poly):
    path = []
    node = []
    path_length = 0

    # print(f'this is {start}, {goal}')
    # created nodes
    for i, row in enumerate(grid):
        col = []
        for j, val in enumerate(row):
            col.append(Node(i, j, val, 0))
        node.append(col)
    # 2d array to keep track of visited/seen
    seen = [[False for _ in range(len(node))]
            for _ in range(len(grid))]

    # initialize starting node
    node[start[0]][start[1]].cost = 0
    current_node = not_seen_min_cost(node, seen)
    # print(current_node.x, current_node.y)
    seen[current_node.x][current_node.y] = True

    while current_node:
        # goal condition
        if current_node.x == goal[0] and current_node.y == goal[1]:
            path.append(goal)
            # print(current_node.x, current_node.y)
            while path[-1] != start:
                path_length += distance(current_node, current_node.parent)
                path.append((current_node.parent.x, current_node.parent.y))
                current_node = current_node.parent
            path.reverse()
            break

        neighbor_locations = travel(current_node, grid, seen)
        for neighbor in neighbor_locations:

            # relax condition
            if node[neighbor[0]][neighbor[1]].cost is None or node[neighbor[0]][
                neighbor[1]].cost > current_node.cost + 1:
                node[neighbor[0]][neighbor[1]].cost = current_node.cost + 1
                node[neighbor[0]][neighbor[1]].parent = current_node

        current_node = not_seen_min_cost(node, seen)  # get current node with min cost
        if current_node is None:  # if no path exists
            break
        seen[current_node.x][current_node.y] = True  # won't visit this element again


    return path, path_length

def astar(grid, start, goal, poly):
    path = []
    node = []
    for i, row in enumerate(grid):
        col = []
        for j, val in enumerate(row):
            col.append(Node(i, j, val, 0))
        node.append(col)

    for i in range(len(node)):
        for j in range(len(node[i])):
            node[i][j].h = abs(goal[0] - node[i][j].x) + abs(goal[1] - node[i][j].y)  # h value for each node

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
            while path[-1] != start:
                path_length += distance(current_node, current_node.parent)
                path.append((current_node.parent.x, current_node.parent.y))
                current_node = current_node.parent
            path.reverse()
            break

        neighbor_locations = travel(current_node, grid, seen)
        for neighbor in neighbor_locations:
            # check f value
            if node[neighbor[0]][neighbor[1]].cost is None or \
                    node[neighbor[0]][neighbor[1]].cost > current_node.g + 1 + node[neighbor[0]][neighbor[1]].h:
                node[neighbor[0]][neighbor[1]].g = current_node.g + 1  # 1 is the cost to move
                node[neighbor[0]][neighbor[1]].cost = node[neighbor[0]][neighbor[1]].g + node[neighbor[0]][
                    neighbor[1]].h
                node[neighbor[0]][neighbor[1]].parent = current_node

        current_node = not_seen_min_cost(node, seen)
        # if no path exists
        if current_node is None:
            break
        seen[current_node.x][current_node.y] = True  # won't visit this node again

    return path, path_length

def animate(i, path):
    x = []
    y = []
    for i in path:
        x.append(i[0])
        y.append(i[1])
    ax.clear()
    ax.plot(x, y)



def draw_polygon(ax, n, lim_x, lim_y):

    x = np.random.randint(0, lim_x, n)
    y = np.random.randint(0, lim_y, n)

    # computing the (or a) 'center point' of the polygon
    center_point = [np.sum(x)/n, np.sum(y)/n]

    angles = np.arctan2(x-center_point[0],y-center_point[1])

    # sorting the points:
    sort_tups = sorted([(i, j, k) for i, j, k in zip(x, y, angles)], key=lambda t: t[2])

    # making sure that there are no duplicates:
    if len(sort_tups) != len(set(sort_tups)):
        raise Exception('two equal coordinates -- exiting')

    x, y, angles = zip(*sort_tups)
    x = list(x)
    y = list(y)

    # appending first coordinate values to lists:
    x.append(x[0])
    y.append(y[0])

    ax.plot(x, y, label='{}'.format(n))
    poly = list(zip(x, y))
    return poly

def TSP(graph, v, current_pos, num_nodes, count, cost):   # brute force with recursion
    if count == num_nodes and graph[current_pos][0]:
        answer.append(cost + graph[current_pos][0])
        return

    for i in range(num_nodes):
        if v[i] == False and graph[current_pos][i]:
            v[i] = True
            TSP(graph, v, i, num_nodes, count+1, cost + graph[current_pos][i])
            v[i] = False

def held_karp(paths):
    n = len(paths)
    C = {}
    for k in range(1, n):
        pass

def brute_force(graph):   # naive brute force
    vertex = []
    for i in range(len(graph)):
        if i != 0:
            vertex.append(i)
    # print(vertex)

    min_path = sys.maxsize
    next_permutation = permutations(vertex)
    # next_permutation = [perm for perm in permutations(vertex)]
    # print(next_permutation)
    for i in next_permutation:
        # print(i)
        current_length = 0
        k = 0
        for j in i:
            current_length += graph[k][j]
            k = j
        current_length += graph[k][0]

        min_path = min(min_path, current_length)

    return min_path

# def pairwise(it):
#     it = iter(it)
#     while True:
#         try:
#             yield next(it), next(it)
#         except StopIteration:
#             # no more elements in the iterator
#             return


if __name__ == '__main__':
    fig, ax = plt.subplots()
    lim_x, lim_y = 40, 40
    poly = draw_polygon(ax, 10, lim_x,  lim_y)
    vor = Voronoi(poly)
    vor_int, vor_x, vor_y, vor_inside_poly = [], [], [], []
    for i in vor.vertices:
        if i[0] < 0 or i[1] < 0 or i[0] > lim_x or i[1] > lim_y:
            continue
        x = int(i[0])
        y = int(i[1])
        vor_int.append((x, y))
    polygon = Polygon(poly)
    grid = grid(lim_x, lim_y, polygon)
    for i in vor_int:
        if polygon.contains(Point(i)):
            vor_x.append(int(i[0]))
            vor_y.append(int(i[1]))
            vor_inside_poly.append(i)
    print(len(vor_inside_poly))
    track = {}
    back_track = []
    paths = []

    s = time.time()
    for i in range(len(vor_inside_poly)):
        track[vor_inside_poly[i]] = {}
        f = []
        for j in range(len(vor_inside_poly)):
            # path, path_length = search(grid, vor_inside_poly[i], vor_inside_poly[j], polygon)
            path, path_length = astar(grid, vor_inside_poly[i], vor_inside_poly[j], polygon)
            f.append(path_length)
            track[vor_inside_poly[i]][vor_inside_poly[j]] = path_length
            back_track.append(path)
        paths.append(f)
    # print(track)
    # print(back_track)
    # print(vor_inside_poly)
    v = [False for i in range(len(vor_inside_poly))]
    v[0] = True
    # print(paths)
    # print(paths)
    TSP(paths, v, 0, len(vor_inside_poly), 1, 0)
    print(brute_force(paths))
    # plotTSP(paths, vor_inside_poly, len(paths))
    print(answer)
    print(min(answer), answer.index(min(answer)))
    # print(track)
    # print()
    # print(back_track)
    # print(time.time()-s)
    # point_list = [i for i in combinations(vor_inside_poly, 2)]
    # print(point_list)
    for i in back_track:
        if len(i) == 1:
            continue
        point_list = [pair for pair in pairwise(i)]
        for j in range(len(point_list)):
            x_values = [point_list[j][0][0], point_list[j][1][0]]
            y_values = [point_list[j][0][1], point_list[j][1][1]]
            plt.plot(x_values, y_values, c='grey')
    # ani = FuncAnimation(fig, animate(back_track), frames=30, interval=500, repeat=False)
    plt.scatter(vor_x, vor_y, c="red")
    plt.show()

    plt.scatter(vor_x, vor_y, c='blue')

