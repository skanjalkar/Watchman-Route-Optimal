import itertools
import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry.polygon import Polygon
import time
import Shrink_Polygon_Watchman_Route
from astar import *
from dijkstra import *
from TSP_search import *

def get_min(polygon):
    min_x = 0
    min_y = 0

    for x, y in polygon:
        if x < min_x:
            min_x = x
        if y < min_y:
            min_y = y

    return min_x, min_y

def translate_poly(polygon, min_x, min_y):
    new_poly = []
    for x,y in polygon:
        new_poly.append((x+abs(min_x), y+abs(min_y)))
    return new_poly


if __name__ == '__main__':
    fig, ax = plt.subplots()
    lim_x, lim_y = 150, 150

    final_guards, guard_x, guard_y, watchman_route_pts = [], [], [], []
    guards, old_p = Shrink_Polygon_Watchman_Route.shrink()
    mx, my = get_min(old_p)
    Guards = translate_poly(guards, mx, my)
    P = translate_poly(old_p, mx, my)
    final_guards = tuple(tuple(map(int, tup)) for tup in Guards)
    polygon = Polygon(P)
    grid = grid(lim_x, lim_y, polygon)
    for i in final_guards:
        if polygon.contains(Point(i)):
            guard_x.append(int(i[0])), guard_y.append(int(i[1]))
            watchman_route_pts.append(i)
    print(len(f"Number of guards: {watchman_route_pts}"))
    back_track, paths, back_track_track = [], [], []
    s = time.time()
    for i in range(len(watchman_route_pts)):
        path_length_for_one_pt = []
        path_for_specific_point = []
        for j in range(len(watchman_route_pts)):
            # path, path_length = search(grid, watchman_route_pts[i], watchman_route_pts[j], polygon)
            path, path_length = astar(grid, watchman_route_pts[i], watchman_route_pts[j])
            path_length_for_one_pt.append(path_length)
            back_track.append(path)
            path_for_specific_point.append(path)
        back_track_track.append(path_for_specific_point)
        paths.append(path_length_for_one_pt)
    v = [False for i in range(len(watchman_route_pts))]
    v[0] = True
    print(paths)
    # TSP(paths, v, 0, len(watchman_route_pts), 1, 0)
    print(f'The path lenght from brute force {brute_force(paths)}')
    # plotTSP(paths, watchman_route_pts, len(paths))
    # print(answer)
    # print(min(answer), answer.index(min(answer)))
    pl, final_path, t = held_karp(paths)
    print(f'This is the order of traversal {final_path}')
    print(f'The path length for held-karp {pl}')

    for i in back_track:
        if len(i) == 1:
            continue
        point_list = [pair for pair in itertools.pairwise(i)]
        for j in range(len(point_list)):
            x_values = [point_list[j][0][0], point_list[j][1][0]]
            y_values = [point_list[j][0][1], point_list[j][1][1]]
            plt.plot(x_values, y_values, c='grey')

    x, y = zip(*P)
    plt.plot(x, y)
    plt.scatter(guard_x, guard_y, c='blue')
    plt.scatter(guard_x, guard_y, c="red")
    plt.show()


    x, y = zip(*P)
    plt.figure()
    plt.plot(x, y)
    plt.scatter(guard_x, guard_y, c='blue')
    final_arrange = []
    final_path.append(0)

    for i in range(len(final_path)-1):
        final_arrange.append(back_track_track[final_path[i]][final_path[i+1]])

    final_pairs = [[j for j in itertools.pairwise(final_pair)] for final_pair in final_arrange]

    for j in final_pairs:
        for i in range(len(j)):
            x_values = [j[i][0][0], j[i][1][0]]
            y_values = [j[i][0][1], j[i][1][1]]
            plt.arrow(x_values[0], y_values[0], x_values[1]-x_values[0], y_values[1]-y_values[0], head_width = 0.3, width = 0.05, ec='red')
    plt.show()
