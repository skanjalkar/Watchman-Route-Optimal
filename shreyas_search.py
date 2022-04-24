import itertools
import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry.polygon import Polygon
import time
import Shrink_Polygon_Watchman_Route
from matplotlib.animation import FuncAnimation
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


if __name__ == '__main__':
    fig, ax = plt.subplots()
    lim_x, lim_y = 150, 150
    # P = draw_polygon(ax, 8, lim_x,  lim_y)
    # vor = Voronoi(P)
    vor_int, vor_x, vor_y, vor_inside_poly = [], [], [], []
    # for i in vor.vertices:
    #     if i[0] < 0 or i[1] < 0 or i[0] > lim_x or i[1] > lim_y:
    #         continue
    #     x = int(i[0])
    #     y = int(i[1])
    #     vor_int.append((x, y))

    guards, old_p = Shrink_Polygon_Watchman_Route.shrink()
    mx, my = get_min(old_p)
    Guards = translate_poly(guards, mx, my)
    P = translate_poly(old_p, mx, my)
    vor_int = tuple(tuple(map(int, tup)) for tup in Guards)
    polygon = Polygon(P)
    print(polygon)
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
    back_track_track = []
    s = time.time()
    print(f'these are the {vor_inside_poly}')
    for i in range(len(vor_inside_poly)):
        track[vor_inside_poly[i]] = {}
        f = []
        zz = []
        for j in range(len(vor_inside_poly)):
            # path, path_length = search(grid, vor_inside_poly[i], vor_inside_poly[j], polygon)
            path, path_length = astar(grid, vor_inside_poly[i], vor_inside_poly[j])
            f.append(path_length)
            # track[vor_inside_poly[i]][vor_inside_poly[j]] = path_length
            back_track.append(path)
            zz.append(path)
        back_track_track.append(zz)
        paths.append(f)
    # print(track)
    # print(back_track)
    # print(vor_inside_poly)
    v = [False for i in range(len(vor_inside_poly))]
    v[0] = True
    print(paths)
    # print(paths)
    # TSP(paths, v, 0, len(vor_inside_poly), 1, 0)
    print(brute_force(paths))
    # plotTSP(paths, vor_inside_poly, len(paths))
    # print(answer)
    # print(min(answer), answer.index(min(answer)))
    pl, final_path, t = held_karp(paths)
    print(final_path)
    # print(track)
    # print()
    # print(back_track)
    # print(time.time()-s)
    # point_list = [i for i in combinations(vor_inside_poly, 2)]
    # print(point_list)
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
    plt.scatter(vor_x, vor_y, c='blue')
    plt.scatter(vor_x, vor_y, c="red")
    plt.show()


    x, y = zip(*P)
    plt.figure()
    plt.plot(x, y)
    plt.scatter(vor_x, vor_y, c='blue')
    # print(back_track_track)
    final_arrange = []
    final_path.append(0)
    for i in range(len(final_path)-1):
        # print(back_track_track[final_path[i]][final_path[i+1]])
        final_arrange.append(back_track_track[final_path[i]][final_path[i+1]])
    final_pairs = [[j for j in itertools.pairwise(final_pair)] for final_pair in final_arrange]
    print(final_pairs)
    for j in final_pairs:
        for i in range(len(j)):
            x_values = [j[i][0][0], j[i][1][0]]
            y_values = [j[i][0][1], j[i][1][1]]
            # plt.plot(x_values, y_values, c='red')
            plt.arrow(x_values[0], y_values[0], x_values[1]-x_values[0], y_values[1]-y_values[0], head_width = 0.3, width = 0.05, ec='red')
            # plt.annotate("Final path", ())
    plt.show()
