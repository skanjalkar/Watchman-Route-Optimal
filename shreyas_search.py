import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry.polygon import Polygon
import time
from matplotlib.animation import FuncAnimation
from astar import *
from dijkstra import *
from TSP_search import *


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


if __name__ == '__main__':
    fig, ax = plt.subplots()
    lim_x, lim_y = 100, 100
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
    print(paths)
    # print(paths)
    TSP(paths, v, 0, len(vor_inside_poly), 1, 0)
    print(brute_force(paths))
    # plotTSP(paths, vor_inside_poly, len(paths))
    # print(answer)
    print(min(answer), answer.index(min(answer)))
    pl, final_path = held_karp(paths)
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
    # ani = FuncAnimation(fig, animate(back_track), frames=30, interval=500, repeat=False)
    plt.scatter(vor_x, vor_y, c="red")
    plt.show()

    plt.scatter(vor_x, vor_y, c='blue')
    for i in final_path:
        pass

