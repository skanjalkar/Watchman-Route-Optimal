import itertools
import Shrink_Polygon_AGP
import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry.polygon import Polygon
import time
from matplotlib.animation import FuncAnimation
from astar import *
from dijkstra import *
from TSP_search import *


if __name__ == '__main__':
    fig, ax = plt.subplots()
    lim_x, lim_y = 100, 100
    # poly = draw_polygon(ax, 7, lim_x,  lim_y)


    SP = Shrink_Polygon_AGP
    P = [(4, 4), (8, 4), (8, 0), (14, -5), (20, 0), (20, 6), (15, 6), (15, 10), \
         (20, 10), (20, 14), (16, 14), (16, 16), (10, 16), (10, 14), (6, 14), (6, 16), (2, 16) \
        , (2, 14), (0, 14), (-5, 7), (0, 0), (2, -2), (4, 0)]
    row,col = zip(*P)
    row = list(row)
    col = list(col)
    row.append(row[0])
    col.append(row[0])
    ax.plot(row, col)
    P.reverse()  # already clockwise
    P.append(P[0])
    points = SP.shrink(P)
    points.append(points[0])
    Pb = SP.create_point_pair(P)
    Yx = SP.non_intersecting_diag(points, P, Pb)
    Yn = SP.mini_chk_pts(Pb, points, P, Yx)
    Final_Diagonals = SP.clean_up_final(Yn)
    guards = SP.Guards(Final_Diagonals)
    Guards = tuple(tuple(map(int, tup)) for tup in guards)
    print(Guards)
    polygon = Polygon(P)
    grid = grid(lim_x, lim_y, polygon)
    # for i in vor_int:
    #     if polygon.contains(Point(i)):
    #         vor_x.append(int(i[0]))
    #         vor_y.append(int(i[1]))
    #         vor_inside_poly.append(i)
    # print(len(vor_inside_poly))
    track = {}
    back_track = []
    paths = []
    back_track_track = []
    s = time.time()
    for i in range(len(Guards)):
        track[Guards[i]] = {}
        f = []
        zz = []
        for j in range(len(Guards)):
            # path, path_length = search(grid, vor_inside_poly[i], vor_inside_poly[j], polygon)
            path, path_length = astar(grid, Guards[(i)], Guards[(j)], polygon)
            f.append(path_length)
            # track[vor_inside_poly[i]][vor_inside_poly[j]] = path_length
            back_track.append(path)
            zz.append(path)
        back_track_track.append(zz)
        paths.append(f)
    # print(track)
    # print(back_track)
    # print(vor_inside_poly)
    v = [False for i in range(len(Guards))]
    v[0] = True
    print(paths)
    # print(paths)
    TSP(paths, v, 0, len(Guards), 1, 0)
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
    # plt.scatter(vor_x, vor_y, c="red")
    plt.show()

    # plt.scatter(vor_x, vor_y, c='blue')
    print(back_track_track)
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
            plt.arrow(x_values[0], y_values[0], x_values[1]-x_values[0], y_values[1]-y_values[0], head_width = 0.2, width = 0.05, ec='red')
            # plt.annotate("Final path", ())
    plt.show()

