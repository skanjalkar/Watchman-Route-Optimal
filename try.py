# from shapely import geometry
# import matplotlib.pyplot as plt

# # your variables
# coords = [(0, 0), (0, 100), (20, 100), (30, 60), (40, 100), (60, 100), (60, 0), (40, 10), (40, 40), (20, 40), (20, 10)]
# lines = [[coords[i-1], coords[i]] for i in range(len(coords))]

# # your factor of 10%
# # Note: with 20% the polygon becomes a multi-polygon, so a loop for plotting would be needed.
# factor = 0.1

# # code from nathan
# xs = [i[0] for i in coords]
# ys = [i[1] for i in coords]
# x_center = 0.5 * min(xs) + 0.5 * max(xs)
# y_center = 0.5 * min(ys) + 0.5 * max(ys)

# min_corner = geometry.Point(min(xs), min(ys))
# max_corner = geometry.Point(max(xs), max(ys))
# center = geometry.Point(x_center, y_center)
# shrink_distance = center.distance(min_corner)*factor

# assert abs(shrink_distance - center.distance(max_corner)) < 0.0001

# my_polygon = geometry.Polygon(coords)
# my_polygon_shrunken = my_polygon.buffer(-shrink_distance)


# x, y = my_polygon.exterior.xy
# plt.plot(x,y)
# x, y = my_polygon_shrunken.exterior.xy
# plt.plot(x,y)

# # to net let the image be distorted along the axis
# plt.axis('equal')

# plt.show()


A = [1,11,11,5,4,6,3,7]
print(A.index(max(A)))