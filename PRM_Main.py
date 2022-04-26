# Main and helper function

from PIL import Image
import numpy as np
from RRT import RRT
from PRM import PRM
import time
import matplotlib.pyplot as plt


def load_map(file_path, resolution_scale):
    ''' Load map from an image and return a 2D binary numpy array
        where 0 represents obstacles and 1 represents free space
    '''
    # Load the image with grayscale
    img = Image.open(file_path).convert('L')
    # Rescale the image
    size_x, size_y = img.size
    new_x, new_y  = int(size_x*resolution_scale), int(size_y*resolution_scale)
    img = img.resize((new_x, new_y), Image.ANTIALIAS)

    map_array = np.asarray(img, dtype='uint8')

    # Get bianry image
    threshold = 127
    map_array = 1 * (map_array > threshold)

    # Result 2D numpy array
    return map_array

Start = time.time() #starting the time
print("The start time is:",Start)
if __name__ == "__main__":
    # Load the map
    # start = (200, 75)
    # goal  = (30, 250)
    # map_array = load_map("D:\Educational\A WPI Assignments and Materials\Motion Planning\Assignments\Assignment 3\Standard Search Algorithms\WPI_map.jpg", 0.3)

    # # Planning class
    # PRM_planner = PRM(map_array)
    # RRT_planner = RRT(map_array, start, goal)

    # Search with PRM
    points = [(63,159),(121,220),(273,324),(63,159)]
    points = [(147,130),(228,90),(265,322),(125,265),(147,130)] # Polygon with holes
    print(points)
    for i in range(len(points)-1):
        start = points[i]
        goal = points[i+1]
        map_array = load_map("D:\Educational\A WPI Assignments and Materials\Motion Planning\Project\Colored Polygons\PH_BW1.png",1)
        # RRT_planner = RRT(map_array, start, goal)
        PRM_planner = PRM(map_array)
        # PRM_planner.sample(n_pts=1000, sampling_method="uniform")
        # PRM_planner.search(start, goal)
        # PRM_planner.sample(n_pts=1000, sampling_method="random")
        # PRM_planner.search(start, goal)
        # PRM_planner.sample(n_pts=5000, sampling_method="gaussian")
        # PRM_planner.search(start, goal)
        PRM_planner.sample(n_pts=20000, sampling_method="bridge")
        PRM_planner.search(start, goal)
    End = time.time()
    print("The End time is:",End),print("The runtime is:",(End-Start))
    plt.show()
    # Search with RRT and RRT*
    # RRT_planner.RRT(n_pts=1000)
    # RRT_planner.RRT_star(n_pts=2000)
