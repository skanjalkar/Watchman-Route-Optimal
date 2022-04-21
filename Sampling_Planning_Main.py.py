# Main and helper function

from re import A
from PIL import Image
import numpy as np
from RRT import RRT

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

if __name__ == "__main__":
    # Load the map
    start = (200, 75)
    goal  = (30, 250)
    map_array = load_map("D:\Educational\A WPI Assignments and Materials\Motion Planning\Assignments\Assignment 5\Advanced-Search-Algorithms\informed_RRT\WPI_map.jpg", 0.3)

    # start = (40.0, 139.59183673469389)
    # goal = (100.0, 58.96668365142222)
    # map_array = load_map("D:\Educational\NTU_Research_Work\A_Research_Internship\Blank.png",0.3)
    
    # Planning class
    RRT_planner = RRT(map_array, start, goal)

    # Search with RRT and RRT*
    # RRT_planner.RRT(n_pts=1000)
    # RRT_planner.RRT_star(n_pts=2000)
    RRT_planner.informed_RRT_star(n_pts=2000)
