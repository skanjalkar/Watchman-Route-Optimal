from typing import Final
import numpy as np
import matplotlib.pyplot as plt
import Shrink_Polygon_AGP

SP = Shrink_Polygon_AGP

#Keep the polygon same in both the codes (this code, and shrink_polygon_AGP.py)

# P = [(24970,19250),(23600,19250),(20740,22110),(22790,24160),(19395,27554),\
#     (17345,25504),(15560,27289),(15560,30215),(11165,30215),(11165,27915),\
#     (12435,27915),(15220,24415),(12445,21630),(16865,17210),(19650,19995),\
#     (23600,16045),(24970,16045)]
P = [(20,140),(0,140),(-50,70),(0,0),(20,-20),(40,0),(40,40),\
    (40,40),(80,40),(80,0),(140,-50),(200,0),(200,60),(150,60),(150,100),\
    (200,100),(200,140),(160,140),(160,160),(100,160),(100,140),(60,140),(60,160),(20,160)] 

P = [(12,3),(8,8),(14,8),(9,11),(4,8),\
       (6,5),(2,8.5),(4,11),(0.5,11),(-2,7.5),(2,6),(-3,6),\
       (-2,3),(3,0),(4,2),(8,0),(14,0),(16,2),(15,8),(13,7)]

P.reverse()   #already clockwise
P.append(P[0])
points = SP.shrink(P)
points.append(points[0])
points = SP.sampling_points(points)
AAP = points
Pb = SP.create_point_pair(P)
Yx = SP.non_intersecting_diag(points,P,Pb)
Yn = SP.mini_chk_pts(Pb,points,P,Yx)
Final_Diagonals = SP.clean_up_final(Yn)
Guards = SP.Guards(Final_Diagonals)
print(Guards)
''' Below is the code for Watchmen Route'''
SP.plt_plot(P,Yn)

''' Add the watchman code below'''
# [(40.0, 139.59183673469389), (100.0, 58.96668365142222)]