from typing import Final
import numpy as np
import matplotlib.pyplot as plt
import Shrink_Polygon_AGP

SP = Shrink_Polygon_AGP

# P = [(24970,19250),(23600,19250),(20740,22110),(22790,24160),(19395,27554),\
#     (17345,25504),(15560,27289),(15560,30215),(11165,30215),(11165,27915),\
#     (12435,27915),(15220,24415),(12445,21630),(16865,17210),(19650,19995),\
#     # (23600,16045),(24970,16045)]
P = [(4,4),(8,4),(8,0),(14,-5),(20,0),(20,6),(15,6),(15,10),\
    (20,10),(20,14),(16,14),(16,16),(10,16),(10,14),(6,14),(6,16),(2,16)\
       ,(2,14),(0,14),(-5,7),(0,0),(2,-2),(4,0)]
P.reverse()   #already clockwise
P.append(P[0])
points = SP.shrink(P)
points.append(points[0])
Pb = SP.create_point_pair(P)
Yx = SP.non_intersecting_diag(points,P,Pb)
Yn = SP.mini_chk_pts(Pb,points,P,Yx)
Final_Diagonals = SP.clean_up_final(Yn)
Guards = SP.Guards(Final_Diagonals)
print(f'These are the guards {Guards}')


