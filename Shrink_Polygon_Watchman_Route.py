import numpy as np
import matplotlib.pyplot as plt
import Shrink_Polygon_AGP

SP = Shrink_Polygon_AGP

P = [(24970,19250),(23600,19250),(20740,22110),(22790,24160),(19395,27554),\
    (17345,25504),(15560,27289),(15560,30215),(11165,30215),(11165,27915),\
    (12435,27915),(15220,24415),(12445,21630),(16865,17210),(19650,19995),\
    (23600,16045),(24970,16045)]
# P = P.reverse   #already clockwise
SP_Points = SP.shrink(P)
SP_Points.append(SP_Points[0])
Pb = SP.create_point_pair(P)

