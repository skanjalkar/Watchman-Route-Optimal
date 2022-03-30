import numpy as np
from scipy.spatial import  Voronoi, voronoi_plot_2d
# import sys
# sys.path.append('/Educational/NTU_Research_Work/A_Research_Internship/Python_Codes/')

import Final_Voronoi_New_and_Short

Poly = [(24970,19250),(23600,19250),(20740,22110),(22790,24160),(19395,27554),\
     (17345,25504),(15560,27289),(15560,30215),(11165,30215),(11165,27915),\
     (12435,27915),(15220,24415),(12445,21630),(16865,17210),(19650,19995),\
     (23600,16045),(24970,16045)]

'''Poly = [(24970,19250),(23600,19250),(20740,22110),(22790,24160),(19395,27554)\
    ,(17345,25504),(15560,27289),(15560,30215),(16490,30215),(16490,31500)\
    ,(20670,31500),(20670,33700),(23370,33700),(23370,31150),(25785,31150)\
    ,(25785,41415),(16740,41416),(16740,39400),(10060,39400),(10060,41415)\
    ,(4315,41415),(4315,39400),(1300,39400),(1300,31300),(3545,31300)\
    ,(3545,34300),(6245,34300),(6245,29085),(4785,29085),(4785,26570)\
    ,(2085,26570),(2085,28615),(0,28615),(0,21110),(11925,21110),(12445,21630)\
    ,(16865,17210),(14600,14946),(16407,13139),(14498,11230),(12691,13036)\
    ,(9085,9430),(11430,7085),(11430,4800),(19590,4800),(19590,9250),(26255,9250)\
    ,(26255,7085),(32000,7085),(32000,9720),(36510,9720),(36510,15050)\
    ,(34330,15050),(34330,12850),(31430,12850),(31430,19250),(34330,19250)\
    ,(34330,17050),(37480,17050),(37480,23430),(34330,23430),(34330,26060)\
    ,(28385,26060),(28385,24260),(24970,24260)]'''

points = np.array(Poly)
vor = Voronoi(points)
VG = Final_Voronoi_New_and_Short
regions, vertices = VG.voronoi_finite_polygons_2d(vor)
fig = voronoi_plot_2d(vor)
vert = []

for i in range(len(vertices)):
    Vx = [] ;Vy = []
    Vx = (vertices[i][0])
    Vy = (vertices[i][1])
    V = (Vx,Vy)
    vert.append(V)

P = Poly; AP = P; P.append(P[0]); Pc = vert; Pc.append(Pc[0]); AAP = Pc
print(P)

Pb = VG.create_point_pair(P)
Yx = VG.non_intersecting_diag(Pc,P,Pb)
Yn = VG.mini_chk_pts(Pb,Pc,P,Yx)

Pfinal = Yn; final = []; R = [];r = []
for i in Pfinal:
    if not i in final:
       final.append(i)
for p in range(len(final)): #solution for adjecent points
    for q in range(len(final)): #this is a big change!!!!!!!!!
        for r in range(len(Pc)-1):
            if (final[p][0][0] or final[p][1][0]) == Pc[r]:
                if (Pc[r+1] or Pc[r-1])==(final[q][0][1] or final[q][1][1]):
                    R.append(final[q])
for r in range(len(R)):#PREVENT REPITITION
    if R[r] in final:
       final.remove(R[r])
Yn = final

print(Yn)
