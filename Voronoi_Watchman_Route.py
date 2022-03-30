import numpy as np
from scipy.spatial import  Voronoi, voronoi_plot_2d
import sys
sys.path.append('/Educational/NTU_Research_Work/A_Research_Internship/Python_Codes/')

from AGP_Codes import Final_Voronoi_New_and_Short

Poly = [(24970,19250),(23600,19250),(20740,22110),(22790,24160),(19395,27554),\
     (17345,25504),(15560,27289),(15560,30215),(11165,30215),(11165,27915),\
     (12435,27915),(15220,24415),(12445,21630),(16865,17210),(19650,19995),\
     (23600,16045),(24970,16045)]

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

P = Poly; AP = P; P.append(P[0])
Pc = vert; Pc.append(Pc[0]); AAP = Pc
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