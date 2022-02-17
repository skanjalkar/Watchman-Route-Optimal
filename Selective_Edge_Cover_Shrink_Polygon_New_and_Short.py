import matplotlib.pyplot as plt
import math
import pyclipper
from shapely.geometry import Point,Polygon #used to chk pt in or out of poly
from time import process_time
T_start = process_time()
X = [];Y = [];Pi = [];PS = [];Xn = []; S = [];Yx = [];Yn = []; Yy = []; Pout = []
MP = []; Ym = [];Yp = [];Poly = []; YN = [];m = []
'Poly = [(24970, 16045), (23600, 16045), (19650, 19995), (16865, 17210), \
        (12445, 21630), (15220, 24415), (12435, 27915), (11165, 27915), \
        (11165, 30215), (15560, 30215), (15560, 27289), (17345, 25504), \
        (19395, 27554), (22790, 24160), (20740, 22110), (23600, 19250), \
        (24970, 19250)]' #clockwise 
Poly = [(24970, 24260), (28385, 24260), (28385, 26060), (34330, 26060), \
        (34330, 23430), (37480, 23430), (37480, 17050), (34330, 17050), \
        (34330, 19250), (31430, 19250), (31430, 12850), (34330, 12850), \
        (34330, 15050), (36510, 15050), (36510, 9720), (32000, 9720), \
        (32000, 7085), (26255, 7085), (26255, 9250), (19590, 9250), \
        (19590, 4800), (11430, 4800), (11430, 7085), (9085, 9430), \
        (12691, 13036), (14498, 11230), (16407, 13139), (14600, 14946), \
        (16865, 17210), (12445, 21630), (11925, 21110), (0, 21110), \
        (0, 28615), (2085, 28615), (2085, 26570), (4785, 26570), (4785, 29085),\
        (6245, 29085), (6245, 34300), (3545, 34300), (3545, 31300), \
        (1300, 31300), (1300, 39400), (4315, 39400), (4315, 41415), \
        (10060, 41415), (10060, 39400), (16740, 39400), (16740, 41416), \
        (25785, 41415), (25785, 31150), (23370, 31150), (23370, 33700), \
        (20670, 33700), (20670, 31500), (16490, 31500), (16490, 30215), \
        (15560, 30215), (15560, 27289), (17345, 25504), (19395, 27554),\
        (22790, 24160), (20740, 22110), (23600, 19250), (24970, 19250)] #clockwise
#Poly = [(7, 7), (6, 5), (8, 2), (5, 0), (0, 0), (4, 10), (0, 13), (3, 13),\
#       (5, 16), (6, 12), (9, 12), (10, 10), (11, 8), (9, 6), (8, 8)] #clockwise
Poly =[(2, 6), (-2, 7.5), (0.5, 11), (4, 11), (2, 8.5), (5, 5), (4, 8), (9, 11),\
        (14, 8), (8, 8), (12, 3),(13, 7), (15, 8), (16, 2), (14, 0), (8, 0), \
        (5, 2), (3, 0), (-2, 3),(-3, 6)] #clockwise
#Poly = [(5,5),(5,7),(7,6),(9,9),(8,2),(7,4),(6,0),(4,3),(1,8)]
#Enter the edges in Req_Pb in clockwise manner, considering the order of poly

'Req_Pb = [[(24970,16045),(23600,16045)],\
          [(16865,17210),(12445,21630)],[(12445,21630),(15220,24415)],\
          [(11165, 27915),(11165, 30215)],[(15560,30215),(15560,27289)],\
          [(17345, 25504),(19395, 27554)]]'
Req_Pb = [[(24970,19250),(24970,24260)],[(24970, 24260), (28385, 24260)],\
          [(28385,24260),(28385,26060)],[(28385, 26060), (34330, 26060)],\
          [(34330,26060),(34330,23430)],[(34330, 23430), (37480, 23430)],\
          [(37480,23430),(37480,17050)],\
          [(19590,9250),(19590,4800)],[(19590, 4800), (11430, 4800)],\
          [(11430, 4800), (11430, 7085)],[(11430, 7085), (9085, 9430)],\
          [(9085,9430),(12691,13036)],[(3545, 34300), (3545, 31300)],\
          [(3545, 31300),(1300, 31300)],[(1300, 31300), (1300, 39400)],\
          [(16740, 39400), (16740, 41416)],[(25785, 41415), (25785, 31150)],\
          [(25785, 31150), (23370, 31150)]]
Req_Pb = [[(12, 3),(13, 7)],[(13, 7), (15, 8)],[(15, 8), (16, 2)],\
         [(16, 2), (14, 0)],[(14, 0), (8, 0)],[(0.5, 11), (4, 11)],\
          [(4, 11),(2, 8.5)],[(15, 8), (16, 2)],[(16, 2), (14, 0)],\
          [(14, 0), (8, 0)]]
#Req_Pb = [[(8, 2), (5, 0)],[(5, 0), (0, 0)],\
#          [(10, 10), (11, 8)],[(11, 8), (9, 6)],[(9, 6), (8, 8)]]
#Req_Pb = [[(5,5),(5,7)],[(9,9),(8,2)],[(4,3),(1,8)]]
print("The Req_Pb",Req_Pb)

'''Poly = list();Xp = list();Yp = list()  
while True:
    Vx = input("Enter the x coordinates:")
    if Vx == "done": break
    try: Vx = float(Vx)
    except: print("Invalid Input");continue
    Xp.append(Vx)
    Vy = input("Enter the y coordinates:")
    if Vy == "done": break
    try: Vy = float(Vy)
    except: print("Invalid Input");continue
    Yp.append(Vy)
Poly = [(Xp[i],Yp[i]) for i in range(0,len(Xp))]'''
#Poly.reverse()
P = Poly
AP = P
P.append(P[0])
def det(a, b): #readymade function taken from the net
        return a[0] * b[1] - a[1] * b[0]
def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1]) #Typo was here
    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')
    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y
def shrink(Poly):
# how much the coordinates are moved as an absolute value
    shrink_x = 0.005
    shrink_y = 0.005
    # coords must be clockwise
    lines = [[Poly[i-1], Poly[i]] for i in range(len(Poly))]
    new_lines = []
    for i in lines:
        dx = i[1][0] - i[0][0]
        dy = i[1][1] - i[0][1]
    # this is to take into account slopes
        if (dx*dx + dy*dy)==0:
            continue
        else:
            factor = 1 / (dx*dx + dy*dy)**0.5
            new_dx = dy*shrink_x * factor
            new_dy = dx*shrink_y * factor
            new_lines.append([(i[0][0] + new_dx, i[0][1] - new_dy),
                            (i[1][0] + new_dx, i[1][1] - new_dy)])
# find position of intersection of all the lines
    new_polygon = []
    for i in range(len(new_lines)):
        new_polygon.append((line_intersection(new_lines[i-1], new_lines[i])))
    return new_polygon
Pc = shrink(Poly)
Pc.append(Pc[0])
#print("The Pc is:",Pc)
AAP = Pc
#print("Bhai P kya Hai",P)

def Sorting(lst): 
    lst2 = sorted(lst, key=len, reverse = True) 
    return lst2 
''' orientation function: To check the orientation on points (x1,y1),(x2,y2),(x3,y3)'''
def orientation(x1,y1,x2,y2,x3,y3):
        val = (float((y2-y1)*(x3-x2)))-(float((x2-x1)*(y3-y2)))
        if (val>0):
            return 1 #clockwise
        elif (val<0):
            return 2 #counterclockwise
        else:
            return 0 #collinear
''' point_in_seg_area function: To check if the point lies in segment area'''
def point_in_seg_area(x1,y1,x2,y2,x3,y3):  
        if ((x2<=max(x1,x3)) and (x2>=min(x1,x3))\
                and (y2<=max(y1,y3)) and (y2>=min(y1,y3))):
            return True
        return False
''' check_intersection function: To check if the line formed by points (x1,y1) and (x2,y2) intersects line
      formed by (x3,y3) and (x4,y4)'''
def check_intersection(x1,y1,x2,y2,x3,y3,x4,y4):
        o1 = orientation(x1,y1,x2,y2,x3,y3)
        o2 = orientation(x1,y1,x2,y2,x4,y4)
        o3 = orientation(x3,y3,x4,y4,x1,y1)
        o4 = orientation(x3,y3,x4,y4,x2,y2)
        if ((o1 == 0) and point_in_seg_area(x1,y1,x3,y3,x2,y2)): #both are neede to tell if the point is on the segment 
            return False
        if ((o2 == 0) and point_in_seg_area(x1,y1,x4,y4,x2,y2)):
            return False
        if ((o3 == 0) and point_in_seg_area(x3,y3,x1,y1,x4,y4)):
            return False
        if ((o4 == 0) and point_in_seg_area(x3,y3,x1,y1,x4,y4)):
            return False
        if ((o1!=o2) and (o3!=o4)):
            return True
        return  False
def create_point_pair(P):
    Pb = []
    for i in range(len(P)-1):
        Pa = []
        Pa.append(P[i])
        Pa.append(P[i+1])
        Pb.append(Pa)
    return Pb

Pb = create_point_pair(P)
#print("The pair of points are:",create_point_pair(P))
Yx = []

def non_intersecting_diag(Pc,P, Pb):
    for i in range(len(Pc)-1):
        S = []
        for j in range(len(P)-1):
            Pi = []
            Pi.append(Pc[i])
            Pi.append(P[j])
            S.append(Pi)
        PS.append(S)
    #print("Bhai PS:",PS)
    for n in range(len(PS)):
        for k in range(len(PS[n])):
            Xn = []
            for l in range(len(P)-1):
                if  check_intersection(PS[n][k][0][0],PS[n][k][0][1],PS[n][k][1][0]\
                    ,PS[n][k][1][1],P[l][0],P[l][1],P[l+1][0],P[l+1][1])==True:
                      continue
                else:
                      Xn.append(PS[n][k][0])
                      Xn.append(PS[n][k][1])
            Y = []
            if len(Xn) == 2*(len(P)-1): #no intersection with any polygon side
               Y.append(Xn[0])
               Y.append(Xn[1])
            if Y == []:
                continue
            else:
                Yx.append(Y)
    #print("Yx is:",Yx)
    for m in range(len(Yx)):
        px = float((Yx[m][0][0]+Yx[m][1][0])/2)
        py = float((Yx[m][0][1]+Yx[m][1][1])/2)
        mp = (px,py)
        if not (Point(mp).within(Polygon(AP))): #chk point in or out
               Pout.append(Yx[m])
        MP.append(mp)
    #print("The list of outer lines:",Pout)
    for n in range(len(Pout)):
            if Pout[n] in Yx:
                Yx.remove(Pout[n])
    return Yx
Yx = non_intersecting_diag(Pc,P,Pb)


def mini_chk_pts(Req_Pb,Pc,P,Yx):
    Yn=[];M=[];Ys1=[];Yk1=[];Yy1=[];Yf1 = [];Ye1 = []; R = []
    for r in range(len(Pc)-1):#this is important for arranging the diagonals.
        Yy1 = []
        for s in range(len(Yx)):
            if Pc[r] == Yx[s][0]:
               Yy1.append(Yx[s])
        if not Yy1 == []:
               Yy1.append(Yy1[0])
        Ys1.append(Yy1)
    Yk1 = Sorting(Ys1)  #sorting in descending order of  length of sub-list.
    #print("The list Yk1 is:",Yk1)
    for b in range(len(Yk1)):
            Yg = []
            for c in range(len(Yk1[b])-1):
                 for a in range(len(Req_Pb)): #look for alternate*******
                     Yf = []
                     if ((Req_Pb[a][0] == Yk1[b][c][1]) and (Req_Pb[a][1] == Yk1[b][c+1][1])):
                         Yf.append(Yk1[b][c])
                         Yf.append(Yk1[b][c+1])
                         Yg.append(Yf)
            if not Yg == []:
                Ye1.append(Yg)
    Yf1 = Sorting(Ye1)
    F = Req_Pb
    Yf2 = []
    while F != []:
        Yy = []; Ys = []; M = []
        for a in range(len(Yf1)):
            Yy = []
            for b in range(len(Yf1[a])):
                for c in range(len(F)):
                    if (F[c][0] in Yf1[a][b][0]) and (F[c][1] in Yf1[a][b][1])\
                       and (Yf1[a][b][0][1] in F[c]) and Yf1[a][b][1][1] in F[c]:
                        Yy.append(Yf1[a][b])
            Ys.append(Yy)
        Yf2 = Sorting(Ys)
        A2 = Yf2[0]
        for i in range(len(Yf2[0])):
            Yn.append(Yf2[0][i])
        Yf2.remove(Yf2[0])
        for j in range(len(F)):
            for k in range(len(A2)):
                if (F[j][0] == A2[k][0][1]) and (F[j][1] == A2[k][1][1]):
                        M.append(F[j])
                else:
                        continue
        F2 = []
        for l in range(len(F)):
                if not F[l] in M:
                        F2.append(F[l])
                else:
                        continue
        Yf1 = Yf2
        F = F2
    return Yn

Pfinal = mini_chk_pts(Req_Pb,Pc,P,Yx)
#print("The Pfinal is:",Pfinal)
final = []
R = []
for i in Pfinal:
    if not i in final:
       final.append(i)
r = []
for p in range(len(final)): #solution for adjecent points
    for q in range(len(final)): 
        for r in range(len(Pc)-1):
            if (final[p][0][0] or final[p][1][0]) == Pc[r]:
                if (Pc[r+1] or Pc[r-1])==(final[q][0][1] or final[q][1][1]):
                    R.append(final[q])
#print("Bro R is:",R)
for r in range(len(R)):
    if R[r] in final:
       final.remove(R[r])   
Yn = final
print("The co-ordinates are:",Yn)

def plt_plot(P,Yn,Req_Pb):
    Rx = [];Ry = []
    Px = [];Py = [];Dx = [];Dy = [];Sx = [];Sy = [];APx = [];APy = []
    for h in range(len(AAP)):
        APx.append(AAP[h][0])
        APy.append(AAP[h][1])
    for j in range(len(Yn)):    
        Dx=[];Dy=[]
        Dx.append(Yn[j][0][0][0])
        Dy.append(Yn[j][0][0][1])
        Dx.append(Yn[j][0][1][0])
        Dy.append(Yn[j][0][1][1])
        Dx.append(Yn[j][1][0][0])
        Dy.append(Yn[j][1][0][1])
        Dx.append(Yn[j][1][1][0])
        Dy.append(Yn[j][1][1][1])
        Sx.append(Yn[j][0][0][0])
        Sy.append(Yn[j][0][0][1])
        plt.plot(Dx,Dy, color = 'g')
    for i in range(len(P)):
        Px.append(P[i][0])
        Py.append(P[i][1])
    plt.plot(Px,Py,color = 'b')
    for i in range(len(Req_Pb)):
        Rx = [];Ry = []
        Rx.append(Req_Pb[i][0][0])
        Ry.append(Req_Pb[i][0][1])
        Rx.append(Req_Pb[i][1][0])
        Ry.append(Req_Pb[i][1][1])
        plt.plot(Rx,Ry,color = 'r')
   # plt.plot(APx,APy,color = 'r')
    plt.scatter(Sx,Sy,s = 700,marker = '.',color = 'k')
    return plt.show()
plt_plot(P,Yn,Req_Pb)
T_stop = process_time()
print("Final Elasped time:", T_stop, T_start)
print("Final Elasped time during the whole program in seconds:",((T_stop) - (T_start)))



'''for m in range(len(Yx)):
        px = float((Yx[m][0][0]+Yx[m][1][0])/2)
        py = float((Yx[m][0][1]+Yx[m][1][1])/2)
        mp = (px,py)
        px1 = float(((Yx[m][0][0])*(2/3) + (Yx[m][1][0])*(1/3)))
        py1 = float(((Yx[m][0][1])*(2/3) + (Yx[m][1][1])*(1/3)))
        mp1 = (px1,py1)
        px2 = float(((Yx[m][0][0])*(1/3) + (Yx[m][1][0])*(2/3)))
        py2 = float(((Yx[m][0][1])*(1/3) + (Yx[m][1][1])*(2/3)))
        mp2 = (px2,py2)
        px3 = float(((Yx[m][0][0])*(1/4) + (Yx[m][1][0])*(3/4)))
        py3 = float(((Yx[m][0][1])*(1/4) + (Yx[m][1][1])*(3/4)))
        mp3 = (px3,py3)
        px4 = float(((Yx[m][0][0])*(1/5) + (Yx[m][1][0])*(4/5)))
        py4 = float(((Yx[m][0][1])*(1/5) + (Yx[m][1][1])*(4/5)))
        mp4 = (px4,py4)
        px5 = float(((Yx[m][0][0])*(1/6) + (Yx[m][1][0])*(5/6)))
        py5 = float(((Yx[m][0][1])*(1/6) + (Yx[m][1][1])*(5/6)))
        mp5 = (px5,py5)
        if not (Point(mp).within(Polygon(AP)) or Point(mp1).within(Polygon(AP))\
                or Point(mp2).within(Polygon(AP)) or Point(mp3).within(Polygon(AP))\
                or Point(mp4).within(Polygon(AP)) or Point(mp5).within(Polygon(AP))): #chk point in or out
               Pout.append(Yx[m])
        MP.append(mp)'''



''' for b in range(len(Yk1)):
            Yg = []
            for c in range(len(Yk1[b])-1):
                 for a in range(len(Req_P)-1): #look for alternate*******
                     Yf = []
                     if ((Req_P[a] == Yk1[b][c][1]) and (Req_P[a+1] == Yk1[b][c+1][1])):
                         Yf.append(Yk1[b][c])
                         Yf.append(Yk1[b][c+1])
                         Yg.append(Yf)
            if not Yg == []:
                Ye1.append(Yg)
    Yf1 = Sorting(Ye1)'''
