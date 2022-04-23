import matplotlib.pyplot as plt
import math
import pyclipper
from shapely.geometry import Point, Polygon  # used to chk pt in or out of poly
import time

from sympy import Max

# Start = time.time() #starting the time
# print("The start time is:",Start)
X = [];
Y = [];
Pi = [];
PS = [];
Xn = [];
S = [];
Yx = [];
Yn = [];
Yy = [];
Pout = []
MP = [];
Ym = [];
Yp = [];
Poly = [];
YN = [];
m = []
Poly = [(24970, 19250), (23600, 19250), (20740, 22110), (22790, 24160), (19395, 27554), \
        (17345, 25504), (15560, 27289), (15560, 30215), (11165, 30215), (11165, 27915), \
        (12435, 27915), (15220, 24415), (12445, 21630), (16865, 17210), (19650, 19995), \
        (23600, 16045), (24970, 16045)]

Poly = [(24970, 19250), (23600, 19250), (20740, 22110), (22790, 24160), (19395, 27554) \
    , (17345, 25504), (15560, 27289), (15560, 30215), (16490, 30215), (16490, 31500) \
    , (20670, 31500), (20670, 33700), (23370, 33700), (23370, 31150), (25785, 31150) \
    , (25785, 41415), (16740, 41416), (16740, 39400), (10060, 39400), (10060, 41415) \
    , (4315, 41415), (4315, 39400), (1300, 39400), (1300, 31300), (3545, 31300) \
    , (3545, 34300), (6245, 34300), (6245, 29085), (4785, 29085), (4785, 26570) \
    , (2085, 26570), (2085, 28615), (0, 28615), (0, 21110), (11925, 21110), (12445, 21630) \
    , (16865, 17210), (14600, 14946), (16407, 13139), (14498, 11230), (12691, 13036) \
    , (9085, 9430), (11430, 7085), (11430, 4800), (19590, 4800), (19590, 9250), (26255, 9250) \
    , (26255, 7085), (32000, 7085), (32000, 9720), (36510, 9720), (36510, 15050) \
    , (34330, 15050), (34330, 12850), (31430, 12850), (31430, 19250), (34330, 19250) \
    , (34330, 17050), (37480, 17050), (37480, 23430), (34330, 23430), (34330, 26060) \
    , (28385, 26060), (28385, 24260), (24970, 24260)]

'''
Poly = [(8,8),(9,6),(11,8),(10,10),(9,12),(6,12),(5,16),(3,13),(0,13),(4,10)\
       ,(0,0),(5,0),(8,2),(6,5),(7,7)]
Poly = [(-3,6),(-2,3),(3,0),(5,2),(8,0),(14,0),(16,2),(15,8),(13,7),(12,3),(8,8),(14,8)\
       ,(9,11),(4,8),(5,5),(2,8.5),(4,11),(0.5,11),(-2,7.5),(2,6)]
Poly = [(8000,8000),(9000,6000),(11000,8000),(9000,12000),(6000,12000)\
        ,(5000,16000),(3000,13000),(0,13000),(4000,10000),\
        (0,0),(5000,0),(8000,2000),(6000,5000),(7000,7000)]
Poly = [(0,0),(10,0),(10,1),(10,5),(9,5),(8,1),(8,5),(7,5)\
        ,(6,1),(6,5),(5,5),(4,1),(4,5),(3,5),(2,1),(2,5),(1,5),(0,1)]
Poly = [(0,0),(4,0),(4,4),(0,4)]
Poly = [(10,10),(8,6),(7,8),(5,6),(4,7),(2,3)\
       ,(0,4),(0,2),(1,1),(6,0),(8,0),(10,2),(9,4)]
Poly = [(0,6),(1,1),(3,0),(7,2),(5,4),(7,5),(6,8),(4,7),(2,11)]
Poly = [(1,1),(2,2),(1,5),(4,3),(5,7),(5,5),(10,10),(10,-5),(6,-1),(3,-3)]
Poly = [(0,0),(2,2),(0,4),(3,4),(3,0)]
Poly = [(1,2),(0,0),(2,-3),(5,-3),(7,1),(6,2),(7,4),(4,6),(3,6),(-3,4)]
Poly = [(0,0),(10,0),(10,1),(10,5),(8,1),(8,5),(6,1),(6,5),(4,1),(4,5),(2,1),(2,5),(0,1)]
Poly = [(1,8),(4,3),(6,0),(7,4),(8,2),(9,9),(7,6),(5,7),(5,5)]
Poly = [(3,16),(4,13),(5,11),(6,10),(8,9),(7,7),(5,4),(6,2),(7,4),(7,0),\
       (4,0),(1,0),(1,1),(2,1),(3,3),(2,3),(2,6),(1,3),(0,3),(1,6)]
'''
'''Poly = [(4000,4000),(8000,4000),(8000,0000),(14000,-5000),(20000,0),(20000,6000),(15000,6000),(15000,10000)\
        ,(20000,10000),(20000,14000),(16000,14000),(16000,16000),(10000,16000),(10000,14000),(6000,14000),(6000,16000),(2000,16000)\
       ,(2000,14000),(0,14000),(-5000,7000),(0,0),(2000,-2000),(4000,0)]'''
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
# Poly = [(4, 4), (8, 4), (8, 0), (14, -5), (20, 0), (20, 6), (15, 6), (15, 10), \
#         (20, 10), (20, 14), (16, 14), (16, 16), (10, 16), (10, 14), (6, 14), (6, 16), (2, 16) \
#     , (2, 14), (0, 14), (-5, 7), (0, 0), (2, -2), (4, 0)]

# Poly = [(2000,14000),(0,14000),(-5000,7000),(0,0),(2000,-2000),(4000,0),(4000,4000),\
#     (4000,4000),(8000,4000),(8000,0),(14000,-5000),(20000,0),(20000,6000),(15000,6000),(15000,10000),\
#         (20000,10000),(20000,14000),(16000,14000),(16000,16000),(10000,16000),(10000,14000),(6000,14000),(6000,16000),(2000,16000)]
# Poly = [(1,8),(4,3),(6,0),(7,4),(8,2),(9,9),(7,6),(5,7),(5,5)]

Poly = [(120,30),(80,80),(140,80),(90,110),(40,80),\
           (60,50),(20,85),(40,110),(5,110),(-20,75),(20,60),(-30,60),\
           (-20,30),(30,0),(40,20),(80,0),(140,0),(160,20),(150,80),(130,70)]


Poly.reverse();
P = Poly;
AP = P;
P.append(P[0])


def det(a, b):  # readymade function taken from the net
    return a[0] * b[1] - a[1] * b[0]


def point_of_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])  # Typo was here
    div = det(xdiff, ydiff)
    if div == 0:
        raise Exception('lines do not intersect')
    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y


def shrink(Poly):
    # how much the coordinates are moved as an absolute value
    shrink_x = 5
    shrink_y = 5
    # coords must be clockwise
    lines = [[Poly[i - 1], Poly[i]] for i in range(len(Poly))]
    new_lines = []
    for i in lines:
        dx = i[1][0] - i[0][0]
        dy = i[1][1] - i[0][1]
        # this is to take into account slopes
        if (dx * dx + dy * dy) == 0:
            continue
        else:
            factor = 1 / (dx * dx + dy * dy) ** 0.5
            new_dx = dy * shrink_x * factor
            new_dy = dx * shrink_y * factor
            new_lines.append([(i[0][0] + new_dx, i[0][1] - new_dy),
                              (i[1][0] + new_dx, i[1][1] - new_dy)])
    # find position of intersection of all the lines
    new_polygon = []
    for i in range(len(new_lines)):
        new_polygon.append((point_of_intersection(new_lines[i - 1], new_lines[i])))
    return new_polygon


Pc = shrink(Poly)
Start = time.time()  # starting the time
# print("The start time is:",Start)
Pc.append(Pc[0])
AAP = Pc


def Sorting(lst):
    lst2 = sorted(lst, key=len, reverse=True)
    return lst2


''' orientation function: To check the orientation on points (x1,y1),(x2,y2),(x3,y3)'''


def orientation(x1, y1, x2, y2, x3, y3):
    val = (float((y2 - y1) * (x3 - x2))) - (float((x2 - x1) * (y3 - y2)))
    if (val > 0):
        return 1  # clockwise
    elif (val < 0):
        return 2  # counterclockwise
    else:
        return 0  # collinear


''' point_in_seg_area function: To check if the point lies in segment area'''


def point_in_seg_area(x1, y1, x2, y2, x3, y3):
    if ((x2 <= max(x1, x3)) and (x2 >= min(x1, x3)) \
            and (y2 <= max(y1, y3)) and (y2 >= min(y1, y3))):
        return True
    return False


''' check_intersection function: To check if the line formed by points (x1,y1) and (x2,y2) intersects line
      formed by (x3,y3) and (x4,y4)'''


def check_intersection(x1, y1, x2, y2, x3, y3, x4, y4):
    o1 = orientation(x1, y1, x2, y2, x3, y3)
    o2 = orientation(x1, y1, x2, y2, x4, y4)
    o3 = orientation(x3, y3, x4, y4, x1, y1)
    o4 = orientation(x3, y3, x4, y4, x2, y2)
    if ((o1 == 0) and point_in_seg_area(x1, y1, x3, y3, x2,
                                        y2)):  # both are neede to tell if the point is on the segment
        return False
    if ((o2 == 0) and point_in_seg_area(x1, y1, x4, y4, x2, y2)):
        return False
    if ((o3 == 0) and point_in_seg_area(x3, y3, x1, y1, x4, y4)):
        return False
    if ((o4 == 0) and point_in_seg_area(x3, y3, x1, y1, x4, y4)):
        return False
    if ((o1 != o2) and (o3 != o4)):
        return True
    return False


def create_point_pair(P):
    Pb = []
    for i in range(len(P) - 1):
        Pa = []
        Pa.append(P[i])
        Pa.append(P[i + 1])
        Pb.append(Pa)
    return Pb


# Code perfect................................................#
Pb = create_point_pair(P)
# print("The pair of points are:",create_point_pair(P))
Yx = []


# Code perfect................................................#
def non_intersecting_diag(Pc, P, Pb):
    for i in range(len(Pc) - 1):
        S = []
        for j in range(len(P) - 1):
            Pi = []
            Pi.append(Pc[i])
            Pi.append(P[j])
            S.append(Pi)
        PS.append(S)
    # print("Bhai PS:",PS)
    for n in range(len(PS)):
        for k in range(len(PS[n])):
            Xn = []
            for l in range(len(P) - 1):
                if check_intersection(PS[n][k][0][0], PS[n][k][0][1], PS[n][k][1][0] \
                        , PS[n][k][1][1], P[l][0], P[l][1], P[l + 1][0], P[l + 1][1]) == True:
                    continue
                else:
                    Xn.append(PS[n][k][0])
                    Xn.append(PS[n][k][1])
            Y = []
            if len(Xn) == 2 * (len(P) - 1):  # no intersection with any polygon side
                Y.append(Xn[0])
                Y.append(Xn[1])
            if Y == []:
                continue
            else:
                Yx.append(Y)
    # print("Yx is:",Yx)
    for m in range(len(Yx)):
        px = float((Yx[m][0][0] + Yx[m][1][0]) / 2)
        py = float((Yx[m][0][1] + Yx[m][1][1]) / 2)
        mp = (px, py)
        if not (Point(mp).within(Polygon(AP))):  # chk point in or out
            Pout.append(Yx[m])
        MP.append(mp)
    # print("The list of outer lines:",Pout)
    for n in range(len(Pout)):
        if Pout[n] in Yx:
            Yx.remove(Pout[n])
    return Yx


# Code perfect................................................#
Yx = non_intersecting_diag(Pc, P, Pb)


def mini_chk_pts(Pb, Pc, P, Yx):
    Yn = [];
    M = [];
    Ys1 = [];
    Yk1 = [];
    Yy1 = [];
    Yf1 = [];
    Ye1 = [];
    R = []
    for r in range(len(Pc) - 1):  # this is important for arranging the diagonals.
        Yy1 = []
        for s in range(len(Yx)):
            if Pc[r] == Yx[s][0]:
                Yy1.append(Yx[s])
        if not Yy1 == []:
            Yy1.append(Yy1[0])
        Ys1.append(Yy1)
    Yk1 = Sorting(Ys1)  # sorting in descending order of  length of sub-list.
    # print("The list Yk1 is:",Yk1)
    for b in range(len(Yk1)):
        Yg = []
        for c in range(len(Yk1[b]) - 1):
            for a in range(len(P) - 1):
                Yf = []
                if ((P[a] == Yk1[b][c][1]) and (P[a + 1] == Yk1[b][c + 1][1])):
                    Yf.append(Yk1[b][c])
                    Yf.append(Yk1[b][c + 1])
                    Yg.append(Yf)
        if not Yg == []:
            Ye1.append(Yg)
    Yf1 = Sorting(Ye1)
    F = Pb
    Yf2 = []

    while F != []:
        Yy = [];
        Ys = [];
        M = []

        for a in range(
                len(Yf1)):  # So this loop is find the guards which guard maximum of the unguarded edges (first case is similar to yf1)
            Yy = []
            for b in range(len(Yf1[a])):
                for c in range(len(F)):
                    if (F[c][0] in Yf1[a][b][0]) and (F[c][1] in Yf1[a][b][1]) \
                            and (Yf1[a][b][0][1] in F[c]) and Yf1[a][b][1][1] in F[c]:
                        Yy.append(Yf1[a][b])
            Ys.append(Yy)
        Yf2 = Sorting(Ys)

        '''...........................................................'''
        '''This part of code compares the distances of the guards with the previous guards, in the hope of binding them closer'''

        Yf2_len = []
        for i in Yf2:
            Yf2_len.append(len(i))
        print(Yf2_len)

        high = []
        for i in Yf2:
            if len(i) == len(Yf2[0]):
                high.append(Yf2.index(i))

        Dist = []
        for i in high:
            if Yn == []:
                continue
            else:
                a = Yn[len(Yn) - 1][0][0]  # Because the first element has no one to compare with
                b = Yf2[i][0][0][0]  # current elements list
                dist = math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)
                Dist.append(dist)

        '''..........................................................'''
        # if Yn == []:
        #     A2 = Yf2[0]
        # else:
        #     A2 = Yf2[Dist.index(min(Dist))]

        A2 = Yf2[0]
        for i in range(len(A2)):
            Yn.append(A2[i])
        # print(Yn)
        Yf2.remove(A2)
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


# Code perfect................................................#
Yn = mini_chk_pts(Pb, Pc, P, Yx)
print(Yn)


def clean_up_final(Yn):
    final = [];
    R = [];
    r = []
    for i in Yn:  # avoiding repetition
        if not i in final:
            final.append(i)
    for p in range(len(final)):  # solution for adjecent points
        for q in range(len(final)):  # this is a big change!!!!!!!!!
            for r in range(len(Pc) - 1):
                if (final[p][0][0] or final[p][1][0]) == Pc[r]:
                    if (Pc[r + 1] or Pc[r - 1]) == (final[q][0][1] or final[q][1][1]):
                        R.append(final[q])
    for r in range(len(R)):
        if R[r] in final:
            final.remove(R[r])
    Yn = final
    return Yn


# Code perfect................................................#
Final_Diagonals = clean_up_final(Yn)
Yn = Final_Diagonals


def Guards(Final_Diagonals):
    Guards = []
    for i in range(len(Final_Diagonals)):
        if not Final_Diagonals[i][0][0] in Guards:
            Guards.append(Final_Diagonals[i][0][0])
    return Guards


# print(Guards(Final_Diagonals))
# print("The co-ordinates are:",Yn)
# Code perfect................................................#
def plt_plot(P, Yn):
    Px = [];
    Py = [];
    Dx = [];
    Dy = [];
    Sx = [];
    Sy = [];
    APx = [];
    APy = []
    for h in range(len(AAP)):
        APx.append(AAP[h][0])
        APy.append(AAP[h][1])
    for i in range(len(P)):
        Px.append(P[i][0])
        Py.append(P[i][1])
    for j in range(len(Yn)):
        Dx = [];
        Dy = []
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
        plt.plot(Dx, Dy, color='g')
    plt.plot(Px, Py, color='b')
    plt.plot(APx, APy, color='r')
    plt.scatter(Sx, Sy, s=700, marker='.', color='k')
    End = time.time()
    return plt.show(), print("The End time is:", End), print("The runtime is:", (End - Start))


plt_plot(P, Yn)



