# RBE550Project

## Motion Planning Group Project
### Shreyas Kanjalkar, [Rutwik Bonde](https://github.com/Rubo12345)
Before you run the code, to make sure you have all the requirements met do:

```pip install -r requirements.txt```

## Problem Statement: 
Try to find the Watchman Route, when the area of map is given, by using the optimal positions of static guards.

## Approach

For a given "realistic" polygon with or without holes, first find the ideal static guard locations by using the 
geometry of the polygon. In this case, we check the visibility of the edges of the polygon for each vertex and 
choose the vertex which has the highest visibility, and keep on recursively repeating this process
till all the edges of the polygon are guarded. Here our guards are considered *Omnipotent* meaning they have 
360 degree vision and no limit on how far they can look. There are also certain constraints where if the guard 
can see < 10 degrees of the edge, then it is not considered to be guarded by that particular guard. Given the positions
of these guards, if the watchman were to visit all these guard locations at least once then it would mean that the
watchman has guarded all the edges of the polygon and thus successfully completed a watchman route. 

### Art Gallery Problem - Scan Locations:
Art Gallery Problem is a problem to determine the minimum number of scan locations that are required or are sufficient to cover or see every point in the interior of an indoor environment. An indoor environment can be viewed as a polygon with or without holes with a total of n vertices; and scanners as points in the polygon, or on the vertex of the polygon. or on the edge of the polygon. Any point P in the polygon is said to be visible from a scanner G if the line segment joining P and G does not intersect the exterior of the polygon. 

#### Proposed Algorithm to solve the Art Gallery Problem:
  1) Create a polygon
  2) Find a vertex(Vi) in the polygon which scans maximum no. of edges of the polygon
  3) Now search for edges that remain unscanned by the previous vertex (Vi)
  4) Find another vertex (Vj) on the polygon which scans maximum of the remaining unscanned edges.
  5) Continue step 3 and 4 until all th edges are scans
  
  The figure below gives an idea of how the algorithm works:
  ![Proposed Algorithm Flow](https://user-images.githubusercontent.com/79450753/166077173-88e61091-e632-4abb-84a0-ed738710028f.png)

Run the following python code to get the solution for Art Gallery Problem:
- The Shrink_Polygon_AGP.py gives the output as a list of co-ordinates of the scan locations (on the vertices of the polygon)

  `Shrink_Polygon_AGP.py`
- The Final_Voronoi_New_and_Short.py gives the output as a list of co-ordinates of the scan locations (on the vertices of the voronoi diagram of the polygon)

  `Final_Voronoi_New_and_Short.py`
  
  (The detailed instructions for the above codes are given in code itself)
  
### Travelling Salesman Problem - Watchman Route:
Watchman Route is considered a Travelling Salesman Problem, which is defined as Given a list of cities and the distances
between each pair of cities, what is the shortest possible route that visits each city exactly once and returns to the 
origin city?". In our case the cities are replaced with guard locations. In normal TSP, you get the path distance between two
cities directly by connecting them in a straight line. However, for our case, the points are constrained in the polygon
and subsequently subjected to discretization. In order to connect the points, we used Dijkstra and A* search algorithm 
by setting heuristic as eucledian distance between the points. Here we have to consider if a point is inside the polygon 
or not before exploring those points.

After we get the adjacency matrix of path lengths for each guard, we use them as input for the three algorithms we have 
tried to implement: 
1. Brute Force
2. Held-Karp Algorithm
3. Genetic Algorithm

To run the code, please type in 

`python search.py --holes`

to run the code without holes and

`python serach.py --no-holes`

This will print the path length for each Algorithm and the time it took for that algorithm, and the order of traversal.

  The figure below gives an idea of how the TSP works:
  ![TSP_Path](https://user-images.githubusercontent.com/79450753/166120028-24daafd8-d80e-4687-9a23-b80037624b10.png)

Referral Papers:
1. https://math.mit.edu/~apost/courses/18.204_2018/Nicole_Chesnokov_paper.pdf
2. https://members.loria.fr/HEverett/McGill-ISA/Bellairs-2005/rafa.pdf
3. https://arxiv.org/pdf/1708.01461.pdf
4. https://arxiv.org/pdf/1611.05092.pdf
5. https://arxiv.org/pdf/1612.03246.pdf
6. https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.105.1879&rep=rep1&type=pdf
7. http://www.ams.sunysb.edu/~jsbm/papers/watchman.pdf
8. https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6320878
9. https://link.springer.com/content/pdf/10.1007/BF02574671.pdf (WRP)
10. https://link.springer.com/content/pdf/10.1007/PL00009467.pdf (WRP)
11. https://en.wikipedia.org/wiki/Travelling_salesman_problem
12. https://en.wikipedia.org/wiki/Held%E2%80%93Karp_algorithm
