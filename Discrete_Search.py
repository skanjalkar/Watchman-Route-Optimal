# Basic searching algorithms

# Class for each node in the grid

from operator import mod

class Node:
    def __init__(self, row, col, is_obs, h):
        self.row = row        # coordinate
        self.col = col        # coordinate
        self.is_obs = is_obs  # obstacle?
        self.g = None         # cost to come (previous g + moving cost)
        self.h = h            # heuristic
        self.cost = None      # total cost (depend on the algorithm)
        self.parent = None    # previous node

def graph(grid):
    white = []; nodes = []; parent_nodes = []
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            node = []
            if grid[i][j]==0:
                node.append(i)
                node.append(j)
                white.append(node)
                nodes.append(node)
                parent_nodes.append(None)
    return nodes, parent_nodes, white 

def backprop(start,goal,nodes,parent_nodes):
    path = []
    path.append(goal)
    while (goal != start):
        if goal == None:     
            found = False
            break
        else:
            goal = parent_nodes[nodes.index(goal)]
            path.append(goal)
            goal = goal
    path.reverse()
    return path

def adjecent_nodes(u):
    return [[u[0]+0,u[1]+1],[u[0]-1,u[1]+1],[u[0]+1,u[1]+0],[u[0]-1,u[1]-1],[u[0]+0,u[1]-1],[u[0]+1,u[1]-1],[u[0]-1,u[1]+0],[u[0]+1,u[1]+1]]

def Steps(black,goal):
    steps = 0
    for i in black:
        if i == goal:
            steps += 1
            found = True
            break
        else:
            steps += 1 
    return steps

def bfs(grid, start, goal):
    '''Return a path found by BFS algorithm 
       and the number of steps it takes to find it.

    arguments:
    grid - A nested list with datatype int. 0 represents free space while 1 is obstacle.
           e.g. a 3x3 2D map: [[0, 0, 0], [0, 1, 0], [0, 0, 0]]   
    start - The start node in the map. e.g. [0, 0]
    goal -  The goal node in the map. e.g. [2, 2]

    return:
    path -  A nested list that represents coordinates of each step (including start and goal node), 
            with data type int. e.g. [[0, 0], [0, 1], [0, 2], [1, 2], [2, 2]]
    steps - Number of steps it takes to find the final solution, 
            i.e. the number of nodes visited before finding a path (including start and goal node)

    >>> from main import load_map
    >>> grid, start, goal = load_map('test_map.csv')
    >>> bfs_path, bfs_steps = bfs(grid, start, goal)
    It takes 10 steps to find a path using BFS
    >>> bfs_path
    [[0, 0], [1, 0], [2, 0], [3, 0], [3, 1]]
    '''
    ### YOUR CODE HERE ###
    path = []
    steps = 0
    found = False
    white = []   #Undiscovered
    gray = []    #Discovered
    black = []   #Visited
    
    nodes, parent_nodes, white = graph(grid)

    gray.append(start)
    white.remove(start)
    parent_nodes[0] = None
    Q = []
    Q.append(start)

    while (Q != []):
        u = Q.pop(0)                   #popping the first node
        Adj_U = adjecent_nodes(u)
        for v in Adj_U:               # for loop to check each neighbor
            if v in white:
                gray.append(v)
                parent_nodes[nodes.index(v)] = u
                white.remove(v)
                Q.append(v)
        black.append(u)
        gray.remove(u)

    if goal not in black:
        found = False
    else:
        path = backprop(start,goal,nodes,parent_nodes)  #back-propagation to get the final path
        found = True
        steps = Steps(black,goal)
    
    if found:
        print(f"It takes {steps} steps to find a path using BFS")
    else:
        print("No path found")
    return path, steps

def dfs_visit(u,nodes,white,black,gray,goal,steps,parent_nodes):
    gray.append(u)
    white.remove(u)
    Adj_U = adjecent_nodes(u)
    if u != goal:
        black.append(u)
    if u == goal:
        found = True
        black.append(goal)
    for v in Adj_U:                    #for loop to check each neighbor
        if v in white:
            parent_nodes[nodes.index(v)] = u
            dfs_visit(v,nodes,white,black,gray,goal,steps,parent_nodes)
        
def dfs(grid, start, goal):
    '''Return a path found by DFS algorithm 
       and the number of steps it takes to find it.

    arguments:
    grid - A nested list with datatype int. 0 represents free space while 1 is obstacle.
           e.g. a 3x3 2D map: [[0, 0, 0], [0, 1, 0], [0, 0, 0]]
    start - The start node in the map. e.g. [0, 0]
    goal -  The goal node in the map. e.g. [2, 2]

    return:
    path -  A nested list that represents coordinates of each step (including start and goal node), 
            with data type int. e.g. [[0, 0], [0, 1], [0, 2], [1, 2], [2, 2]]
    steps - Number of steps it takes to find the final solution, 
            i.e. the number of nodes visited before finding a path (including start and goal node)

    >>> from main import load_map
    >>> grid, start, goal = load_map('test_map.csv')
    >>> dfs_path, dfs_steps = dfs(grid, start, goal)
    It takes 9 steps to find a path using DFS
    >>> dfs_path
    [[0, 0], [0, 1], [0, 2], [1, 2], [2, 2], [2, 3], [3, 3], [3, 2], [3, 1]]
    '''
    ### YOUR CODE HERE ###
    path = []
    steps = 0
    found = False
    white = []   #Undiscovered
    gray = []    #Discovered
    black = []   #Visited
    
    nodes, parent_nodes, white = graph(grid)
    
    for u in nodes:
        if u in white:
            dfs_visit(u,nodes,white,black,gray,goal,steps,parent_nodes)   #dfs visit function
    
    if goal not in black:
        found = False
    else:
        path = backprop(start,goal,nodes,parent_nodes)   #back-propagation to get the final path
        found = True
        steps = Steps(black,goal)

    if None in path:
        path = []
        found = False

    if found:
        print(f"It takes {steps} steps to find a path using DFS")
    else:
        print("No path found")
    return path, steps

def dijkstra(grid, start, goal):
    '''Return a path found by Dijkstra alogirhm 
       and the number of steps it takes to find it.

    arguments:
    grid - A nested list with datatype int. 0 represents free space while 1 is obstacle.
           e.g. a 3x3 2D map: [[0, 0, 0], [0, 1, 0], [0, 0, 0]]
    start - The start node in the map. e.g. [0, 0]
    goal -  The goal node in the map. e.g. [2, 2]

    return:
    path -  A nested list that represents coordinates of each step (including start and goal node), 
            with data type int. e.g. [[0, 0], [0, 1], [0, 2], [1, 2], [2, 2]]
    steps - Number of steps it takes to find the final solution, 
            i.e. the number of nodes visited before finding a path (including start and goal node)

    >>> from main import load_map
    >>> grid, start, goal = load_map('test_map.csv')
    >>> dij_path, dij_steps = dijkstra(grid, start, goal)
    It takes 10 steps to find a path using Dijkstra
    >>> dij_path
    [[0, 0], [1, 0], [2, 0], [3, 0], [3, 1]]
    '''
    ### YOUR CODE HERE ###
    path = []
    steps = 0
    found = False
    Min_distance = []
    infinity = 10**5

    nodes, parent_nodes, white = graph(grid)

    S = []; Q = []
    for i in range(len(nodes)):
        Min_distance.append(infinity)      #setting distance for each node as infinity
        Q.append(nodes[i])

    Min_distance[0] = 0
    w = 1    #initializing the weight
    S = []

    while (Q != []):
        u = Q[Min_distance.index(min(Min_distance))]      # Assigning u the node with minimum distance
        S.append(u)
        if u == goal:
            found = True
            break
        else:
            Adj_U = adjecent_nodes(u)
            for v in Adj_U:            #for loop to check each neighbor                  
                if v in Q:    
                    if Min_distance[Q.index(v)] > (Min_distance[Q.index(u)] + w):
                        Min_distance[Q.index(v)] = Min_distance[Q.index(u)] + w
                        parent_nodes[nodes.index(v)] = u
        Q.remove(u)
        Min_distance.remove(min(Min_distance))

    if goal not in S:
        found = False
    else:
        path = backprop(start,goal,nodes,parent_nodes)    # Back-propagation to get final path
        steps = Steps(S,goal)
    
    if None in path:
        path = []
        found = False

    if found:
        print(f"It takes {steps} steps to find a path using Dijkstra")
    else:
        print("No path found")
    return path, steps

def astar(grid, start, goal):
    '''Return a path found by A* algorithm 
       and the number of steps it takes to find it.

    arguments:
    grid - A nested list with datatype int. 0 represents free space while 1 is obstacle.
           e.g. a 3x3 2D map: [[0, 0, 0], [0, 1, 0], [0, 0, 0]]
    start - The start node in the map. e.g. [0, 0]
    goal -  The goal node in the map. e.g. [2, 2]

    return:
    path -  A nested list that represents coordinates of each step (including start and goal node), 
            with data type int. e.g. [[0, 0], [0, 1], [0, 2], [1, 2], [2, 2]]
    steps - Number of steps it takes to find the final solution, 
            i.e. the number of nodes visited before finding a path (including start and goal node)

    >>> from main import load_map
    >>> grid, start, goal = load_map('test_map.csv')
    >>> astar_path, astar_steps = astar(grid, start, goal)
    It takes 7 steps to find a path using A*
    >>> astar_path
    [[0, 0], [1, 0], [2, 0], [3, 0], [3, 1]]
    '''
    ### YOUR CODE HERE ###
    path = []
    steps = 0
    found = False
    open_list = []; closed_list = []
    nodes, parent_nodes, white = graph(grid)

    g = []; h = []; f = []
    for i in range(len(nodes)):
        g.append(10**5)            # g is cost to come
        h.append((abs(nodes[i][0] - goal[0]) + abs(nodes[i][1] - goal[1])))     # h is the heuristic
        f.append((g[i] + h[i]))    # f is the final cost

    open_list.append(start)
    g[nodes.index(start)] = 0
    h[nodes.index(start)] = abs(start[0] - goal[0]) + abs(start[1] + goal[1])
    f[nodes.index(start)] = g[nodes.index(start)] + h[nodes.index(start)]

    while (open_list != []):
        u = nodes[f.index(min(f))]    #Assigning u the node with minimum cost
        if u not in open_list:
            Found = False
            break
        open_list.remove(u)
        if u == goal:
            closed_list.append(goal)
            found = True
            break
        Adj_U = adjecent_nodes(u)
        for v in Adj_U:          #for loop to check each neighbor
            if v in nodes:
                cost = g[nodes.index(u)] + (abs(u[0] - v[0]) + abs(u[1] - v[1]))
                if v in open_list:
                    if g[nodes.index(v)] < cost:
                        continue
                elif v in closed_list:
                    if g[nodes.index(v)] < cost:
                        continue
                else:
                    open_list.append(v)
                h[nodes.index(v)] = abs(v[0] - goal[0]) + abs(v[1] - goal[1])
                g[nodes.index(v)] = cost
                f[nodes.index(v)] = g[nodes.index(v)] + h[nodes.index(v)]
                parent_nodes[nodes.index(v)] = u    
        f[f.index(min(f))] = 10**5
        closed_list.append(u)

    if goal not in closed_list:
        found = False
    else:
        path = backprop(start,goal,nodes,parent_nodes)   #back-propagation to get the final path
        steps = Steps(closed_list,goal)
    
    if None in path:
        path = []
        found = False

    if found:
        print(f"It takes {steps} steps to find a path using A*")
    else:
        print("No path found")
    return path, steps


# Doctest
if __name__ == "__main__":
    # load doc test
    from doctest import testmod, run_docstring_examples
    # Test all the functions
    testmod()
