import math
from shapely.geometry import Point

def grid(limx, limy, polygon):
    grid = []
    for i in range(0, limx+1):
        row = []
        for j in range(0, limy+1):
            if polygon.contains(Point(i, j)):
                row.append(1)
            else:
                row.append(0)
        grid.append(row)
    return grid



def valid(x, y, grid, visited):
    if (x >= 0 and y >= 0) and (x < len(grid) and y < len(grid[0])):
        if grid[x][y] == 1 and visited[x][y] is False:
            return True
    return False


def travel(node, grid, visited):

    possible_travel = []
    if valid(node.x, node.y+1, grid, visited):
        possible_travel.append([node.x, node.y+1])

    if valid(node.x+1, node.y, grid, visited):
        possible_travel.append([node.x+1, node.y])

    if valid(node.x, node.y-1, grid, visited):
        possible_travel.append([node.x, node.y-1])

    if valid(node.x-1, node.y, grid, visited):
        possible_travel.append([node.x-1, node.y])

    if valid(node.x-1, node.y-1, grid, visited):
        possible_travel.append([node.x-1, node.y-1])

    if valid(node.x+1, node.y-1, grid, visited):
        possible_travel.append([node.x+1, node.y-1])

    if valid(node.x-1, node.y+1, grid, visited):
        possible_travel.append([node.x-1, node.y+1])

    if valid(node.x+1, node.y+1, grid, visited):
        possible_travel.append([node.x+1, node.y+1])
    return possible_travel  # returns all possible node locations that you can travel to


def not_seen_min_cost(node, seen):
    min_node = None  # initialize min node
    for i in range(len(node)):
        for j in range(len(node[i])):
            if seen[i][j] or node[i][j].cost is None:  # if true
                continue
            if min_node is None:  # updating node
                min_node = node[i][j]
            elif min_node.cost > node[i][j].cost:  # cost
                min_node = node[i][j]
            elif min_node.cost == node[i][j].cost:
                if min_node.h > node[i][j].h:  # comparing h when cost is same
                    min_node = node[i][j]

    return min_node

def distance(start, end):
    distance = math.hypot(abs(end.x-start.x), abs(end.y - start.y))
    return distance

