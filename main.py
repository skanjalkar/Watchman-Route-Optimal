import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import csv

from search import dfs, bfs, dijkstra, astar


# Load map, start and goal point.
def load_map(file_path):
    grid = []
    start = [0, 0]
    goal = [0, 0]
    # Load from the file
    with open(file_path, 'r') as map_file:
        reader = csv.reader(map_file)
        for i, row in enumerate(reader):
            # load start and goal point
            if i == 0:
                start[0] = int(row[1])
                start[1] = int(row[2])
            elif i == 1:
                goal[0] = int(row[1])
                goal[1] = int(row[2])
            # load the map
            else:
                int_row = [int(col) for col in row]
                grid.append(int_row)
    return grid, start, goal


# Draw final results
def draw_path(grid, path, title="Path"):
    # Visualization of the found path using matplotlib
    fig, ax = plt.subplots(1)
    ax.margins()
    # Draw map
    row = len(grid)     # map size
    col = len(grid[0])  # map size
    for i in range(row):
        for j in range(col):
            if grid[i][j]: 
                ax.add_patch(Rectangle((j-0.5, i-0.5),1,1,edgecolor='k',facecolor='k'))  # obstacle
            else:          
                ax.add_patch(Rectangle((j-0.5, i-0.5),1,1,edgecolor='k',facecolor='w'))  # free space
    # Draw path
    for x, y in path:
        ax.add_patch(Rectangle((y-0.5, x-0.5),1,1,edgecolor='k',facecolor='b'))          # path
    ax.add_patch(Rectangle((start[1]-0.5, start[0]-0.5),1,1,edgecolor='k',facecolor='g'))# start
    ax.add_patch(Rectangle((goal[1]-0.5, goal[0]-0.5),1,1,edgecolor='k',facecolor='r'))  # goal
    # Graph settings
    plt.title(title)
    plt.axis('scaled')
    plt.gca().invert_yaxis()


if __name__ == "__main__":
    # Load the map
    grid, start, goal = load_map('map.csv')

    # Search
    bfs_path, bfs_steps = bfs(grid, start, goal)
    dfs_path, dfs_steps = dfs(grid, start, goal)
    dij_path, dij_steps = dijkstra(grid, start, goal)
    aster_path, aster_steps = astar(grid, start, goal)

    # # Show result
    draw_path(grid, bfs_path, 'BFS')
    draw_path(grid, dfs_path, 'DFS')
    draw_path(grid, dij_path, 'Dijkstra')
    draw_path(grid, aster_path, 'A*')
    plt.show()
