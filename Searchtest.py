class Search:

    def __init__(self, start, goal, grid):
        self.start = start
        self.goal = goal
        self.grid = grid

        self.OPEN = []
        self.CLOSED = []
        self.PARENT = dict()