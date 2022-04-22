import math


class Node:
    def __init__(self, x, y, inside, h):
        self.x = x
        self.y = y
        self.inside = inside
        self.h = h
        self.g = 0
        self.cost = math.inf
        self.parent = None