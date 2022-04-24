'''This file is for testing don't delete '''
from shapely.geometry import Polygon
from shapely.geometry import Point
# Poly = [(4,4),(8,4),(8,0),(14,-5),(20,0),(20,6),(15,6),(15,10),\
#     (20,10),(20,14),(16,14),(16,16),(10,16),(10,14),(6,14),(6,16),(2,16)\
#        ,(2,14),(0,14),(-5,7),(0,0),(2,-2),(4,0)]
# Poly = [(0,0),(10,0),(10,1),(10,5),(9,5),(8,1),(8,5),(7,5)\
#         ,(6,1),(6,5),(5,5),(4,1),(4,5),(3,5),(2,1),(2,5),(1,5),(0,1)]
# Poly = [(250,190),(236,192),(207,221),(227,241),(193,275)\
#     ,(173,255),(155,272),(155,302),(164,302),(164,315)\
#     ,(206,315),(206,337),(233,337),(233,311),(257,311)\
#     ,(257,414),(167,414),(167,394),(100,394),(100,414)\
#     ,(43,414),(43,394),(13,394),(13,313),(35,313)\
#     ,(35,343),(62,343),(62,290),(47,290),(47,265)\
#     ,(20,265),(20,286),(0,286),(0,211),(119,211),(124,216)\
#     ,(168,172),(146,149),(164,131),(144,112),(126,130)\
#     ,(90,94),(114,70),(114,48),(195,48),(195,92),(262,92)\
#     ,(262,70),(320,70),(320,97),(365,97),(365,150)\
#     ,(343,150),(343,128),(314,128),(314,192),(343,192)\
#     ,(343,170),(374,170),(374,234),(343,234),(343,260)\
#     ,(283,260),(283,242),(249,242)]


# SP = Shrink_Polygon_AGP
#
# row,col = zip(*P)
# row = list(row)
# col = list(col)
# row.append(row[0])
# col.append(row[0])
# ax.plot(row, col)
# P.reverse()  # already clockwise
# P.append(P[0])
# points = SP.shrink(P)
# points.append(points[0])
# Pb = SP.create_point_pair(P)
# Yx = SP.non_intersecting_diag(points, P, Pb)
# Yn = SP.mini_chk_pts(Pb, points, P, Yx)
# Final_Diagonals = SP.clean_up_final(Yn)
# guards = SP.Guards(Final_Diagonals)
# Guards = tuple(tuple(map(int, tup)) for tup in guards)

import math

from node import Node
from environment import *
import heapq

# Python3 implementation of the above approach
from random import randint

INT_MAX = 2147483647
# Number of cities in TSP
V = 5

# Names of the cities
GENES = "ABCDE"

# Starting Node Value
START = 0

# Initial population size for the algorithm
POP_SIZE = 10

# Structure of a GNOME
# defines the path traversed
# by the salesman while the fitness value
# of the path is stored in an integer


class individual:
	def __init__(self) -> None:
		self.gnome = ""
		self.fitness = 0

	def __lt__(self, other):
		return self.fitness < other.fitness

	def __gt__(self, other):
		return self.fitness > other.fitness


# Function to return a random number
# from start and end
def rand_num(start, end):
	return randint(start, end-1)


# Function to check if the character
# has already occurred in the string
def repeat(s, ch):
	for i in range(len(s)):
		if s[i] == ch:
			return True

	return False


# Function to return a mutated GNOME
# Mutated GNOME is a string
# with a random interchange
# of two genes to create variation in species
def mutatedGene(gnome):
	gnome = list(gnome)
	while True:
		r = rand_num(1, V)
		r1 = rand_num(1, V)
		if r1 != r:
			temp = gnome[r]
			gnome[r] = gnome[r1]
			gnome[r1] = temp
			break
	return ''.join(gnome)


# Function to return a valid GNOME string
# required to create the population
def create_gnome():
	gnome = "0"
	while True:
		if len(gnome) == V:
			gnome += gnome[0]
			break

		temp = rand_num(1, V)
		if not repeat(gnome, chr(temp + 48)):
			gnome += chr(temp + 48)

	return gnome


# Function to return the fitness value of a gnome.
# The fitness value is the path length
# of the path represented by the GNOME.
def cal_fitness(gnome):
	mp = [
		[0, 2, INT_MAX, 12, 5],
		[2, 0, 4, 8, INT_MAX],
		[INT_MAX, 4, 0, 3, 3],
		[12, 8, 3, 0, 10],
		[5, INT_MAX, 3, 10, 0],
	]
	f = 0
	for i in range(len(gnome) - 1):
		if mp[ord(gnome[i]) - 48][ord(gnome[i + 1]) - 48] == INT_MAX:
			return INT_MAX
		f += mp[ord(gnome[i]) - 48][ord(gnome[i + 1]) - 48]

	return f


# Function to return the updated value
# of the cooling element.
def cooldown(temp):
	return (90 * temp) / 100


# Comparator for GNOME struct.
# def lessthan(individual t1,
#			 individual t2)
# :
#	 return t1.fitness < t2.fitness


# Utility function for TSP problem.
def TSPUtil(mp):
	# Generation Number
	gen = 1
	# Number of Gene Iterations
	gen_thres = 10

	population = []
	temp = individual()

	# Populating the GNOME pool.
	for i in range(POP_SIZE):
		temp.gnome = create_gnome()
		temp.fitness = cal_fitness(temp.gnome)
		population.append(temp)

	print("\nInitial population: \nGNOME	 FITNESS VALUE\n")
	for i in range(POP_SIZE):
		print(population[i].gnome, population[i].fitness)
	print()

	found = False
	temperature = 10000

	# Iteration to perform
	# population crossing and gene mutation.
	while temperature > 1000 and gen <= gen_thres:
		population.sort()
		print("\nCurrent temp: ", temperature)
		new_population = []

		for i in range(POP_SIZE):
			p1 = population[i]

			while True:
				new_g = mutatedGene(p1.gnome)
				new_gnome = individual()
				new_gnome.gnome = new_g
				new_gnome.fitness = cal_fitness(new_gnome.gnome)

				if new_gnome.fitness <= population[i].fitness:
					new_population.append(new_gnome)
					break

				else:

					# Accepting the rejected children at
					# a possible probability above threshold.
					prob = pow(
						2.7,
						-1
						* (
							(float)(new_gnome.fitness - population[i].fitness)
							/ temperature
						),
					)
					if prob > 0.5:
						new_population.append(new_gnome)
						break

		temperature = cooldown(temperature)
		population = new_population
		print("Generation", gen)
		print("GNOME	 FITNESS VALUE")

		for i in range(POP_SIZE):
			print(population[i].gnome, population[i].fitness)
		gen += 1


if __name__ == "__main__":

	mp = [
		[0, 2, INT_MAX, 12, 5],
		[2, 0, 4, 8, INT_MAX],
		[INT_MAX, 4, 0, 3, 3],
		[12, 8, 3, 0, 10],
		[5, INT_MAX, 3, 10, 0],
	]
	TSPUtil(mp)



''' 
This code is used for random polygon testing
    # P = draw_polygon(ax, 8, lim_x,  lim_y)
    # vor = Voronoi(P)
    # for i in vor.vertices:
    #     if i[0] < 0 or i[1] < 0 or i[0] > lim_x or i[1] > lim_y:
    #         continue
    #     x = int(i[0])
    #     y = int(i[1])
    #     vor_int.append((x, y))
'''
'''
lim_x, lim_y = 150, 150
    xr = random.randint(5,10)
    P = draw_polygon(ax, xr, lim_x,  lim_y)
    vor = Voronoi(P)
    final_guards = []
    for i in vor.vertices:
        if i[0] < 0 or i[1] < 0 or i[0] > lim_x or i[1] > lim_y:
            continue
        x = int(i[0])
        y = int(i[1])
        final_guards.append((x, y))
        '''

# print(track)
    # print()
    # print(back_track)
    # print(time.time()-s)
    # point_list = [i for i in combinations(watchman_route_pts, 2)]
    # print(point_list)
