import numpy as np
import math
import random
import matplotlib.pyplot as plt
import plotly.plotly as py
import sys

title = raw_input()

title = "Random Segregation: Simulation " + title

# Constants
A = 0
ALPHA = 1

# Mitochondria per cell
initA = 25.0
initAlpha = 25.0

# Helper Functions

# grow(n, p)
def grow(cell):
	factor = (initA + initAlpha) / (cell[A] + cell[ALPHA])
	cell[A] = round(cell[A]*factor)
	cell[ALPHA] = round(cell[ALPHA]*factor)
	return cell

def percentHomoplasmic(cells):
	totalCells = 0.0
	homoplasmicCells = 0.0
	for cell in cells:
		totalCells += 1
		if (cell[A] == 0 or cell[ALPHA] == 0):
			homoplasmicCells += 1
	return homoplasmicCells / totalCells

def percentAlpha(cells):
	percentAlphaPerCell = []
	for cell in cells:
		percentAlphaPerCell.append(cell[ALPHA] / (cell[ALPHA] + cell[A]))
	return percentAlphaPerCell

# Dilute cells down to k, where len(cell) >= k
def dilute(cells, k):
	return random.sample(cells, k)

def cellSize(cells):
	sizes = []
	for cell in cells:
		sizes.append(cell[A] + cell[ALPHA])
	return sizes

buddingProportions = np.linspace(0.01, 0.5)

generationsToHomoplasmy = []

print "Budding Proportion,Generations to Homoplasmy"

for buddingProportion in buddingProportions:

	percentHomoplasmicInGeneration = [0]
	                                 
	# Number mitochondria DNA from of a and alpha [a, alpha] 
	cells = [[initA, initAlpha]]

	#  Dilute after the number of cells is greater
	#  than or equal to 10^5
	diluteAfter = 100000

	# Number of cells to dilute to 10^3
	diluteTo = 1000

	# Max Homplasmy to simulation
	homplasmicMax = 1

	generation = 1;

	# For each generation
	while True:

		cellsInNextGeneration = [];

		# For each cell in generation
		for motherCell in cells:
			# run binomial distribution to see how many alpha are selected
			daughterAlphas = np.random.binomial(motherCell[ALPHA], buddingProportion)
			daughterAs = np.random.binomial(motherCell[A], buddingProportion)

			# Create a budding daughter cell
			daughterCell = [daughterAs, daughterAlphas]

		 	# Update the mother for A and Alpha
			motherCell[A] -= daughterCell[A]
			motherCell[ALPHA] -= daughterCell[ALPHA]

			# Check for "dead" (alpha + a = 0) for both mother and daughter
			if (motherCell[ALPHA] + motherCell[A] != 0):
				motherCell = grow(motherCell);
				cellsInNextGeneration.append(motherCell)
			
			if (daughterCell[ALPHA] + daughterCell[A] != 0):
				daughterCell = grow(daughterCell)
				cellsInNextGeneration.append(daughterCell)


		# Update the cells
		cells = cellsInNextGeneration

		# Dilute if necessary
		if (len(cells) >= diluteAfter):
			cells = dilute(cells, diluteTo)

		percent = percentHomoplasmic(cells)

	 	percentHomoplasmicInGeneration.append(percent)



	 	sys.stderr.write('{3} | Generation: {0} | Homoplasmic: {1:.2%} | Mean Cell Size {2:.2f}\r'.format(generation, percent, np.mean(cellSize(cells)), buddingProportion))
		sys.stderr.flush()

	 	# print "Generation %d: Percent Homoplasmic %.2f, Number of Cells %.2f" % (generation, percent, len(cells))

	 	if (percent >= homplasmicMax):
	 		break

	 	generation += 1

	generationsToHomoplasmy.append(generation)
	 	
	print "%.2f,%d" % (buddingProportion, generation)

mpl_fig = plt.figure()

plt.scatter(buddingProportions, generationsToHomoplasmy)
plt.xlabel('Budding Proportion')
plt.ylabel('Generations to Homoplasmy')
plt.title(title)

unique_url = py.plot_mpl(mpl_fig, filename=title)	