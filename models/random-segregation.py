import numpy as np
import math
import random
import matplotlib.pyplot as plt
#import plotly.plotly as py

# Constants
A = 0
ALPHA = 1

# Mitochondria per cell
initA = 50
initAlpha = 50

# Helper Functions

# grow(n, p)
def grow(p, cell):
	cell[A] = round(1/p*cell[A])
	cell[ALPHA] = round(1/p*cell[ALPHA])
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

# Proportion of mitochondria passed to daughter
buddingProportion = 0.25

percentHomoplasmicInGeneration = []

# Number mitochondria DNA from of a and alpha [a, alpha] 
cells = [[initA, initAlpha]]

# Number of generation to dilute after
diluteAfter = 20

# Number of cells to dilute to
diluteTo = math.pow(2, 10)

# Max Homplasmy to simulation
homplasmicMax = 0.5

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

		# Grow the cells
		motherCell = grow(1 - buddingProportion, motherCell);
		daughterCell = grow(buddingProportion, daughterCell);

		# Check for "dead" (alpha + a = 0) for both mother and daughter
		if (motherCell[ALPHA] + motherCell[A] != 0):
			cellsInNextGeneration.append(motherCell)
		
		if (daughterCell[ALPHA] + daughterCell[A] != 0):
			cellsInNextGeneration.append(daughterCell)


	# Update the cells
	cells = cellsInNextGeneration

	# Dilute if necessary
	if ((generation % diluteAfter) == 0 and generation != 0):
		cells = dilute(cells, int(diluteTo))	

 	percentHomoplasmicInGeneration.append(percentHomoplasmic(cells))

 	if (percentHomoplasmicInGeneration[generation - 1] > homplasmicMax):
 		break

 	generation += 1



# Graph Histogram of Last Generation of Alpha/(Alpha + A)

print generation

# numBins = round(math.sqrt(len(cells)))
# percentAlphaInLastGeneration = percentAlpha(cells)


# mpl_fig = plt.figure()

# plt.hist(percentAlphaInLastGeneration, numBins, normed=0, facecolor='green', alpha=0.5)
# unique_url = py.plot_mpl(mpl_fig, filename="Mitochondria: Random Segregation")

# print unique_url
