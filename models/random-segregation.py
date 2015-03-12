import numpy as np

# Constants
A = 0
ALPHA = 1

# Mitochondria per cell
mtDNAPerCell = 100;

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
		if (cell[A] == mtDNAPerCell or cell[ALPHA] == mtDNAPerDaughter):
			homoplasmicCells += 1
	return homoplasmicCells / totalCells

# Number of Generations
generations = 30;

# Proportion of mitochondria passed to daughter
buddingProportion = 0.25

percentHomoplasmicInGeneration = generations*[0];

# Mitochondria per Daughter cell
mtDNAPerDaughter = round(buddingProportion * mtDNAPerCell)

# The probability Alpha will be in daughter
probabilityOfAlphaInDaughter = 0.5

# Number mitochondria DNA from of a and alpha [a, alpha] 
cells = [2*[mtDNAPerCell*0.5]]

# For each generation
for generation in range(generations):

	cellsInNextGeneration = [];

	# For each cell in generation
	for motherCell in cells:
		# run binomial distribution to see how many alpha are selected
		daughterAlphas = np.random.binomial(mtDNAPerDaughter, probabilityOfAlphaInDaughter)


		# Create a budding daughter cell
		daughterCell = [0, 0]


	 	# If the number of potiential Alphas is greater than 
	 	# the number of Alphas in the mother 
	 	if (daughterAlphas > motherCell[ALPHA]):
	 		# Select all the Alphas from the Mother cell
	 		daughterCell[ALPHA] = motherCell[ALPHA]

	 		# Then, the number of As in the daughter is 
	 		# the total number of mtDNA - the number of Alphas in the daughter
	 		daughterCell[A] = mtDNAPerDaughter - daughterCell[ALPHA]

	 	else:
	 		daughterCell[ALPHA] = daughterAlphas

	 		# Potential number of As in the daughter cell
	 		daughterAs = mtDNAPerDaughter - daughterCell[ALPHA]

	 		# If the number of potential As is greater than
	 		# the number of As in the mother, we can take what 
	 		# we have in the mother, and get the rest from the
	 		# Alphas
	 		daughterCell = \
	 		[motherCell[A], mtDNAPerDaughter - motherCell[A]] if \
	 		(daughterAs > motherCell[A]) else \
	 		[daughterAs, daughterCell[ALPHA]]

	 	# Update the mother for A and Alpha
		motherCell[A] -= daughterCell[A]
		motherCell[ALPHA] -= daughterCell[ALPHA]

		# Grow the cells up to mtDNAPerCell
		motherCell = grow(1 - buddingProportion, motherCell);
		daughterCell = grow(buddingProportion, daughterCell);

		cellsInNextGeneration.append(motherCell)
		cellsInNextGeneration.append(daughterCell)

	# Update the cells
	cells = cellsInNextGeneration
	percentHomoplasmicInGeneration[generation] = percentHomoplasmic(cells)

print percentHomoplasmicInGeneration;