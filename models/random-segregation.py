import numpy

cells = numpy.zeros((2048, 2))

# Number of Generations
generations = 10;

# Proportion of mitochondria passed to daughter
buddingProportion = 0.25

# Mitochondria per cell
mtDNAPerCell = 50;

# Number mitochondria DNA from of a and alpha [a, alpha] 
cells[0] = [50,50]

# double for every generation
for generation in range(generations): 
	# iterate through the existing cells:
	for i in range (0, math.pow(2,g)):
		# run binomial distribution to see how many alpha are selected
		numAlpha = np.random.binomial(buddingPercent * mtDNAPerCell, buddingPercent)
		numBeta = buddingPercent * cells[i,0] - numAlpha
		cells[2*i] =  
		numpy.multiply(2,cells)