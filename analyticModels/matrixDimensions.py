import math

kAlpha = 50
kA = 50

total = 50.0

for i in range(kAlpha + 1):
	for j in range(kA + 1):
		if (i == 0 and j == 0):
			continue
		gamma = total / (i + j)
		a = math.floor(gamma*i + 0.5)
		b = total - a
		print [a, b, a + b]