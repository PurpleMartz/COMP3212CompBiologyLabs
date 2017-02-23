def loadBlosum50():
	with open('blosum50.txt') as f:
		content = f.readlines()
	content = [x.strip() for x in content]
	
	blosumCosts = []
	for x in content:
		lineArray = x.split(' ')
		lineArray = [int(x) for x in lineArray if x]
		blosumCosts.append(lineArray)
	return blosumCosts

blosumCosts = loadBlosum50()
blosumAxis = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 
              'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
def subCost(a, b):
	a = a.upper()
	b = b.upper()
	if a in blosumAxis and b in blosumAxis:
		return blosumCosts[blosumAxis.index(a)][blosumAxis.index(b)]
	else:
		return -123456789 - 42

insertionCost = -8
def nwAlgo(b, a):
	b = " " + b
	a = " " + a
	costMatrix = []
	
	for i in range(len(a)):
		costMatrix.append([])
		#forwards algorithm
		for j in range(len(b)):
			vals = []
			if i > 0 and j > 0:
				vals.append(costMatrix[i-1][j-1] + subCost(a[i],b[j]))
			
			if i > 0:
				vals.append(costMatrix[i-1][j] + insertionCost)
			
			if j > 0:
				vals.append(costMatrix[i][j-1] + insertionCost)
			
			if not vals:
				vals.append(0)
			
			costMatrix[i].append(max(vals))
	
	#backwards algorithm
	done = False
	i = len(a)-1
	j = len(b)-1
	a_match = ""
	b_match = ""
	blank = "-"
	while not done:
		print("i and j:" + str(i) + " " + str(j))
		currentVal = costMatrix[i][j]
		if currentVal is costMatrix[i-1][j-1] + subCost(a[i],b[j]):
			a_match+=a[i]
			b_match+=b[j]
			i-=1
			j-=1
		elif currentVal is costMatrix[i-1][j] + insertionCost:
			a_match+=a[i]
			b_match+=blank
			i-=1
		elif currentVal is costMatrix[i][j-1] + insertionCost:
			a_match+=blank
			b_match+=b[j]
			j-=1
		else:
			print("Mistakes were made - check your code!")
		if not i and not j:
			done = True
	print("a: " + str(b_match[::-1]))
	print("b: " + str(a_match[::-1]))
	return costMatrix

def swAlgo(b, a):
	b = " " + b
	a = " " + a
	costMatrix = []
	
	for i in range(len(a)):
		costMatrix.append([])
		#forwards algorithm
		for j in range(len(b)):
			vals = []
			if i > 0 and j > 0:
				vals.append(costMatrix[i-1][j-1] + subCost(a[i],b[j]))
			
			if i > 0:
				vals.append(costMatrix[i-1][j] + insertionCost)
			
			if j > 0:
				vals.append(costMatrix[i][j-1] + insertionCost)
			
			vals.append(0)
			
			costMatrix[i].append(max(vals))
	
	#backwards algorithm
	done = False
	max_index_j = []
	max_vals = []
	for row in costMatrix:
		max_index_j.append(row.index(max(row)))
		max_vals.append(max(row))
	
	i = max_vals.index(max(max_vals))
	j = max_index_j[i]
	
	a_match = ""
	b_match = ""
	blank = "-"
	while not done:
		print("i and j:" + str(i) + " " + str(j))
		currentVal = costMatrix[i][j]
		if currentVal is costMatrix[i-1][j-1] + subCost(a[i],b[j]):
			a_match+=a[i]
			b_match+=b[j]
			i-=1
			j-=1
		elif currentVal is costMatrix[i-1][j] + insertionCost:
			a_match+=a[i]
			b_match+=blank
			i-=1
		elif currentVal is costMatrix[i][j-1] + insertionCost:
			a_match+=blank
			b_match+=b[j]
			j-=1
		else:
			print("Mistakes were made - check your code!")
		if not costMatrix[i][j]:
			done = True
	print("a: " + str(b_match[::-1]))
	print("b: " + str(a_match[::-1]))
	return costMatrix

matrix = swAlgo("HEAGAWGHEE", "PAWHEAE")

for row in matrix:
	print(row)
