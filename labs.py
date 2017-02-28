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
		if currentVal == costMatrix[i-1][j-1] + subCost(a[i],b[j]):
			a_match+=a[i]
			b_match+=b[j]
			i-=1
			j-=1
		elif currentVal == costMatrix[i-1][j] + insertionCost:
			a_match+=a[i]
			b_match+=blank
			i-=1
		elif currentVal == costMatrix[i][j-1] + insertionCost:
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
		if currentVal == costMatrix[i-1][j-1] + subCost(a[i],b[j]):
			a_match+=a[i]
			b_match+=b[j]
			i-=1
			j-=1
		elif currentVal == costMatrix[i-1][j] + insertionCost:
			a_match+=a[i]
			b_match+=blank
			i-=1
		elif currentVal == costMatrix[i][j-1] + insertionCost:
			a_match+=blank
			b_match+=b[j]
			j-=1
		else:
			print("Mistakes were made - check your code!")
			return costMatrix
		if not costMatrix[i][j]:
			done = True
	print("a: " + str(b_match[::-1]))
	print("b: " + str(a_match[::-1]))
	return costMatrix
	
class HiddenMarkovModel(object):
	"""A Hidden Markov Model
	
	Attributes
		outputs: A list of possible outputs of the HMM
		states: A list of State objects that have each state's information
	"""
	class State(object):
		""" A State in a Hidden Markov Model
		
		Attributes:
			transitions: A list of probabilities that show which will be the next state
			outputs: A list of probabilites the show the likelihood of an output
		"""
		def __init__(self, outputs, transitions):
			tolerance = 0.01
			self.transitions = transitions
			self.outputs = outputs
			count = 0
			for output in outputs:
				count += output
			if 1+tolerance < count or count < 1-tolerance:
				#print a warning if probabilities are too far off 1
				print("WARNING: Sum of output probabilities is " + str(count))
				
		def __str__(self):
			return ("Outputs: " + str(self.outputs) + "\n" 
			     + "Transitions: " + str(self.transitions))
	
	def __init__(self, outputs, states):
		self.outputs = outputs
		self.states = states
		for i in range(len(states)):
			state = states[i]
			if len(outputs) != len(state.outputs):
				print("Error: State " + str(i) + " has " + str(len(state.outputs))
				      + " outputs, not " + str(len(outputs)) + " outputs.")
			if len(states) != len(state.transitions):
				print("Error: State " + str(i) + " has " + str(len(state.transitions))
				      + " transitions, not " + str(len(states)) + " transitions.")
	
	def __str__(self):
		string = "Outputs: " + str(self.outputs) + "\n"
		for state in self.states:
			string+= "States: " + str(state) + "\n"
		return string
		
	def transitionProb(self, fromState, toState):
		return self.states[fromState].transitions[toState]
		
	def emissionProb(self, state, emission):
		return self.states[state].outputs[self.outputs.index(emission)]
		
	def viterbi(self, sequence):
		#forwards algoirithm
		initialState = 0
		probabilities = []
		pointer = []
		for i in range(len(self.states)):
			pointer.append([None]) #first pointer points to nothing
			if i == initialState:
				probabilities.append([1])
			else:
				probabilities.append([0])	

		for t in range(1, len(sequence)):
			for i in range(len(self.states)):
				vals = []
				for j in range(len(self.states)):
					vals.append(probabilities[j][t-1] * self.transitionProb(i, j))
				probabilities[i].append(self.emissionProb(i, sequence[t]) * max(vals))
				pointer[i].append(vals.index(max(vals)))
		
		#backwards algorithm
		finalValues = []
		for probs in probabilities:
			finalValues.append(probs[len(sequence)-1])
		state = finalValues.index(max(finalValues))
		stateSequence = ""
		for t in range(len(sequence)-1,-1, -1):
			stateSequence += str(state)
			state = pointer[state][t]
		print("Sequence: " + sequence)
		print("  States: " + stateSequence[::-1])
		       
hmm = HiddenMarkovModel(("123456"), 
                        (HiddenMarkovModel.State((1/6, 1/6, 1/6, 1/6, 1/6, 1/6),
                                                 (9/10, 1/10)),
                         HiddenMarkovModel.State((1/10, 1/10, 1/10, 1/10, 1/10, 1/2),
                                                 (1/10, 9/10))))
hmm.viterbi( "5453525456666664365666635661416626365666621166211311155566351166565663466653642535666662541345464155")        
