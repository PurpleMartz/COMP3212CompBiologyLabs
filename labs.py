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
              
class AminoAcidMutation(object):
	def __init__(self, blosumCostMatrix):
		self.blosumAxis = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 
              'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
		self.blosumCosts = blosumCostMatrix
		self.insertionCost = -8
		
	def subCost(self, a, b):
		a = a.upper()
		b = b.upper()
		if a in self.blosumAxis and b in self.blosumAxis:
			return self.blosumCosts[self.blosumAxis.index(a)][self.blosumAxis.index(b)]
		else:
			return -123456789 - 42 #some tiny number in case we match a "-"
			
	def __forwardsAlgorithm(self, isSW, a, b):
		costMatrix = []
	
		for i in range(len(a)):
			costMatrix.append([])
			for j in range(len(b)):
				vals = []
				if i > 0 and j > 0:
					vals.append(costMatrix[i-1][j-1] + self.subCost(a[i],b[j]))
			
				if i > 0:
					vals.append(costMatrix[i-1][j] + self.insertionCost)
			
				if j > 0:
					vals.append(costMatrix[i][j-1] + self.insertionCost)
			
				if isSW or not vals:
					vals.append(0)
			
				costMatrix[i].append(max(vals))
				
		return costMatrix
		
	def __backwardsAlgorithm(self, isSW, a, b, costMatrix):
		done = False
		a_match = ""
		b_match = ""
		blank = "-"
		i = j = 0
		
		if isSW:
			max_index_j = []
			max_vals = []
			for row in costMatrix:
				max_index_j.append(row.index(max(row)))
				max_vals.append(max(row))
	
			i = max_vals.index(max(max_vals))
			j = max_index_j[i]		
		else:	
			i = len(a)-1
			j = len(b)-1
		
		print("Probability of this match is " + str(costMatrix[i][j]))
		
		while not done:
			currentVal = costMatrix[i][j]
			if currentVal == costMatrix[i-1][j-1] + self.subCost(a[i],b[j]):
				a_match+=a[i]
				b_match+=b[j]
				i-=1
				j-=1
			elif currentVal == costMatrix[i-1][j] + self.insertionCost:
				a_match+=a[i]
				b_match+=blank
				i-=1
			elif currentVal == costMatrix[i][j-1] + self.insertionCost:
				a_match+=blank
				b_match+=b[j]
				j-=1
			else:
				print("Mistakes were made - check your code!")
			if i==0 and j==0 or (isSW and costMatrix[i][j]==0):
				done = True
		print("a: " + str(b_match[::-1]))
		print("b: " + str(a_match[::-1]))
		
	def __matchAlgo(self, isSW, a, b):
		print("Matching " + b + " with " + a)
		b = " " + b
		a = " " + a
		
		costMatrix = self.__forwardsAlgorithm(isSW, a, b)	
		self.__backwardsAlgorithm(isSW, a, b, costMatrix)

	def needlemanWunsch(self, b, a):
		isSW = False
		self.__matchAlgo(isSW, a, b)

	def smithWaterman(self,b, a):
		isSW = True
		self.__matchAlgo(isSW, a, b)
		
	
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

blosumCosts = loadBlosum50()
match = AminoAcidMutation(blosumCosts)
print("Needleman-Wunsch: ")
match.needlemanWunsch("HEAGAWGHEE", "PAWHEAE")
match.needlemanWunsch("PQPTTPVSSFTSGSMLGRTDTALTNTYSAL",
                      "PSPTMEAVEASTASHPHSTSSYFATTYYHL")
print("")
print("Smith-Waterman: ")
match.smithWaterman("HEAGAWGHEE", "PAWHEAE")
match.smithWaterman("MQNSHSGVNQLGGVFVNGRPLPDSTRQKIVELAHSGARPCDISRILQVSNGCVSKILGRY",
                    "TDDECHSGVNQLGGVFVGGRPLPDSTRQKIVELAHSGARPCDISRI")
print("")
print("Hidden Markov Model")
dishonestCasino = HiddenMarkovModel(("123456"), 
                        (HiddenMarkovModel.State((1/6, 1/6, 1/6, 1/6, 1/6, 1/6),
                                                 (9/10, 1/10)),
                         HiddenMarkovModel.State((1/10, 1/10, 1/10, 1/10, 1/10, 1/2),
                                                 (1/10, 9/10))))
dishonestCasino.viterbi("5453525456666664365666635661416626365666621166211311155566351166565663466653642535666662541345464155")
