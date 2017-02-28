from __future__ import division
def loadBlosum50():
	"""Loads the BLOSUM50 Matrix from a text file
	"""
	with open('blosum50.txt') as file:
		content = file.readlines()
	content = [x.strip() for x in content] #removes \n character
	
	blosumCosts = []
	for x in content:
		lineArray = x.split(' ')
		lineArray = [int(x) for x in lineArray if x]
		blosumCosts.append(lineArray)
	return blosumCosts
              
class AminoAcidMutation(object):
	"""Class used for comparing whether protein sequences are from the same species
	"""
	def __init__(self, blosumCostMatrix):
		"""Initializes the class
		Keyword arguments:
		blosumCostMatrix -- the blosum cost matrix to use
		"""
		self.__blosumAxis = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 
              'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
		self.__blosumCosts = blosumCostMatrix
		self.__insertionCost = -8
		
	def __subcost(self, a, b):
		""" Returns the cost between the chars, a and b, using the
		    blosumCostMatrix
		"""
		a = a.upper()
		b = b.upper()
		if a in self.__blosumAxis and b in self.__blosumAxis:
			return self.__blosumCosts[self.__blosumAxis.index(a)][self.__blosumAxis.index(b)]
		else:
			return -123456789 - 42 #some tiny number in case we match a "-"
			
	def __forwardsAlgorithm(self, isSW, a, b):
		""" The Forwards Algorithm for a Dynamic Programming Problem
		Returns the produced cost matrix
		
		Keyword arguments:
		isSW -- if True, this is for the Smith-Waterman algorithm, else Needleman-Wunsch
		a -- The first sequence to match
		b -- The second sequence to match
		"""
		costMatrix = []
	
		for i in range(len(a)):
			costMatrix.append([])
			for j in range(len(b)):
				vals = []
				if i > 0 and j > 0:
					vals.append(costMatrix[i-1][j-1] + self.__subcost(a[i],b[j]))
			
				if i > 0:
					vals.append(costMatrix[i-1][j] + self.__insertionCost)
			
				if j > 0:
					vals.append(costMatrix[i][j-1] + self.__insertionCost)
			
				if isSW or not vals:
					vals.append(0)
			
				costMatrix[i].append(max(vals))
				
		return costMatrix
		
	def __backwardsAlgorithm(self, isSW, a, b, costMatrix):
		""" The Backwards Algorithm for a Dynamic Programming Problem
		
		Keyword arguments:
		isSW -- if True, this is for the Smith-Waterman algorithm, else Needleman-Wunsch
		a -- The first sequence to match
		b -- The second sequence to match
		costMatrix -- The costMatrix made by the forward algorithm
		"""
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
			if currentVal == costMatrix[i-1][j-1] + self.__subcost(a[i],b[j]):
				a_match+=a[i]
				b_match+=b[j]
				i-=1
				j-=1
			elif currentVal == costMatrix[i-1][j] + self.__insertionCost:
				a_match+=a[i]
				b_match+=blank
				i-=1
			elif currentVal == costMatrix[i][j-1] + self.__insertionCost:
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
		""" Returns the probability the fromState will transition to toState """
		return self.states[fromState].transitions[toState]
		
	def emissionProb(self, state, emission):
		""" Returns the probability that emission is emitted while in state """
		return self.states[state].outputs[self.outputs.index(emission)]
		
	def viterbi(self, sequence):
		""" Finds the most likely sequence of states in the given sequence """
		
		sequence = " " + sequence #adds a leading space for initial state
		
		#forwards algoirithm
		initialState = 0
		probabilities = []
		pointer = []
		
		for i in range(len(self.states)):
			pointer.append([None]) #first pointer points to nothing
			if i == initialState: 
				#sets the initialState to have probability of 100%
				probabilities.append([1])
			else:
				#all other states start at 0%
				probabilities.append([0])	

		for t in range(1, len(sequence)):
			#t starts at 1 since first prob is predefined
			for i in range(len(self.states)):
				vals = []
				for j in range(len(self.states)):
					vals.append(probabilities[j][t-1] * self.transitionProb(i, j))
				probabilities[i].append(self.emissionProb(i, sequence[t]) * max(vals))
				#pointer to which state was the previous state
				pointer[i].append(vals.index(max(vals)))
		
		#backwards algorithm
		finalValues = []
		for probs in probabilities:
			finalValues.append(probs[len(sequence)-1])
		#finds the state which ends with the highest probability
		state = finalValues.index(max(finalValues))
		
		stateSequence = ""
		for t in range(len(sequence)-1,-1, -1):
			#go in reverse order ie for(t = size-1; i != -1; i--)
			stateSequence += str(state)
			state = pointer[state][t]
			
		stateSequence = stateSequence[::-1]
		printSequence = ""
		
		#the following loop adds colors to the sequence showing the state
		for t in range(1,len(stateSequence)):
			#t starts at 1 since first char is a space
			if t <= 1 or stateSequence[t-1] != stateSequence[t]:
				# if state changed, changes background color of text
				printSequence+= '\033[4'+stateSequence[t]+';97m'
			printSequence += sequence[t]
		printSequence+='\033[39;49m' #resets colors
		print("Sequence: " + printSequence)

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

phaseLamda = open("phaseLambda.fasta", "r").readlines()[1:] #load all but the first line
phaseLamda = "".join(phaseLamda)
phaseLamda = [x.strip() for x in phaseLamda] #removes \n characters
phaseLamda = "".join(phaseLamda)

cgRich = HiddenMarkovModel(("ATCG"),
                     (HiddenMarkovModel.State((0.2698, 0.3237, 0.2080, 0.1985),
                                              (0.9998, 0.0002)),
                     (HiddenMarkovModel.State((0.2459, 0.2079, 0.2478, 0.2984),
                                              (0.0003, 0.9997)))))
cgRich.viterbi(phaseLamda)
