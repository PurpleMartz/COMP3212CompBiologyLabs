import tensorflow as tf
import numpy as np
import random
from datetime import datetime

class SecondaryStructFASTA(object):
	"""Class used to store a FASTA object that has
	secondary protein structure information.
	
	Attributes:
		header - Contains the Header of the FASTA file
		sequence - Contains the sequence in the FASTA file
		secondaryStruct - Contains the secondary protein structure information
	
	TODO: This code doesn't meet the FASTA specs:
		see <https://en.wikipedia.org/wiki/FASTA_format> for more info
	"""
	def __init__(self, header, sequence, secondaryStruct):
		"""Constructor, fileName is a string showing the fileName to load """
		self.header = header
		self.sequence = sequence
		self.secondaryStruct = secondaryStruct
	
	aminoAcids = "ARNDCQEGHILKMFPSTWYVBZ"
	secondaryStructKey = "XCEH"
	 
	def createInputOutput(self, windowsSize = 17):
	
		backSize = int( (windowsSize - 1)/2 )
		
		paddedSequence = ((self.sequence[0]*backSize) + self.sequence 
		                 + (self.sequence[-1]*backSize))
		
		X = []
		Y = []
		
		for i in range(backSize, len(self.sequence) + backSize):
		
			Xi = np.zeros( (windowsSize, len(self.aminoAcids)), dtype=np.bool_)
			
			subSequence = paddedSequence[i-backSize: i+backSize+1]
			
			for j in range(windowsSize):
				try:
					aa = self.aminoAcids.index(subSequence[j])
				except ValueError as e:
					print(subSequence[j] + " was not found")
					raise
				Xi[j,aa] = True
			
			X.append(Xi.flatten())
			
			Yi = np.zeros(len(self.secondaryStructKey), dtype=np.bool_)
			ss = self.secondaryStructKey.index(self.secondaryStruct[i-backSize])
			Yi[ss] = True
			
			Y.append(Yi)
			
		X = np.array(X)
		Y = np.array(Y)
		return X, Y

def outputWithoutSecondaryStruct(SecondaryStructFASTAs, fileName):
	""" Writes the SecondaryStructFASTAs without the secondaryStruct data.

		This is so for testing online algorithms.
		
		Args:
			SecondaryStructFASTAs: a list of SecondaryStructFASTA objects.
			fileName: the file to save the data
	"""
	with open(fileName, "w") as f:
		for SecondaryStructFASTA in SecondaryStructFASTAs:
			f.write(">" + SecondaryStructFASTA.header + "\n")
			f.write(SecondaryStructFASTA.sequence + "\n")

def loadSecondaryStructFASTAs(fileName, limit = 2000, seqLengthLimit = 2000):
	""" Loads a list of FASTA objects with secondary structure information.
	
	The file should have the following format:
	>header
	SEQUENCE
	SECONDARY_STRUCTUE_INFORMATION
	>header2
	... etc
	
	Args:
		fileName: the file-name of the list of FASTA objects.
		limit: The max amount of FASTA objects to load.
		seqLengthLimit: the max length of a sequence (some algorithms have limits)
		
	Returns a list of SecondaryStructFASTA objects.
	"""
	secondaryStructFASTAs = []
	
	with open(fileName, "r") as f:
		numberFound = 0
		while True:
			line = f.readline()
			
			if not line: #finished reading
				return secondaryStructFASTAs
			if line[0] != ">": #not the beginning of a FASTA
				print("Didn't find a header, file may be incorrect.")
				continue
			
			header = line[1:].strip() #remove leading >
			sequence = f.readline().strip()
			secondaryStruct = f.readline().strip()
			
			if(len(sequence) >= seqLengthLimit):
				continue #too long so can't use
				
			secondaryStructFASTA = SecondaryStructFASTA(header, sequence,
			                                            secondaryStruct)
			secondaryStructFASTAs.append(secondaryStructFASTA)
			numberFound+=1
			
			if numberFound >= limit:
				return secondaryStructFASTAs

class SecondaryStructSolver(object):
	def runMLP(self, trainingFASTAs, testingFASTAs):
		learning_rate = 0.001
		training_epochs = 15
		display_step = 1
	
		# Network Parameters
		n_hidden_1 = 500 # 1st layer number of features
		n_hidden_2 = 1000 # 2nd layer number of features
		n_input = 17*22 # size of the input (window size 17, 21 entries)
		n_classes = 4 # total classes (H, C, E, and X)
	
		# tf Graph input
		x = tf.placeholder(tf.float32, [None, n_input])
		y = tf.placeholder(tf.float32, [None, n_classes])
	
		# Create model
		def multilayer_perceptron(x, weights, biases):
			# Hidden layer with RELU activation
			layer_1 = tf.add(tf.matmul(x, weights['h1']), biases['b1'])
			layer_1 = tf.nn.relu(layer_1)
			# Hidden layer with RELU activation
			layer_2 = tf.add(tf.matmul(layer_1, weights['h2']), biases['b2'])
			layer_2 = tf.nn.relu(layer_2)
			# Output layer with linear activation
			out_layer = tf.matmul(layer_2, weights['out']) + biases['out']
			return out_layer
	
		# Store layers weight & bias
		weights = {
			'h1': tf.Variable(tf.random_normal([n_input, n_hidden_1])),
			'h2': tf.Variable(tf.random_normal([n_hidden_1, n_hidden_2])),
			'out': tf.Variable(tf.random_normal([n_hidden_2, n_classes]))
		}
		biases = {
			'b1': tf.Variable(tf.random_normal([n_hidden_1])),
			'b2': tf.Variable(tf.random_normal([n_hidden_2])),
			'out': tf.Variable(tf.random_normal([n_classes]))
		}
	
		# Construct model
		pred = multilayer_perceptron(x, weights, biases)

		# Define loss and optimizer
		cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))
		optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)

		# Initializing the variables
		init = tf.global_variables_initializer()
		
		# Launch the graph
		with tf.Session() as sess:
			sess.run(init)

			# Training cycle
			print("Started training at: " + str(datetime.now()))
			print("Press CTRL+C to stop early")
			try:
				for epoch in range(training_epochs):
					avg_cost = 0.
					total_batch = 100
					# Loop over all batches
					for i in random.sample(range(len(trainingFASTAs)), total_batch):
						batch_x, batch_y = trainingFASTAs[i].createInputOutput()
						# Run optimization op (backprop) and cost op (to get loss value)
						_, c = sess.run([optimizer, cost], feed_dict={x: batch_x,
								                                      y: batch_y})
						# Compute average loss
						avg_cost += c / total_batch
					# Display logs per epoch step
					if epoch % display_step == 0:
						print("Epoch:", '%04d' % (epoch+1), "cost=", \
							  "{:.9f}".format(avg_cost))
					
			except KeyboardInterrupt:
				pass #if ctrl+c is pressed stop early
				 
			print("Optimization Finished!")

			# Test model
			correct_prediction = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
			# Calculate accuracy
			accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
		
			test_x, test_y =  testingFASTAs[0].createInputOutput()
			testSize = 100
			for i in range(1, testSize):
				test_x_i, test_y_i = testingFASTAs[i].createInputOutput()
				test_x = np.concatenate( (test_x, test_x_i) )
				test_y = np.concatenate( (test_y, test_y_i) )
			test_x = np.array(test_x)
			test_y = np.array(test_y)
		
			print("Accuracy:", accuracy.eval({x: test_x, y: test_y}))
			
trainingFASTAs = loadSecondaryStructFASTAs("data/seq+ss_train.txt", 200)
testingFASTAs = loadSecondaryStructFASTAs("data/seq+ss_test1199.txt", 200)
secondaryStructSolver = SecondaryStructSolver()
secondaryStructSolver.runMLP(trainingFASTAs, testingFASTAs)
