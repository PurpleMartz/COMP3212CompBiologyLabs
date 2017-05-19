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

def secondaryStructSolverMLP(trainingFASTAs, testingFASTAs, layerSizes = [50, 50, 50],
           training_epochs = 100, batch_size = 100, learning_rate = 0.1,
           seperateValidationAndTest = True, validationPercentage = 0.2,
           weightCost = 0.05):
	"""Runs a Multi-layer Perceptron Neural Network.
	
	Trains with the trainingFASTAs data, and tests it with
	the testingFASTAs data.
	
	Code loosely based on the code found here:
		https://github.com/aymericdamien/TensorFlow-Examples/blob/master/examples/3_NeuralNetworks/multilayer_perceptron.py
		
	Args:
		trainingFASTAs: A list of SecondaryStructFASTA objects to train with
		testingFASTAs: A list of SecondaryStructFASTA objects to test with
		layerSizes: A list of the number of neurons wanted in each layer
		training_epochs: How long training should go on for
		batch_size: How many sequences will be trained each epoch
		learning_rate: The learning rate of the ADAM optimizer
		seperateValidationAndTest: if True, creates a validation from training 
		                           data, otherwise uses testing data
		validationPercentage: percent of batch_size that should be used
		weightCost: how much having large weights is punished
	"""
	display_step = 1
	
	validationFASTAs = []
	if (seperateValidationAndTest):
		validationEnd = int( len(trainingFASTAs) * validationPercentage )
		shuffle = random.sample(range(len(trainingFASTAs)), len(trainingFASTAs))
		# splits training data into training and validation data
		validationFASTAs = trainingFASTAs[:validationEnd]
		trainingFASTAs = trainingFASTAs[validationEnd:]
	else:
		validationFASTAs = testingFASTAs

	# Network Parameters
	n_input = 17*22 # size of the input (window size 17, 21 entries)
	n_classes = 4 # total classes (H, C, E, and X)
	n_hidden = layerSizes

	# tf Graph input
	x = tf.placeholder(tf.float32, [None, n_input])
	y = tf.placeholder(tf.float32, [None, n_classes])

	# Create model
	def multilayer_perceptron(x, weights, biases):
		prev_layer = x
		for i in range(len(n_hidden)):
			# Hidden layer with RELU activation
			layer = tf.add(tf.matmul(prev_layer, weights[i]), biases[i])
			prev_layer = tf.nn.relu(layer)

		# Output layer with linear activation
		out_layer = tf.matmul(prev_layer, weights[len(n_hidden)]) + biases[len(n_hidden)]
		return out_layer
	
	# Create Weights and Biases for the layers
	weights = []
	biases = []
	prev_layerSize = n_input
	for layerSize in (n_hidden + [n_classes]):
		# weights are from previous layer to this layer
		weights.append(tf.Variable(tf.random_normal([prev_layerSize, layerSize])))
		# biases are for each neuron
		biases.append( tf.Variable(tf.random_normal([layerSize])) )
		prev_layerSize = layerSize

	# Construct model
	pred = multilayer_perceptron(x, weights, biases)

	# Define loss and optimizer
	cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))
	
	# penalizes large weights
	l2_weights_loss = tf.add_n([tf.nn.l2_loss(w) for w in weights])
	weight_loss_gain = weightCost / len(weights)
	cost = cost + weight_loss_gain * l2_weights_loss
	
	optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)

	# Initializing the variables
	init = tf.global_variables_initializer()
	
	# Launch the graph
	with tf.Session() as sess:
		sess.run(init)
		
		# Test model
		correct_prediction = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
		# Calculate accuracy
		accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))

		# Training cycle
		print("MLP Neural Net with Neurons: " + str(n_hidden) 
		      + " weightCost: " + str(weightCost))
		print("Started training at: " + str(datetime.now()))
		print("Press CTRL+C to stop early")
		try:
			prev_valid_cost = 0.0
			for epoch in range(training_epochs):
				avg_cost = 0.
				# Loop over all batches
				for i in random.sample(range(len(trainingFASTAs)), batch_size):
					batch_x, batch_y = trainingFASTAs[i].createInputOutput()
					# Run optimization op (backprop) and cost op (to get loss value)
					_, c = sess.run([optimizer, cost], feed_dict={x: batch_x,
							                                      y: batch_y})
					# Compute average loss
					avg_cost += c / batch_size
					
				# Validate Data
				valid_avg_percentage = 0.
				for i in random.sample(range(len(validationFASTAs)),
				                       int( batch_size*validationPercentage ) ):
					batch_x, batch_y = validationFASTAs[i].createInputOutput()
					c, = sess.run([accuracy], feed_dict={x: batch_x,
					                                y: batch_y})
					valid_avg_percentage += c / int(batch_size*validationPercentage)
					
				# Display logs per epoch step
				if epoch % display_step == 0:
					print("Epoch:", '%04d' % (epoch+1), " train cost =", \
						  "{:.9f}".format(avg_cost))
					print("      correct validations =", \
						  "{:.9f}".format(100*valid_avg_percentage) + "%")
						  
				if(prev_valid_cost >= valid_avg_percentage):
					print("Overfitting so maybe you should stop")
					
				prev_valid_cost = valid_avg_percentage
				
		except KeyboardInterrupt:
			pass #if ctrl+c is pressed stop early
			 
		print("Optimization Finished!")
	
		test_x, test_y =  testingFASTAs[0].createInputOutput()
		testSize = 500
		for i in range(1, testSize):
			test_x_i, test_y_i = testingFASTAs[i].createInputOutput()
			test_x = np.concatenate( (test_x, test_x_i) )
			test_y = np.concatenate( (test_y, test_y_i) )
		test_x = np.array(test_x)
		test_y = np.array(test_y)
	
		print("Accuracy:", accuracy.eval({x: test_x, y: test_y}))
			
trainingFASTAs = loadSecondaryStructFASTAs("data/seq+ss_train.txt", 1000)
testingFASTAs = loadSecondaryStructFASTAs("data/seq+ss_test1199.txt", 500)
secondaryStructSolverMLP(trainingFASTAs, testingFASTAs, [500, 1000])