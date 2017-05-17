from numpy import *
from scipy import integrate, optimize
import pylab as p
import random

class ODESystemSimulator(object):
	"""Class used to simulate an N-chemical system
	
	Attributes:
		__systemFunction: A function consisting of all the ODEs
		__systemSize: The number of varaibles in the system
		
	"""
	def __init__(self, function, numVariables):
		"""Constructor
		Args:
			listOfODEs: A list of ODEs, one for each variable.
				They should be in the format: func dX_dt(X, t) return X
		"""
		self.__systemFunction = function
		self.__systemSize = numVariables
	
	def getOutput(self, initialConditions, t=0):
		"""Returns the output of the system at times t, when the system has the
		given initialConditions
		
		Args:
			initialConditions: A list of the same size of the system stating all
			                   the initial conditions of the ODEs.
			t: A list of all the times where an output is wanted. The first time
			   is when the initial conditions are used.
		
		Returns:
			A list of the outputs at times t.
		"""
		if(len(initialConditions) != self.__systemSize):
			print("Inputed " + initialConditions + 
			       " initial conditions but there is only " + self.__systemSize)
		X = integrate.odeint(self.__systemFunction, initialConditions, t)
		return X
		
	def concentrations(self, endTime = 1000, Xs = [0, 1]):
		"""Shows the concentrations of the system over time
		
		Args:
			endTime: Time to end the simulation, 100000 time steps will be used
			Xs: The variables to plot
		"""
		t = linspace(0, endTime,  100000)
		X = self.getOutput([0] * self.__systemSize, t)
		graph = p.figure()
		for Xnum in Xs:
			p.plot(t, X[:, Xnum], ls="--", label=r'$X_' + str(Xnum) + r'$')
		p.title('System Concentrations over time')
		p.xlabel(r"Time, $t$")
		p.ylabel("Variable Concentrations")
		p.legend()
		p.show()
		
	def trajectories(self, listInitialConditions = [[0, 0], [8, 8], [4, 23]], samplingPoints = 20, endTime = 1000, X1 = 0, X2 = 1):
		"""Shows the trajectories the system depending on the initial conditions
		
		Args:
			listInitialConditions: a list of all the intialConditions, ie a list
			                       of arrays.
			samplingPoints: the number of arrows to show in both the x-direction
			                and the y-direction.
			endTime: the time to end the simulation. 100000 time points will be
			         used.
			X1: the x-axis variable
			X2: the y-axis variable
		"""
		fig = p.figure()
		
		t = linspace(0, endTime,  100000)
		
		lineThickness = 1.5**(len(listInitialConditions))
		
		for initialConditions in listInitialConditions:
			lineThickness /= 1.5
			X = self.getOutput(initialConditions, t)
			p.plot( X[:,X1], X[:,X2], lw=lineThickness, ls='--', label='Initial Conditions = (%.f, %.f)' % ( initialConditions[X1], initialConditions[X2]) )

		# define a grid and find the direction at each point
		ymax = p.ylim(ymin=0)[1] # get axis limits
		xmax = p.xlim(xmin=0)[1] # no need to make arrows past axis limits

		x = linspace(0, xmax, samplingPoints)
		y = linspace(0, ymax, samplingPoints)

		X, Y = meshgrid(x, y) # create grid
		DX, DY = self.__systemFunction([X, Y]) # compute derivative of the grid
		M = (hypot(DX, DY)) # Norm of the derivative 
		M[ M == 0] = 1. # Make sure there are no M==1, so no divide by 0 error 
		DX /= M # Normalize each arrow
		DY /= M

		# Draw Directions
		p.title('Trajectories and Directions')
		Q = p.quiver(X, Y, DX, DY, M, pivot='mid')
		p.xlabel(r'$X_' + str(X1) + r'$')
		p.ylabel(r'$X_' + str(X2) + r'$')
		p.legend()
		p.grid()
		p.xlim(0, xmax)
		p.ylim(0, ymax)
		p.show()
		
class Reaction:
	"""A class that stores the information for one reaction
	
	Variables:
		propensity: A function that takes an array of all the variables as
			        an input and outputs the propensity
		changes: An array that shows how all the variables change if the
			     reaction occurs.
	"""
	def __init__(self, propensity, changes):
		"""Constructor
		
		Args:
			propensity: A function that takes an array of all the variables as
			            an input and outputs the propensity
			changes: An array that shows how all the variables change if the
			         reaction occurs.
		"""
		self.propensity = propensity
		self.changes = changes


def gillespieAlgorithm(numberOfVariables, reactions, initialConditions, steps = 500):
	"""A function that calculates and plots the results of Gillespie's Algorithm
	
	Args:
		numberOfVariables: The number of variables in the system. This must
		                   match the expected input of the reactions
		reactions: A list of Reaction objects that show the system reactions
		initialConditions: the initial conditions
		steps: The number of simulation steps that should be used.
	"""
	t = [0]
	X = [initialConditions]
	
	if len(initialConditions) != numberOfVariables:
		print("Length of initialConditions: %.i != numberOfVariables: %.i" %
		      (initialConditions, numberOfVariables))

	for iteration in range(steps):
		currentVals = X[iteration]

		prospensities = []
		for reaction in reactions:
			prospensities.append(reaction.propensity(currentVals))

		totalProspensity = sum(prospensities)
		t.append(-1/totalProspensity * log(1-random.random()) + t[iteration])

		prospensities /= totalProspensity
		subtotalProspensity = 0

		chosenProspensity = random.random() #number in range [0, 1)

		for i in range(len(prospensities)):
			subtotalProspensity += prospensities[i]
			if(subtotalProspensity >= chosenProspensity):
				X.append(array(reactions[i].changes) + array(currentVals))
				break
	X = array(X)
	graph = p.figure()
	for i in range(numberOfVariables):
		p.step(t, X[:, i], lw=0.5, label=r"Gillespie's $X_" + str(i) + r'$')
	p.title('System Concentrations over time')
	p.xlabel(r"Time, $t$")
	p.ylabel("Variable Concentrations")
	p.legend()
	p.show()

def dX_dt(inputs, t=0):
	X = inputs[0]
	Y = inputs[1]
	dX_dt = 1 + 0.02* X**2 *Y - 2*X - 0.04*X
	dY_dt = 2*X - 0.02* X**2 *Y
	return [dX_dt, dY_dt]
	
reactions = [Reaction(lambda X: 1, [1,0]),
             Reaction(lambda X: 2*X[0], [-1,1]),
             Reaction(lambda X: 0.02* X[0]**2 *X[1], [1,-1]),
             Reaction(lambda X: 0.04*X[0], [-1,0])]
gillespieAlgorithm(2, reactions, [0, 0])

system = ODESystemSimulator(dX_dt, 2)
system.concentrations(500)
system.trajectories()
