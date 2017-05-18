#' % Computational Biology Lab 3
#' % Alois Klink
#' % 18 May 2017

#' # Converting Reaction Equations to a ODE

#' To convert many reaction equations to one ODE, one must first find the propensity
#' and the changes of each reaction.

#' The Reaction class takes a lambda function of the propensity and the change matrix
#' as inputs.

from simulateChemicals import *

#' $ \varnothing \xrightarrow{1} X $

reactions = [Reaction(lambda X: 1, [1,0]),

#' $ X \xrightarrow{2} Y $

             Reaction(lambda X: 2*X[0], [-1,1]),
             
#' $ 2 X + Y \xrightarrow{0.02} 3 X $

             Reaction(lambda X: 0.02* X[0]**2 *X[1], [1,-1]),
             
#' $ X \xrightarrow{0.04} \varnothing $
             
             Reaction(lambda X: 0.04*X[0], [-1,0])]

#' # Displaying the ODE

#' The figure below shows how the system described by the above reactions
#' behaves when modelled with an ODE. Notice that as X is being created, it is
#' immediatly turned into Y. However, once Y passes a threshold point, and starts
#' combining with X, due to the X^2 factor, X dramatically jumps up, rapidly
#' converting all Y to X. Once Y runs out, X slowly begins to degrade to
#' an equilibrium position.

#+ trajectories, caption='ODE Simulation Trajectories from 0 initialConditions'
system = ChemicalReactionsSystem(reactions, 2)
system.trajectories()

#' # Gillespie's Algorithm

#' The code below shows how the reaction changes when Gillespie's algorithm is
#' used to simulate the reactions. Gillespie's algorithm can be reduced by a 
#' factor to increase the accuracy of the algorithm. This technically works by
#' increasing the number of molecules, and speeding up reactions. However, these
#' graphs have molecules split into pieces, which is not possible in the real world.

for r in [1]:
	reactions = Reaction.rescaleReactions(reactions, r)
	system = ChemicalReactionsSystem(reactions, 2)
	
#' This shows how the cocentration changes over time. Notice that due the high
#' randomness of the properties, the threshold point is reached much faster. As
#' r increases, however, and the Gillespie's algorithm is reduced, the variance
#' gets smaller and smaller, so that the threshold point is reached at the same
#' time as the ODE.

	print("r is " + str(r))
	system.gillespieConcentrations(500*r)

#' This shows how the cocentration X changes in relation to the concentrations Y.
#' Notice that due the high randomness of the properties, the path is a lot tighter,
#' and the stable point seems to be a lot lower.
	system.gillespieTrajectories([[0, 0], [8, 8], [4, 23]],
			                     10000*r)
