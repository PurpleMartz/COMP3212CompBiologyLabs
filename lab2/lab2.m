%% COMP3212 Lab 2

%% Mutant Alle Gene Frequency
% Simulates the gene frequency of a mutant allele in a population of size
% 100, 
% with a selective advantage of s = 0.01
% (s=0 is no advantage s=1 produces twices as many children), 
% a reverse mutation rate (i.e. from the mutant state to the original state) of v = 0.001,
% and a forward mutation rate of u = 0.001.
simulation(100, 0.01, 0.001, 0.001);
%%
% The mutant doesn't take over the population if S = 0 and v >= u, however
% in all other cases it does. Please note that this assumes infinite time,
% which the above simulation does not use.
%% Markov Chain Analysis
W=transition_matrix(10,0.1,0.001, 0.001);
sum(W)
%%
% W is a left stochastic matrix as all the columns sum up to one.
sum(W')
%%
% However, W is not a right stochastic matrix, as all the rows do NOT sum
% up to one.

%%
% The two figures below show the difference between the simulated
% probabilities and the estimated probabilites. The simulated probabilities
% is noisy in comparison to the estimated probabilites, although they are
% roughly the same.
multisim(100, 0.1, 0.001, 0.001, 100, 100);
title('Simulated');
view(22,37);
%%
markov(100, 0.1, 0.001, 0.001, 100);
title('Expected');
view(22,37);

%%
% The steady state distribution (ie, the probability of each state as t-> infinity)
% is given by the eigenvector with eigenvalue 1.
W = transition_matrix(100, 0.1, 0.01, 0.01);
[V,L]=eig(W);
plot(0:100,V(:,1)/sum(V(:,1)));
xlabel('Number of mutants, n');
ylabel('P(n)');
title('Steady-State Probabilites');

%% Diffusion Analysis
% Diffusion approximation works reasonably well, although if values get too
% extreme it tends to go off course.
steady_state(100,0.01,0.01,0.01);
%%
% For P < 50, the approximate diffusion breaks down
steady_state(49,0.01,0.01,0.01);
%%
% Markov Model Analysis requires an Eigenvector of a PxP
% vector to be calculated, the complexity of which is $O(P^3)$ in practice, ie
% it should only be used for small values of P. Diffusion Approximation has
% a linear cost, and therefore can be used for large values of P, although
% as the approximation breaks down for small values of P, Markov Analysis
% should be used then.

%% Real Populations
% The following shows the results of a "real" genetic algorithm
garesults = ga(100,0.1,0.01,0.01,100,100);
%%
% Below is the previous markov model approximation
markovModelResults = multisim(100,0.1,0.01,0.01,100,100);
view(180,0);
camroll(-90)
%%
% As you can see from the below graph, the results are pretty similar
plot(0:100,garesults,0:100,markovModelResults)
xlabel('Number of mutants, n');
ylabel('P(n)');
title('"Real" Genetic Algorithm vs Markov Model');
legend('Genetic Algorithm', 'Markov Model');
%% Function to find the expected probability of a markov model
function markov(P, s, u, v, T)
  average = zeros(T,1);
  p = zeros(P+1,1);
  p(1) = 1;               % initialise probability distribution
  W = transition_matrix(P, s, u, v);
  x = 0:P;
  clf;
  hold on;
  for t = 1:T
    average(t) = mean(x*p);
    p = W*p;
    plot3(x,p,zeros(length(x),1)+t, 'Color', 'blue');
  end
  plot3(average, zeros(length(average),1), [1:T], 'Color' , 'red');
  xlabel('Number of mutants, n');
  ylabel('P(n)');
  zlabel('t');
  title('Expected Value of a Markov Model');
  hold off
end

