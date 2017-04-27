function av_n = multisim(P, s, u, v, NoCopies, T)
  n = zeros(NoCopies,1);
  average = zeros(T+1,1);
  clf
  hold on
  xlabel('Number of mutants, n')
  ylabel('P(n)')
  zlabel('t')
  title('Multi-Simulated Value of a Markov Model')
  for t=0:T
    % plot distribution
    av_n(t+1) = mean(n);
    [h,x] = hist(n,[0:P]);
    h = h/NoCopies;
    plot3(x,h,zeros(length(x),1)+t, 'Color', 'blue')
    drawnow

    % selection
    p_s = (1+s)*n./(P+s*n);
    % mutations
    p_sm = (1-v)*p_s + u*(1.0-p_s);
    % sampling
    t = t + 1;
    parfor i=1:NoCopies
      n(i) = binomial_rnd(P, p_sm(i));
    end
  end
  % plot average
  plot3(av_n, zeros(length(average),1), [0:T], 'Color' , 'red');
  hold off

