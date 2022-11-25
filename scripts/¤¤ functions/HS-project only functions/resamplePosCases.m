function Jresample = resamplePosCases(J,K)
% This function resamples from positive-class in a way that minimizes
% randomness. It first resamples every positive element deterministically,
% and only resamples stochasticly when it can no longer resample each
% positive element one more timewithout exceeding K. Outputs a
% linear-index-vector with resampling indeces. J contains positions of
% positive cases, and K is the number of times to resample from the
% positive class.

% example input: J = JposCases; K = Nbalance;
%%
Nsample = numel(J);
% compute how many times each obs. gets sampled alteast once:
NatleastOnce = floor(K/Nsample);
Nresample = K - NatleastOnce*Nsample;
% sample without replacement:
Jsample = randsample(J,Nresample);
% create string of copies of J
Jcopies = repmat(J,NatleastOnce,1);
Jresample = [Jcopies;Jsample];
% compute how many times to resample randomly:

end