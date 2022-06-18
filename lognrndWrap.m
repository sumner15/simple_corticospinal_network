function dataOut = lognrndWrap(m,v,N)
%this is a wrapper function for the lognrnd log-normal distribution matlab
%function. This wrapper will convert normal mean and standard deviation for
%use in the lognrnd funciton using the equations defined below:
%   to generate data from a lognormal distribution with mean M and
%   Variance V, use
%
%      MU = log(M^2 / sqrt(V+M^2))
%      SIGMA = sqrt(log(V/M^2 + 1))
% 
% inputs:
% m = mean
% v = variance
% N is the sampled number of results wanted. 

mu = log(m^2 / sqrt(v+m^2));
sigma = sqrt(log(v/m^2 + 1));
dataOut = lognrnd(mu,sigma,[1,N]); 

end