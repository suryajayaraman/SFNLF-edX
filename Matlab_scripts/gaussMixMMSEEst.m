function [ xHat ] = gaussMixMMSEEst( w, mu, sigma2 )
%GAUSSMIXMMSEEST calculates the MMSE estimate from a Gaussian mixture
%density with multiple components.
%
%Input
%   W           Vector of all the weights
%   MU          Vector containing the means of all components
%   SIGMA2      Vector containing the variances of all components
%
%Output
%   xHat        MMSE estimate

% https://stats.stackexchange.com/questions/309622/calculate-moments-of-a-weighted-mixture-of-normal-distributions
% https://stats.stackexchange.com/questions/445231/compute-mean-and-variance-of-mixture-of-gaussians-given-mean-variance-of-compone
% https://stats.stackexchange.com/questions/207178/estimate-weighted-variances-in-mixture-models

xHat = sum(w.* mu);
disp(xHat);
%YOUR CODE HERE
end

