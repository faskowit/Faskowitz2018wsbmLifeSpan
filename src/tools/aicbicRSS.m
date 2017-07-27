function [ aic , bic ] = aicbicRSS( RSS , n , k) 
% get the aic and bic from the residual sum of squares (RSS) 
% minimize these values in theory for better fits
%
% RSS = residual sum squared (yhat - y).^2
% n = number datapoints for x  
% k = free parameters
%
% dependent on gaussian special case for noise:
% https://en.wikipedia.org/wiki/Akaike_information_criterion#Comparison with least squares
% https://en.wikipedia.org/wiki/Bayesian_information_criterion#Gaussian special case

% note form AIC wiki page: "Leave-one-out cross-validation is 
% asymptotically equivalent to the AIC, for ordinary linear regression 
% models."

aic = 2 * k + n * log(RSS) ;
bic = n * log(RSS/n) + k * log(n) ;