function [y] = sigmoid(x)
%% Sigmoid function
% logistic function
% based upon: https://en.wikipedia.org/wiki/Sigmoid_function
%
% Omer Bar 2020
%%
y = pi./(1+exp(-x)) - pi/2;
end