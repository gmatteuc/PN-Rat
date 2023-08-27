function [c,p] = predict(theta, X)
% p = predict(theta, X)
% Predict whether the label is 0 or 1 using learned logistic 
% regression parameters theta.
% Computes the predictions for X using a threshold at 0.5 
% (i.e. if sigmoid(theta'*x) >= 0.5, predict 1)

% Initialize some useful values
m = size(X, 1); % number of training examples

% Compute probability of label=1
p=sigmoid((X*theta));

% Treshold it
c=round(p);

end
