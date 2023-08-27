function [J, grad] = costFunction(theta, X, y)
%   [J, grad] = costFunction(theta, X, y)
%   Compute cost and gradient for logistic regression.
%   Computes the cost of using theta as the
%   parameter for logistic regression and the gradient of the cost
%   with respect to the parameters.

% Initialize some useful values
m = length(y); % number of training examples

% Compute cost function
J=(1/(m))*(-y'*log(sigmoid((X*theta)))-(1-y)'*log(1-sigmoid((X*theta))));

% Compute gradient
grad = (1/m)*((sigmoid((X*theta))-y)'*X)';

end
