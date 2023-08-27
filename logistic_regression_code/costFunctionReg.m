function [J, grad] = costFunctionReg(theta, X, y, lambda)
%   [J, grad] = costFunctionReg(theta, X, y, lambda)
%   Compute cost and gradient for logistic regression with 
%   regularization. Computes the cost of using theta as the
%   parameter for regularized logistic regression and the gradient of the cost
%   with respect to the parameters.

% Initialize some useful values
m = length(y); % number of training examples

% Compute cost function
J=(1/(m))*(-y'*log(sigmoid((X*theta)))-(1-y)'*log(1-sigmoid((X*theta))))+(lambda/(2*m))*sum(theta(2:end).^2);

% Compute gradient
for l=1:length(theta)
    if l==1
       grad(l) = (1/m)*((sigmoid((X*theta))-y)'*X(:,l))';         
    else
       grad(l) = (1/m)*((sigmoid((X*theta))-y)'*X(:,l))' + (lambda/m)*theta(l); 
    end
end

end

