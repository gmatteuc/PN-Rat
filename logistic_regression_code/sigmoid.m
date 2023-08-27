function g = sigmoid(z)
%   g = sigmoid(z) Compute sigmoid function.
%   Computes the sigmoid of z (z can be a matrix, vector or scalar).

% Set steepness variable
k=10; % originally at 1 /// 10000 working well with bool_hemicircle=0 /// 10 working well with bool_hemicircle=0

% Initialize output variable
g = zeros(size(z));

% Compute sigmoid
for i=1:size(z,1)
    for j=1:size(z,2)
        g(i,j)=1/(1+exp(-z(i,j)*k));
    end
end

end
