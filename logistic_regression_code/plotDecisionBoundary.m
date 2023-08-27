function plotDecisionBoundary(theta, X, y)
% plotDecisionBoundary(theta, X, y) 
% Plots the data points X and y into a new figure with
% the decision boundary defined by theta.
%   Plots the data points with + for the 
%   positive examples and o for the negative examples. X is assumed to be 
%   a either 
%   1) Mx3 matrix, where the first column is an all-ones column for the 
%      intercept.
%   2) MxN, N>3 matrix, where the first column is all-ones

% Plot Data
plotData(X(:,[2,25]), y);
hold on

    % Only need 2 points to define a line, so choose two endpoints
    plot_x = [min(X(:,2)),  max(X(:,2))];

    % Calculate the decision boundary line
    plot_y = (-1./theta(2)).*(theta(1).*plot_x + theta(1));

    % Plot, and adjust axes for better viewing
    plot(plot_x, plot_y,'--k','Linewidth',3)
    
    % Legend, specific for the exercise
    legend('Right target', 'Left target', 'Decision Boundary')
hold off

end
