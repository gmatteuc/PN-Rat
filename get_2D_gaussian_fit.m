function [fitParams, fitGaussian, bfCost] = get_2D_gaussian_fit(data,regpars)

% define 2D Gaussian function
gaussian2D = @(params, x, y) ...
    exp(...
    -(    (x - params(1)).^2/(2*params(3)^2) + ...
          (y - params(2)).^2/(2*params(4)^2) - ...
          2*tanh(params(5))*(x - params(1)).*(y - params(2))/(2*params(3)*params(4))) ...
    );
% find the initial guess for the parameters
[~, idx] = max(data(:));
[row_idx, col_idx] = ind2sub(size(data), idx);
mu_y_guess = row_idx;
mu_x_guess = col_idx;
% sigma_x_guess = min(std(data(row_idx,:)), 0.33 * size(data, 2));
% sigma_y_guess = min(std(data(:,col_idx)), 0.33 * size(data, 1));
sigma_x_guess = min(std(data(row_idx,:)), 0.1 * size(data, 2));
sigma_y_guess = min(std(data(:,col_idx)), 0.1 * size(data, 1));
row_data = data(row_idx, :);
col_data = data(:, col_idx);
col_data_resized = interp1(1:size(data, 1), col_data, linspace(1, size(data, 1), size(data, 2)));
corr_coeff = corrcoef(row_data, col_data_resized);
cov_guess = atanh(corr_coeff(1, 2));
if isnan(cov_guess)
    cov_guess = 0;
end
% initial parameters
initialParams = [mu_x_guess, mu_y_guess, sigma_x_guess, sigma_y_guess, cov_guess];
% define the objective function for least squares fitting
[X, Y] = meshgrid(1:size(data, 2), 1:size(data, 1));
% penalty function (regularization)
% penalty_fun = @(params, lambda) ...
%     (...
%     lambda(3).*abs(1/(params(3))) + ...
%     lambda(4).*abs(1/(params(4))) + ...
%     lambda(5).*abs((tanh(params(5)))) + ...
%     nanmean(lambda).*abs(mu_x_guess-params(1))+abs(mu_y_guess-params(2)) ...
%     );
penalty_fun = @(params, lambda) ...
    (...
    lambda(3).*abs((params(3)-2)).^2 + ...
    lambda(4).*abs((params(4)-2)).^2 + ...
    lambda(5).*abs((tanh(params(5)))).^2 + ...
    nanmean(lambda).*abs(mu_x_guess-params(1))+abs(mu_y_guess-params(2)).^2 ...
    );
% overall cost function
objFun_gen = @(params, lambda) ( nansum(nansum((gaussian2D(params, X, Y) - data).^2)) + penalty_fun(params, lambda) );
objFun =  @(params) objFun_gen(params, regpars);

% perform the nonlinear constrained optimization
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');
lb = [1-0.5, 1-0.5, 0, 0, -Inf]; % set lower bounds for the parameters
ub = [size(data, 1)+0.5, size(data, 2)+0.5, size(data, 1), size(data, 2), Inf]; % set upper bounds for the parameters
fitParams = fmincon(objFun, initialParams, [], [], [], [], lb, ub, [], options);

% create the fitted 2D Gaussian function
fitGaussian = @(x, y) gaussian2D(fitParams, x, y);

% evaluate best fit cost
bfCost = objFun_gen(fitParams,zeros(size(regpars)));

end