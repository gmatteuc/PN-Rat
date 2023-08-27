function [ellipseX, ellipseY, ellipseX_center, ellipseY_center] = get_ellipse_from_2D_gaussian_fit(fitParams, numPoints)
% calculate the eigenvalues and eigenvectors of the covariance matrix
covMatrix = [fitParams(3)^2, fitParams(3)*fitParams(4)*tanh(fitParams(5));
             fitParams(3)*fitParams(4)*tanh(fitParams(5)), fitParams(4)^2];
[eigVectors, ~] = eig(covMatrix);
% compute the angle of the major axis
angle = atan2(eigVectors(2,1), eigVectors(1,1));
% extract the major and minor axes lengths from the fitted parameters (3sigma)
majorAxis = 2.355 * 1 * sqrt(fitParams(3)^2); % updated to FWHM on 01/07/23
minorAxis = 2.355 * 1 * sqrt(fitParams(4)^2); % updated to FWHM on 01/07/23
% get center
ellipseX_center = fitParams(1);
ellipseY_center = fitParams(2);
% generate points for the ellipse
t = linspace(0, 2 * pi, numPoints);
ellipseX = ellipseX_center + majorAxis * cos(t) * cos(angle) - minorAxis * sin(t) * sin(angle);
ellipseY = ellipseY_center + majorAxis * cos(t) * sin(angle) + minorAxis * sin(t) * cos(angle);
end