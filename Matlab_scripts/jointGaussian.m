function [mu, Sigma] = jointGaussian(mu_x, sigma2_x, sigma2_r)
%jointGaussian calculates the joint Gaussian density as defined
%in problem 1.3a. 
%
%Input
%   MU_X        Expected value of x
%   SIGMA2_X    Covariance of x
%   SIGMA2_R    Covariance of the noise r
%
%Output
%   MU          Mean of joint density 
%   SIGMA       Covariance of joint density

% https://www.statlect.com/probability-distributions/normal-distribution-linear-combinations
% z    = Ax + By is general form
% mu_z = A * mu_x + B * mu_y

% fprintf('mu_x size to f is %d X %d \n', size(mu_x,1), size(mu_x,2));
% fprintf('sigma2_x size to f is %d X %d \n', size(sigma2_x,1), size(sigma2_x,2));
% fprintf('sigma2_r size to f is %d X %d \n', size(sigma2_r,1), size(sigma2_r,2));

mu_vector = [mu_x; 0.0];
sigma_mat = [sigma2_x, 0; 0, sigma2_r];

A_mat = [1, 0; 1, 1];
B_vec = [0; 0];
[mu, Sigma] = affineGaussianTransform(mu_vector, sigma_mat, A_mat, B_vec);

% disp(mu);
% disp(Sigma);
end