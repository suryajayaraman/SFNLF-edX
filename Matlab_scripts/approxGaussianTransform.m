function [mu_y, Sigma_y, y_s] = approxGaussianTransform(mu_x, Sigma_x, f, N)
%approxGaussianTransform takes a Gaussian density and a transformation 
%function and calculates the mean and covariance of the transformed density.
%
%Inputs
%   MU_X        [m x 1] Expected value of x.
%   SIGMA_X     [m x m] Covariance of x.
%   F           [Function handle] Function which maps a [m x 1] dimensional
%               vector into another vector of size [n x 1].
%   N           Number of samples to draw. Default = 5000.
%
%Output
%   MU_Y        [n x 1] Approximated mean of y.
%   SIGMA_Y     [n x n] Approximated covariance of y.
%   ys          [n x N] Samples propagated through f


if nargin < 4
    N = 5000;
end

% Reference links
% -----------------------------------------------------
% https://www.mathworks.com/help/stats/mvnrnd.html
% https://www.mathworks.com/help/matlab/ref/mean.html
% https://www.mathworks.com/help/matlab/ref/randn.html
% https://www.mathworks.com/help/matlab/matlab_prog/creating-a-function-handle.html
% https://in.mathworks.com/help/stats/normrnd.html#d117e629625
% https://in.mathworks.com/help/matlab/ref/cov.html
% -----------------------------------------------------


% Final solution
% -----------------------------------------------------

fprintf('N is %d\n', N);
fprintf('mu_x size to f is %d X %d \n', size(mu_x,1), size(mu_x,2));

% MU must be a row vector, or must have CASES rows.
% SIGMA must be a symmetric positive semi-definite matrix.
samples = mvnrnd(mu_x', Sigma_x, N)'; 
fprintf('samples size is %d X %d \n', size(samples,1), size(samples,2));

% passing each sample through f(.)
% y0  = f(samples(:,1));
% y_s = y0;
% for i = 2 : N
%     y_s = [y_s, f(samples(:,i))];
% end

% vectorized implementation
[y_s]    = f(samples);

% since cart2pol requries two arguments changing format
% [y_s]    = f(samples(1,:), samples(2,:));

[s1,s2] = size(y_s);
fprintf('y_s size is %d X %d \n', s1, s2);

mu_y = mean(y_s, 2);
fprintf('mu_y size is %d X %d \n', size(mu_y, 1), size(mu_y, 2));
% cov expects each row to be a observation, so passing the transpose
Sigma_y = cov(y_s');
fprintf('Sigma_y size is %d X %d \n', size(Sigma_y, 1), size(Sigma_y, 2))
% -----------------------------------------------------


% -----------------------------------------------------
% Other ideas 
% -----------------------------------------------------

% 1. using a random vector to find output shape of f
% random_input = randn(size(mu_x));
% [test_output] = f(random_input(0), random_input(1));

% 2. Using random matrix of 0 mean and unit std; using sigma_x and mu_x to get samples
% samples = (Sigma_x * randn(size(mu_x,1), N)) + mu_x;

% 3. passing the samples as a matrix to function f(,)
% [y_s]    = f(samples);

% 4. Calling function
% state_size = 2;
% mu_x = randn(state_size, 1);
% Sigma_x = [0.8, 0.2; 0.2, 2.5];
% % disp(mu_x);
% % Sigma_x = randn(state_size, state_size); 
% func_handle = @pol2cart;
% N = 2000;
% [a,b,c] = approxGaussianTransform(mu_x, Sigma_x, func_handle, N);

% 5. initialize output variables to zeros of correct size
% mu_y    = zeros(s1, 1);
% Sigma_y = zeros(s1, s1);
% -----------------------------------------------------
end