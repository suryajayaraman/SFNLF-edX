function [x, P] = kalmanFilter(Y, x_0, P_0, A, Q, H, R)
%KALMANFILTER Filters measurements sequence Y using a Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x N] Estimated state vector sequence
%   P           [n x n x N] Filter error convariance
%

%% Parameters
N = size(Y,2);

n = length(x_0);
m = size(Y,1);

%% Data allocation
x = zeros(n,N);
P = zeros(n,n,N);

%Initial state is given by prior estimates
x_prev = x_0;
P_prev = P_0;

for i = 1: N
	[x_hat,  P_hat]  = linearPrediction(x_prev, P_prev, A, Q);
	[x_post, P_post] = linearUpdate(x_hat, P_hat, Y(:, i), H, R);
	
	%stroing the posterior mean and variance in final array
	x(:,i) = x_post;
	P(:,:,i) = P_post;
	
	%replacing x_prev and p_prev with posteriro values
	x_prev = x_post;
	P_prev = P_post;
end   

end


function [x, P] = linearPrediction(x, P, A, Q)
%LINEARPREDICTION calculates mean and covariance of predicted state
%   density using a liear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%
x = A * x;
P = A * P * A' + Q;
end

function [x, P] = linearUpdate(x, P, y, H, R)
%LINEARPREDICTION calculates mean and covariance of predicted state
%   density using a linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] Measurement
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%

% Your code here

vk = y - H * x;
S = H * P * H' + R;
Kalman_gain = P * H' * inv(S);

x = x + Kalman_gain * vk;
P = P - (Kalman_gain * S * Kalman_gain');
end



