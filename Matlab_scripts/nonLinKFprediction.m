function [x, P] = nonLinKFprediction(x, P, f, Q, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%

    switch type
        case 'EKF'
            [fx, jac_fx] = f(x);  
            x_new = fx;
            P_new = jac_fx * P * jac_fx' + Q;
            x = x_new;
            P = P_new;           
            
        case 'UKF'
            
            % generate sigma points
            [SP,W] = sigmaPoints(x, P, 'UKF')
            
            % find the approximate mean          
            x_new = zeros(size(x));
            for index = 1 : size(SP,2)
               [predicted_sigma_pts, jac_fx] = f(SP(:, index));  
               x_new = x_new + (predicted_sigma_pts .* W(index));
            end

            % find the approximate covariance
            P_new = Q;
            for index = 1 : size(SP,2)
               [predicted_sigma_pts, jac_fx] = f(SP(:, index));  
               sigma_pt_variance = (predicted_sigma_pts - x_new) * (predicted_sigma_pts - x_new)';
               P_new = P_new +  sigma_pt_variance.* W(index);
            end
            
            x = x_new;
            P = P_new;
            
            
        case 'CKF'
            % generate sigma points
            [SP,W] = sigmaPoints(x, P, 'CKF')
            
            % find the approximate mean          
            x_new = zeros(size(x));
            for index = 1 : size(SP,2)
               [predicted_sigma_pts, jac_fx] = f(SP(:, index));  
               x_new = x_new + (predicted_sigma_pts .* W(index));
            end
            
            % find the approximate covariance
            P_new = Q;
            for index = 1 : size(SP,2)
               [predicted_sigma_pts, jac_fx] = f(SP(:, index));  
               sigma_pt_variance = (predicted_sigma_pts - x_new) * (predicted_sigma_pts - x_new)';
               P_new = P_new +  sigma_pt_variance.* W(index);
            end
            
            x = x_new;
            P = P_new;
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end


function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or using cubature.

%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance

%Output:
%   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
%   W           [1 x 2n+1] UKF, [1 x 2n] UKF. Vector with sigma point weights 
%

n   = size(x,1); %state vector size

    switch type        
        case 'UKF'
            %initialise the variables
            SP = zeros(n, 2*n + 1);
            W = zeros(1, 2*n + 1);
            
            % calculating the weights
            W_0 = 1 - (n/3);
            W(1) = W_0;
            W(2:end) = (1 - W_0) / (2*n);
            
            %sigma points
            P_root = chol(P, 'lower');  
            factor = sqrt(n / (1- W_0));
            
            SP(:, 1) = x;            
            for i = 1:n
                SP(:, i+1)   = x  + factor.* P_root(:,i);
                SP(:, i+1+n) = x  - factor.* P_root(:,i);
            end

                
        case 'CKF'
            %initialise the variables
            SP = zeros(n, 2*n);
            
            % calculating the weights
            weight = 1 / (2*n);
            W  = ones(1, 2*n) .* weight;
                        
            %sigma points
            P_root = chol(P, 'lower');  
            factor = sqrt(n);
                       
            for i = 1:n
                SP(:, i)   = x  + factor.* P_root(:,i);
                SP(:, i+n) = x  - factor.* P_root(:,i);
            end
 
        otherwise
            error('Incorrect type of sigma point')
    end

end