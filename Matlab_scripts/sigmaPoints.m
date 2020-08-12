function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
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