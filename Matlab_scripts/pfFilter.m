function [xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, f, Q, h, R, N, bResample, plotFunc)
    % PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
    % state-space model.
    %
    % Input:pfFilter(x_0
    %   x_0         [n x 1] Prior mean or 
    %            OR [n x N] initial particle positions
    %   P_0         [n x n] Prior covariance
    %   Y           [m x K] Measurement sequence to be filtered
    %   f           Handle for process function f(x_k-1)
    %   Q           [n x n] process noise covariance
    %   h           Handle for measurement model function h(x_k)
    %   R           [m x m] measurement noise covariance
    %   N           Number of particles
    %   bResample   boolean false - no resampling, true - resampling
    %   plotFunc    Handle for plot function that is called when a filter
    %               recursion has finished.
    % Output:
    %   xfp         [n x K] Posterior means of particle filter
    %   Pfp         [n x n x K] Posterior error covariances of particle filter
    %   Xp          [n x N x K] Particles for posterior state distribution in times 1:K
    %   Wp          [N x K] Non-resampled weights for posterior state x in times 1:K
    
    n = size(x_0,1);
    K = size(Y,2);

    % allocate memory
    xfp = zeros(n,K);
    Pfp = zeros(n,n,K);
    Xp = zeros(n,N,K);
    Wp = zeros(N,K);
    
    % sample initial particles around prior distribution
    % if x_0 has only one column, x_0 is the mean of the prior
    % otherwise x_0 are the initial particles states
    if size(x_0,2) == 1
        Xp(:,:,1) = mvnrnd(x_0,P_0,N)';
    else
        Xp(:,:,1) = x_0;
    end
    
    Wp(:,1)  = 1/N * ones(1,N);
    
    j = 1:N;
    for k=2:K+1
        
        Xp_km1 = Xp(:,:,k-1);
        Wp_km1 = Wp(:,k-1)';
        % resample
        if bResample
            [Xp_km1, Wp_km1, j] = resampl(Xp_km1, Wp_km1);        
        end
            
        % perform a particle filter step for the next measurement
        [Xp(:,:,k), Wp(:,k)] = pfFilterStep( Xp_km1, Wp_km1, Y(:,k-1), f, Q, h, R);
        % plot particles using function handle
        if ~isempty(plotFunc)
            plotFunc(k-1, Xp(:,:,k), Xp(:,:,k-1), Wp(:,k)', j);
        end
%         % if handle bmap is given: update p(y_k|x(i)_k,M)=p(y_k|x(i)_k)*bmap(M|x(i)_k)
%         if ~isempty(pmap)
%             p_map_x = pmap(Xp(1,:,k),Xp(2,:,k));
%             Wp(:,k) = Wp(:,k) .* p_map_x;
%             Wp(:,k) = Wp(:,k)/sum(Wp(:,k));
%         end
% % % %         % resample
% % % %         if bResample
% % % %             [Xp(:,:,k), Wp(:,k), j] = resampl(Xp(:,:,k), Wp(:,k)');
% % % %         end
        % estimate mean and covariance given the particles
        xfp(:,k)   = sum( Xp(:,:,k).*Wp(:,k)' ,2 );
        Pfp(:,:,k) = Wp(:,k)'.*(Xp(:,:,k) - xfp(:,k))*(Xp(:,:,k) - xfp(:,k))';
    end
    
    % remove prior from vector
    xfp = xfp(:,2:end);
    Pfp = Pfp(:,:,2:end);
    Xp  = Xp(:,:,2:end);
    Wp  = Wp(:,2:end);
end


function [Xr, Wr, j] = resampl(X, W)
%RESAMPLE Resample particles and output new particles and weights.
% resampled particles. 
%
%   if old particle vector is x, new particles x_new is computed as x(:,j)
%
% Input:
%   X   [n x N] Particles, each column is a particle.
%   W   [1 x N] Weights, corresponding to the samples
%
% Output:
%   Xr  [n x N] Resampled particles, each corresponding to some particle 
%               from old weights.
%   Wr  [1 x N] New weights for the resampled particles.
%   j   [1 x N] vector of indices refering to vector of old particles

% Your code here!
    N=size(X,2);
    
    % Generates the segmented numberline from 0 to 1
    segment = [0 cumsum(W)/sum(W)];
    
    % draw samples from uniform distribution on [0,1]
    samples = rand([1 N]);
    
    j=zeros(1,N);
    for i=1:N
        j(i) = find(samples(i) >= segment,1,'last');
    end
    
    Wr = 1/N*ones(1,N);
    Xr = X(:,j);
    
end

function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, yk, proc_f, proc_Q, meas_h, meas_R)



%PFFILTERSTEP Compute one filter step of a SIS/SIR particle filter.
%
% Input:
%   X_kmin1     [n x N] Particles for state x in time k-1
%   W_kmin1     [1 x N] Weights for state x in time k-1
%   y_k         [m x 1] Measurement vector for time k
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%
% Output:
%   X_k         [n x N] Particles for state x in time k
%   W_k         [1 x N] Weights for state x in time k

% Your code here!
   % draw samples p(x(i)_k|x(i)_{k-1})
    X_k = mvnrnd( proc_f(X_kmin1)' , proc_Q )';  
    
    % calculate p(y_k|x(i)_k) 
    Wy = mvnpdf(yk', meas_h(X_k)', meas_R)';
    
    % compute weights
    W_k = W_kmin1 .* Wy;
    W_k = W_k / sum(W_k);

end
