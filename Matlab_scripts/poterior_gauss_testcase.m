%% posterior Gaussian Test cases
%% Fixed test case with public true values
mu_p = 3.3575;
sigma2_p = 2.2732;
sigma2_r = 3.4704;
y = 5.4159;
mu2 = 4.1722;
sigma2 = 1.3735;

% [mu1,sigma1] = posteriorGaussian(mu_p, sigma2_p, y, sigma2_r);
% tol = 1e-4;
% disp(mu1);
% disp(sigma1);

% assert(abs(mu1-mu2) < tol, 'mean is not correct');
% assert(abs(sigma1-sigma2) < tol, 'variance is not correct');

%% Calculate Joint distribution p(x,y) which is proportional to p(x|y)
syms mux sx sr x y;
[mu_xy, Q_xy] = jointGaussian(mux, sx, sr);

%% Express Gaussian distribution p(x,y)= \propto p(x|y) in terms of x (given y)
syms munew snew c;
eqq1 = ([x;y] - mu_xy ).' * Q_xy^-1  * ([x;y] - mu_xy )     ;
eqq2 = (  x   - munew )  * snew^-1  * (  x   - munew )  + c;
sol = solve( coeffs( eqq1 - eqq2, x) , [munew snew c]);

% Show results
simplify(sol.munew);  % -> (mux*sr + sx*y)/(sr + sx)
simplify(sol.snew);   % -> (sr*sx)/(sr + sx)
simplify(sol.c);      % -> (mux - y)^2/(sr + sx)
 
