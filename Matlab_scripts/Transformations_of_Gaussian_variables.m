%% Transformations of Gaussian random variables

% The purpose of this scenario is to get familiar with a few common types of 
% transformations and to practice how to calculate mean and covariance of transformed variables.
% 
% Suppose that we are interested in the position of an object in Cartesian xy-coordinates,
% denoted  x . However, we are not able to measure  x  directly, but rather some 
% related entity  y=h(x) , where  h(?)  describes the relationship between the position 
% of the object and what we observe.
% 
% Further, assume that we have prior knowledge of the position of the object modelled as,
% x?N([73],[0.2008]). 
% 
% As we will see later in the course, it is often of interest to describe what we already 
% know (our prior on  x ) in the same coordinate system as our observed entity  y . 
% That is, we are interested in the transformed density  p(y)=p(h(x)).

% ----------------------------------------------------------------------------------------------
% Now, let's look at another sensor which reports range and angle to the object. That is
% z=h(x)=[||x||atan2(x2,x1)] 
% where  ||x||  is the vector length  x21+x22??????? , and atan2 gives the angle in the interval  (??,?]  for a right triangle with opposite side  x2  and adjacent side  x1  (i.e. similar to arctan(y/x) but also indicating which quadrant the triangle lies in).
% 
% The mean of  z  will have two elements,  [??,??]T .

%%
mu_x = [7; 3];
sigma_x = [0.2, 0.0; 0.0, 8.0];

%function handles of non-linear mapping from cart to polar coordinates
func_handle =  @custom_cart2pol;

%N param for approx gauss function
num_iterations = 5000;

% approximating function by sampling
[mu_y, Sigma_y, y_s] = approxGaussianTransform(mu_x, sigma_x, func_handle, num_iterations);


fprintf('Mu_y is \n');
disp(mu_y);

fprintf('Sigma_y is \n');
disp(Sigma_y);
