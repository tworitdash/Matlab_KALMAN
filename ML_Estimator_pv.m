%
% Pick some numbers for postition and speed of the object. Try several
% numbers to see how things are influenced...
%

p = 50;
v = 10;

x = [p;
     v];

%
% Pick a value for the update intervals between the measurements (in
% this case I selected equal intervals, but you could also try with a 
% variation, like T1 between measurements 1 and 2, and T2
% between measurements 2 and 3. Or... add a fourth measurement T seconds
% separated from measurement 3. 
%
%


T = 2;

t0 = 0;
t1 = T;
t2 = 2*T; 

tb = t1;

%
% Pick some number for the std of the measurement noise. Try different
% numbers to see their effect.
% 

sigma_n = 3;

var_n_t0 = sigma_n^2;
var_n_t1 = sigma_n^2;
var_n_t2 = sigma_n^2;

% The covariance matrix now is. But you could try correleated measurement
% errors as well (as long as you simulate them in the measurements as well!

S = diag([var_n_t0 var_n_t1 var_n_t2]);

H = [1 (t0-tb); 
     1 (t1-tb); 
     1 (t2-tb)];


 mu = 0;
 MSE = 0;
 MSE_n = 0;
 
 mu_n = 0;
 
 Nexp = 10000;
 for k = 1:Nexp
 
% Simulate the measurements  
 

n = sqrt(S)*randn(3,1);

z = H*x + n;

%
% Calcualte the Fisher matrix, the ML-estimate of x and the Cramer-Rao
% lower bound. 
% 

F = H'*inv(S)*H;

x_ML = inv(F)*H'*inv(S)*z;

CR_bound = inv(F);

%
% From the CR lower bound you can compute x-sigma confidence intervals
% where the true x should be in with certain probability... Check if this 
% is the case! 
%

mu = mu + (x-x_ML);
MSE = MSE + (x-x_ML)*(x-x_ML)';

mu_n = mu_n + n;
MSE_n = MSE_n + n*n';

 end
CR_bound = inv(F)

mu = mu/Nexp
MSE = MSE/Nexp

MSE_n = MSE_n/Nexp
 
 ThreeSigma_bounds = 3*sqrt(diag(CR_bound));