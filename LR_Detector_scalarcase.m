%
% This program simulates a simple scalar detection problem (see Lecture Notes SSDP 2-2 
%
clear;

P_H1 = 0.8;
P_H0 = 1-P_H1;

s = 10;

sigma_n = 5;

n = sigma_n*randn;

tau = 2;

H1 = rand < P_H1

if (H1)
    
    z = s + n
    
else 
    
    z = n
    
end
    
detect = z > (s/2 + sigma_n^2*log(tau)/s)