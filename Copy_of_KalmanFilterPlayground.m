% This file simulates a Kalman Filter on a 
% 1D constant velocity linear Gaussian estimation problem. 

% Initial position and speed of an object
clear;
x = 0;
v = 50;

% Measurement Interval

T = .1;

% Observation and process noise variance

varn = 100^2;

varw = 10^2;

K = 100

t = 0:T:K*T;

% Generate a set of position measurements

x_true = x+t*v;
xm = x_true+sqrt(varn)*randn(size(t));

%xm = x_true+sqrt(varn*x_true).*randn(size(t));



% Initial state of the Kalman filter 

sp = [100; 0];
Qp = diag([1E6 1E4]);

% Definition of system matrices: 
% H = Observation matrix, 
% F = System Dynamics Matrix
% W = part of Process Noise Covariance matrix

H = [1 0];

F = [1 T; 0 1];
%C = [T^2/2; T];
W = [T^3/3 T^2/2; T^2/2 T];

for k=1:K
    
    k
    % KF update

    G = Qp*H'*inv(H*Qp*H' + varn);
%    G = Qp*H'*inv(H*Qp*H' + x_true(k)*varn);
    
    sf = sp + G*(xm(k)-H*sp);
    Qf = (eye(2)-G*H)*Qp;
    

   mean_KF(:,k) = sf;
   
   sigmas_KF(:,k) = sqrt(diag(Qf));
   
   figure(1)
   plot ((1:k)*T, mean_KF(1,1:k)-x_true(1:k), 'r+')
   hold 
   plot ((1:k)*T, xm(1:k)-x_true(1:k), 'b+')
   plot ((1:k)*T, 3*sigmas_KF(1,1:k), 'r-.')
   plot ((1:k)*T, -3*sigmas_KF(1,1:k), 'r-.')
   hold
   title ('Position estimation error of KF')
   legend ('KF-error','measurements','3-sigma CR-bound')
   
   axis ([0 K*T -3*sigmas_KF(1,1) 3*sigmas_KF(1,1)])
   
   figure(2)
   plot ((1:k)*T, mean_KF(2,1:k)-v, 'r+')
   hold 
   plot ((1:k)*T, 3*sigmas_KF(2,1:k), 'r-.')
   plot ((1:k)*T, -3*sigmas_KF(2,1:k), 'r-.')
   hold
   title ('Velocity estimation error of KF')
   legend ('KF-error', '3-sigma CR-bound')
   
   axis ([0 K*T -3*sigmas_KF(2,1) 3*sigmas_KF(2,1)])
   

    
    % KF prediction
    
    sp = F*sf;
    Qp = F*Qf*F' + W*varw;
     
    
    
end

figure(10)

plot (x_true)
hold
plot (xm, 'x')
plot (mean_KF(1,:), 'o')

hold
