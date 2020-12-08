% This file simulates a Kalman Filter on a 
% 1D constant velocity linear Gaussian estimation problem. 

% Initial position and speed of an object
clear;
close all;


% Measurement Interval

T = .1;

% Observation and process noise variance

varn = 0.1^2;

varw = 0.1;

L = 100
L_ = 20

x_true = zeros(1, L);
% x = [0 0 0.1 0.1];
v_true = 0.5 .* ones(1, L);
a = 3;


t = 0:T:L*T;

% Generate a set of position measurements



% xm = x_true+sqrt(varn*x_true).*randn(size(t));
for i = 2:L
    v_true(i+1) = v_true(i) + a .* T;
    x_true(i+1) = x_true(i) + v_true(i) .* T + T.^2/2*a;
end

vm = v_true + sqrt(varn)*randn(size(t));


xm = x_true +sqrt(varn)*randn(size(t));

a_true = a * ones(size(t));
am = a_true + sqrt(varn)*randn(size(t));

% xm_ = [xm; vm; am];

% Initial state of the Kalman filter 

sf = [0; 5; 1];
Qf = diag([1 1 1]);

% Definition of system matrices: 
% H = Observation matrix, 
% F = System Dynamics Matrix
% W = part of Process Noise Covariance matrix

% H = [1 0];
H = [1 0 0];

% F = [1 T; 0 1];
F = [1 T T^2/2; 0 1 T; 0 0 1];
%C = [T^2/2; T];

%W = [T^3/3 T^2/2; T^2/2 T];

W = [T^3/3 T^2/2 T]*[T^3/3 T^2/2 T].';

for k=1:L
    
    % KF prediction
    
    sp = F*sf;
    Qp = F*Qf*F' + W*varw;
    
    k
    % KF update
    
    if (k ~= 0) &&  mod(k, L_) == 0 
        
%      xm_ = x_true(k)+sqrt(varn).*randn();
     y = xm(k) - H * sp;
     S = H * Qp * H' + varn;
     K = Qp * H'/S;
     
%      G = Qp*H'*inv(H*Qp*H' + varn);
%    G = Qp*H'*inv(H*Qp*H' + x_true(k)*varn);
    
%      sf = sp + G.*(xm(k)-H*sp);
     sf = sp + K * y;
     
%     sfv = sp + G*(vm(k)-H*sp);
%      Qf = (eye(3)-G*H)*Qp;
     Qf = (eye(3) - K * H) * Qp;
    end
    

%    mean_KF(:,k) = sf;
   mean_KF(:, k) = sf;
   
%    sigmas_KF(:,k) = sqrt(diag(Qf));
   sigmas_KF(:, k) = sqrt(diag(Qf));
   
   figure(1)
   plot ((1:k)*T, mean_KF(1,1:k)-x_true(1:k), 'r+')
   hold 
   plot ((1:k)*T, xm(1:k)-x_true(1:k), 'b+')
   plot ((1:k)*T, 3*sigmas_KF(1,1:k), 'r-.')
   plot ((1:k)*T, -3*sigmas_KF(1,1:k), 'r-.')
   hold
   title ('Position estimation error of KF')
   legend ('KF-error','measurements','3-sigma CR-bound')
   
   axis ([0 L*T -3*sigmas_KF(1,1) 3*sigmas_KF(1,1)])
   
   figure(2)
   plot ((1:k)*T, mean_KF(2,1:k)-v_true(1:k), 'r+')
   hold 
   plot ((1:k)*T, vm(1:k)-v_true(1:k), 'b+')
   plot ((1:k)*T, 3*sigmas_KF(2,1:k), 'r-.')
   plot ((1:k)*T, -3*sigmas_KF(2,1:k), 'r-.')
   hold
   title ('Velocity estimation error of KF')
   legend ('KF-error', 'measurements', '3-sigma CR-bound')
   
   axis ([0 L*T -3*sigmas_KF(2,1) 3*sigmas_KF(2,1)])
   

    
  
     
    
%     pause
end

figure(10)

plot (t, x_true)
hold
plot (t, xm, 'x')
plot (t(1:end-1), mean_KF(1,:), 'o')

hold

figure(11)

plot (t, v_true);
hold
plot (t, vm, 'x')
plot (t(1:end-1), mean_KF(2,:), 'o')

hold

figure(12)

plot (t, a_true);
hold
plot (t, am, 'x')
plot (t(1:end-1), mean_KF(3,:), 'o')

hold
