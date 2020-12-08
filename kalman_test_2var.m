clear;
close all;

%% Initialization 

L = 100;
L_meas = 1;

dt = 0.1;

NUM_OF_VAR = 2;

t = 0:dt:L*dt; % time scale
t_meas = 0:dt*L_meas:dt*L;

meas_variance = 10^2;
noise_variance = 100^2;

%% initial conditions

% a_init = 3;
v_init = 1;
x_init = 0;


%% true values

% a_true = 3;
v_true_0 = 50;
v_true = v_true_0;

x_true_0 = 0;
x_true = x_true_0 + v_true .* t;


%% Available measurements

xm = x_true + randn(size(t)) .* sqrt(meas_variance);

%% State space matrices

F = [1 dt; 0 1];

% G = [0.5*dt^2; dt];

G = [dt^3/3 dt^2/2; dt^2/2 dt];

H = [1 0];

mu_0 = [x_init; v_init];
cov_0 = diag([1 1]);

mus = zeros(length(mu_0), length(t));
covs = zeros(size(cov_0, 1), size(cov_0, 2), length(t));

mus(:, 1) = mu_0;
covs(:, :, 1) = diag([1 1]);

%% Kalman filter

for steps = 1:L
    
%% Predict
    [x_new, P_new] = kf_predict(mus(:, steps), covs(:, :, steps), F, G, noise_variance);
    disp('prediction')
%     disp(x_new)
    disp(P_new)
   
    
    if (steps ~= 0) && (mod(steps, L_meas) == 0)
       [x_new, P_new] = kf_update(NUM_OF_VAR, xm(steps), meas_variance, x_new, P_new, H);
       
       disp('update')
       disp(P_new)
    end
    
    
    mus(:, steps+1) = x_new;
    covs(:, :, steps+1) = P_new;
    
    
%    figure(1)
%    plot ((1:steps)*dt, mus(1,1:steps)-x_true(1:steps), 'r+')
%    hold 
%    plot ((1:steps)*dt, xm(1:steps)-x_true(1:steps), 'b+')
%    plot ((1:steps)*dt, 3*covs(1,1:steps), 'r-.')
%    plot ((1:steps)*dt, -3*covs(1,1:steps), 'r-.')
%    hold
%    title ('Position estimation error of KF')
%    legend ('KF-error','measurements','3-sigma CR-bound')
%    
% %    axis ([0 L*dt -3*covs(1,1) 3*covs(1,1)])
%    
%    figure(2)
%    plot ((1:steps)*dt, mus(2,1:steps)-v_true, 'r+')
%    hold
%    plot ((1:steps)*dt, 3*covs(1,1:steps), 'r-.')
%    plot ((1:steps)*dt, -3*covs(1,1:steps), 'r-.')
%    hold
%    title ('Velocity estimation error of KF')
%    legend ('KF-error', '3-sigma CR-bound')
%    
% %    axis ([0 L*dt -3*covs(2,1) 3*covs(2,1)])
   
end

%% Plot

figure;

plot(t, x_true, 'g');
hold on;
plot(t, xm, 'r');
hold on;
plot(t, mus(1, :), 'b');
hold on;
plot(t, mus(1, :) - 2 .* squeeze(sqrt(covs(1, 1, :))).', 'k-.');
hold on;
plot(t, mus(1, :) + 2 .* squeeze(sqrt(covs(1, 1, :))).', 'k-.');

figure;

plot(t, v_true .* ones(length(t)), 'g');
hold on;
plot(t, mus(2, :), 'b');
hold on;
plot(t, mus(2, :) - 2 .* squeeze(sqrt(covs(2, 2, :))).', 'k-.');
hold on;
plot(t, mus(2, :) + 2 .* squeeze(sqrt(covs(2, 2, :))).', 'k-.');


% figure;
% 
% plot(t, a_true .* ones(length(t)), 'g');
% hold on;
% % plot(t, xm, 'r');
% hold on;
% plot(t, mus(3, :), 'b');
% hold on;
% plot(t, mus(3, :) - 2 .* sqrt(covs(3, :)), 'k-.');
% hold on;
% plot(t, mus(3, :) + 2 .* sqrt(covs(3, :)), 'k-.');


