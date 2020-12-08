clear;
close all;

%% Initialization 

L = 100;
L_meas = 1;

dt = 0.1;

NUM_OF_VAR = 3;

t = 0:dt:L*dt; % time scale
t_meas = 0:dt*L_meas:dt*L;

meas_variance = 10^2;
noise_variance = 100^2;

%% initial conditions

a_init = 3;
v_init = 1;
x_init = 0;


%% true values

a = 3;

bin = zeros(1, length(t));
bin(round(L/2):end) = -10;

bin(1:round(L/2)-1) = 5;

a_true = a .* bin;

v_true_0 = 0.5;

v_true = v_true_0 + a_true .* t;

x_true_0 = 0;
x_true = x_true_0 + v_true .* t + 0.5 * t.^2 .* a_true;


%% Available measurements

xm = x_true + randn(size(t)) .* sqrt(meas_variance);

%% State space matrices

F = [1 dt 0.5*dt^2; 0 1 dt; 0 0 1];

G = [0.5*dt^2; dt; 1];

H = [1 0 0];

mu_0 = [x_init; v_init; a_init];
cov_0 = diag([1 1 1]);

mus = zeros(length(mu_0), length(t));
covs = zeros(size(cov_0, 1), size(cov_0, 2), length(t));

mus(:, 1) = mu_0;
covs(:, :, 1) = cov_0;


%% Kalman filter

for steps = 1:L
    
%% Predict
    [x_new, P_new] = kf_predict(mus(:, steps), (covs(:, :, steps)), F, G, noise_variance);
    
    if (steps ~= 0) && (mod(steps, L_meas) == 0)
       [x_new, P_new] = kf_update(NUM_OF_VAR, xm(steps), meas_variance, x_new, P_new, H);
    end
    
    mus(:, steps+1) = x_new;
    covs(:, :, steps+1) = P_new;
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

plot(t, v_true, 'g');
hold on;
% plot(t, xm, 'r');
hold on;
plot(t, mus(2, :), 'b');
hold on;
plot(t, mus(2, :) - 2 .* squeeze(sqrt(covs(2, 2, :))).', 'k-.');
hold on;
plot(t, mus(2, :) + 2 .* squeeze(sqrt(covs(2, 2, :))).', 'k-.');


figure;

plot(t, a_true, 'g');
hold on;
% plot(t, xm, 'r');
hold on;
plot(t, mus(3, :), 'b');
hold on;
plot(t, mus(3, :) - 2 .* squeeze(sqrt(covs(3, 3, :))).', 'k-.');
hold on;
plot(t, mus(3, :) + 2 .* squeeze(sqrt(covs(3, 3, :))).', 'k-.');


%% error visualization 



