function [x_new, P_new] = kf_predict(x, P, F, G, noise_variance)
    x_new = F * x;
    P_new = F * P * F' + G*noise_variance;
end