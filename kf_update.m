function [x_new, P_new] = kf_update(NUM_OF_VAR, meas_value, meas_variance, x, P, H)
    z = meas_value;
    R = meas_variance;
    
    y = z - H * x;
    S = H*P*H' + R;
    
    K = P*H'*inv(S);
    
    x_new = x + K*y;
    P_new = (eye(NUM_OF_VAR) - K*H)*P;
    
end