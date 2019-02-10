function [y, t, w] = time_domain_sim(system, dt, time_max)
%% simulation
% TIME AXIS AND INPUT VECTOR DEFINITION
    T  = time_max; t = [0:dt:T]; N = length(t);
    nn = zeros(1,N);

    % TURBULENCE INPUTS
    v_g = randn(1,N)/sqrt(dt);

    % INPUT VECTORS    
    w = [nn' nn' nn' nn' v_g'];
    %% From example 8.1

    % RESPONSE to u_g
    y = lsim(system, w, t);
end



