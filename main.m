rng(284);
aircraft = AircraftModel();

% full state space analysis
system = aircraft.state_space();
damp(system)
% damp(system)
K = [0, -0.05, 0, 0; 0, 0, 0, 0];
% damp(system.A - system.B * K)
system.A = system.A - system.B * K;


% Reduced state space analysis phi = 0 and beta = 0;
system_simple = aircraft.augmented_simple_state_space();
damp(system_simple)
K = zeros(5, 8); K(1, 2) = -1;
damp(system_simple.A - system_simple.B * K)

% augmented
system_aug = aircraft.augmented_state_space();
K = zeros(5, 10);
K(1, 2) = -0.05;
system_aug.A = system_aug.A - system_aug.B * K;
damp(system_aug)


%% Time domain simulation 
plot_sim = 0;
% add ay to complete system
% Add a_y to the state space: a_y = V * (beta_dot + psi_dot)
psi_c = zeros(1, 10); psi_c(1, 4) = 2 * aircraft.V / aircraft.b;
C = [system_aug.C; aircraft.V * (system_aug.A(1, :) + psi_c)];
D = [system_aug.D; zeros(1, 5)];
system_aug_ay = ss(system_aug.A, system_aug.B, C, D);

% add ay to simple sys
% psi_dot = r and beta = -psi -> a_y = 0;
a_y_c = zeros(1, 8); a_y_c(1, 2) = 0;
C = [system_simple.C; a_y_c];
system_simple_ay = ss(system_simple.A, system_simple.B, C, 0);

% simulate
time_max =  1500; 
dt = 0.004; 
[y, t, w] = time_domain_sim(system_aug_ay, dt, time_max);
N = length(t);
[y_simple, t, w] = time_domain_sim(system_simple_ay, dt, time_max);
%plotting
if plot_sim
    % PLOT RESULTS
    beta_axis = [0 60 -0.07  0.07];
    phi_axis  = [0 60 -0.15  0.15];
    pb_axis   = [0 60 -1e-2  1e-2];
    rb_axis   = [0 60 -1e-2  1e-2];
    ay_axis = [0 60 -0.6 0.6];

    fig0 = figure;
    plot(t, w) 
    xlabel('time, s'); ylabel('White Gaussian Noise');
    saveas(fig0, 'images/noise.png');

    % RESPONSE TO v_g
    fig1 = figure;
    subplot(2,1,1); plot(t,y(:,1)); axis(beta_axis); 
    xlabel('time, s'); ylabel('\beta [rad]');
    subplot(2,1,2); plot(t,y(:,2)); axis(phi_axis);
    xlabel('time, s'); ylabel('\phi [rad]');
    saveas(fig1, 'images/betaphi.png');

    fig2 = figure;
    subplot(2,1,1); plot(t,y(:,3)); axis(pb_axis); 
    xlabel('time, s'); ylabel('pb/2V');
    subplot(2,1,2); plot(t,y(:,4)); axis(rb_axis);
    xlabel('time, s'); ylabel('rb/2V');
    saveas(fig2, 'images/rp.png');
    %print -depsc2 -r1200 fig8_16b2
    fig3 = figure;
    plot(t, y(:, 5));
    xlabel('time, s'); ylabel('a_y [m/s^2]'); axis(ay_axis);
    saveas(fig3, 'images/ay.png');

    %% Plot result simplified model
    fig4 = figure;
    subplot(2,1,1); plot(t,y_simple(:,1)); axis(beta_axis); 
    xlabel('time, s'); ylabel('\psi [rad]');
    subplot(2,1,2); plot(t,y_simple(:,2)); axis(rb_axis);
    xlabel('time, s'); ylabel('rb/2V');
    saveas(fig4, 'images/rpsi_red.png');

    fig5 = figure;
    plot(t, y_simple(:, 3));
    xlabel('time, s'); ylabel('a_y [m/s^2]'); axis(ay_axis)
    saveas(fig5, 'images/ay_red.png');
end

%% spectral analysis
% analytical
w = logspace(-2, 2, 300);

spectrum_ana = analitycal_spectrum(system_aug_ay, w, 5);
spectrum_red_ana = analitycal_spectrum(system_simple_ay, w, 5);

% experimental

% complete model
beta = y(:, 1);
phi = y(:, 2);
p = y(:, 3);
r = y(:, 4);
a_y = y(:, 5);

beta_w = fft(beta);
phi_w = fft(phi);
p_w = fft(p);
r_w = fft(r);
a_y_w = fft(a_y);

S_beta = dt/N * beta_w .* conj(beta_w);
S_beta_smooth = smoothing_filter(S_beta);
S_phi = dt/N * phi_w .* conj(phi_w);
S_phi_smooth = smoothing_filter(S_phi);
S_p = dt/N * p_w .* conj(p_w);
S_p_smooth = smoothing_filter(S_p);
S_r = dt/N * r_w .* conj(r_w);
S_r_smooth = smoothing_filter(S_r);
S_a_y = dt/N * a_y_w .* conj(a_y_w);
S_a_y_smooth = smoothing_filter(S_a_y);

% reduced model
psi_red = y_simple(:, 1);
r_red = y_simple(:, 2);
a_y_red = y_simple(:, 3);

psi_red_w = fft(psi_red);
r_red_w = fft(r_red);
a_y_red_w = fft(a_y_red);

S_psi_red = dt/N * psi_red_w .* conj(psi_red_w);
S_psi_red_smooth = smoothing_filter(S_psi_red);
S_r_red = dt/N * r_red_w .* conj(r_red_w);
S_r_red_smooth = smoothing_filter(S_r_red);
S_a_y_red = dt/N * a_y_red_w .* conj(a_y_red_w);
S_a_y_red_smooth = smoothing_filter(S_a_y_red);

% frequency axis
omega_0 = 2 * pi/time_max;
omega = omega_0 * (0:1:N/2 - 1);

plot_spectra = 1;
% PLOT ANALYTIC AND ESTIMATED PSDS IN ONE PLOT
if plot_spectra
    clf
    subplot(2,1,1); loglog(w,spectrum_ana(1, :),'--',omega,S_beta(1:N/2), omega, S_beta_smooth(1:N/2), 'b'); 
    axis(10.^[-2 2 -12 -2]); xlabel('\omega [rad/s]'); ylabel('S_{\beta}');
    legend('Analytical', 'Experimental', 'Smoothed');
    subplot(2,1,2); loglog(w,spectrum_ana(2, :),'--',omega,S_phi(1:N/2), omega, S_phi_smooth(1:N/2), 'b');
    axis(10.^[-2 2 -12 -2]); xlabel('\omega [rad/s]'); ylabel('S_{\phi}')
    legend('Analytical', 'Experimental', 'Smoothed');
    saveas(gcf, 'images/spectra_states_0.png');

    fig0bis = figure;
    subplot(2,1,1); loglog(w,spectrum_ana(3, :),'--',omega,S_p(1:N/2), omega, S_p_smooth(1:N/2), 'b');
    axis(10.^[-2 2 -14 -2]); xlabel('\omega [rad/s]'); ylabel('S_{pp}')
    legend('Analytical', 'Experimental', 'Smoothed');
    subplot(2,1,2); loglog(w,spectrum_ana(4, :),'--',omega,S_r(1:N/2), omega, S_r_smooth(1:N/2), 'b');
    axis(10.^[-2 2 -14 -2]); xlabel('\omega [rad/s]'); ylabel('S_{rr}')
    legend('Analytical', 'Experimental', 'Smoothed');
    saveas(gcf, 'images/spectra_states_1.png');

    fig1 = figure; loglog(w,spectrum_ana(5, :),'--',omega,S_a_y(1:N/2), omega, S_a_y_smooth(1:N/2), 'b');
    xlabel('\omega [rad/s]'); ylabel('S_{a_y}')
    legend('Analytical', 'Experimental', 'Smoothed');
    saveas(fig1, 'images/spectra_ay.png');

    fig2 = figure;
    subplot(2,1,1); loglog(w, spectrum_red_ana(1, :), '--', omega, S_psi_red(1: N/2), omega, S_psi_red_smooth(1:N/2), 'b');
    axis(10.^[-2 2 -14 -2]); xlabel('\omega [rad/s]'); ylabel('S_{\psi}');
    subplot(2,1,2); loglog(w, spectrum_red_ana(2, :), '--', omega, S_r_red(1: N/2), omega, S_r_red_smooth(1:N/2), 'b');
    axis(10.^[-2 2 -14 -2]); xlabel('\omega [rad/s]'); ylabel('S_{rr}');
    legend('Analytical', 'Experimental', 'Smoothed');
    saveas(fig2, 'images/spectra_states_red.png');
end

%% Variances
% analytical 
% complete model
var_beta_ana = 1/pi * trapz(w, spectrum_ana(1, :));
var_phi_ana = 1/pi * trapz(w, spectrum_ana(2, :));
var_p_ana = 1/pi * trapz(w, spectrum_ana(3, :));
var_r_ana = 1/pi * trapz(w, spectrum_ana(4, :));
var_ay_ana = 1/pi * trapz(w, spectrum_ana(5, :));

% reduced model
var_psi_red_ana = 1/pi * trapz(w, spectrum_red_ana(1, :));
var_r_red_ana = 1/pi * trapz(w, spectrum_red_ana(2, :));
var_ay_red_ana = 1/pi * trapz(w, spectrum_red_ana(3, :));

% var on experimental pds
% complete
var_beta_exp = 1/pi * trapz(omega, S_beta(1:N/2));
var_phi_exp = 1/pi * trapz(omega, S_phi(1:N/2));
var_p_exp = 1/pi * trapz(omega, S_p(1:N/2));
var_r_exp = 1/pi * trapz(omega, S_r(1:N/2));
var_ay_exp = 1/pi * trapz(omega, S_a_y(1:N/2));

var_beta_exp_smooth = 1/pi * trapz(omega, S_beta_smooth(1:N/2));
var_phi_exp_smooth = 1/pi * trapz(omega, S_phi_smooth(1:N/2));
var_p_exp_smooth = 1/pi * trapz(omega, S_p_smooth(1:N/2));
var_r_exp_smooth = 1/pi * trapz(omega, S_r_smooth(1:N/2));
var_ay_exp_smooth = 1/pi * trapz(omega, S_a_y_smooth(1:N/2));

% reduced
var_psi_red_exp = 1/pi * trapz(omega, S_psi_red(1:N/2));
var_r_red_exp = 1/pi * trapz(omega, S_r_red(1:N/2));
var_red_exp = 1/pi * trapz(omega, S_a_y_red(1:N/2));

var_psi_red_exp_smooth = 1/pi * trapz(omega, S_psi_red_smooth(1:N/2));
var_r_red_exp_smooth = 1/pi * trapz(omega, S_r_red_smooth(1:N/2));
var_red_exp_smooth = 1/pi * trapz(omega, S_a_y_red_smooth(1:N/2));

% var from time series
var_beta_time = var(beta);
var_phi_time = var(phi);
var_p_time = var(p);
var_r_time = var(r);
var_a_y_time = var(a_y);

% reduced
var_psi_red_time = var(psi_red);
var_r_red_time = var(r_red);
var_a_y_red_time = var(a_y_red);

% find ensemble variance
num_ensembles = 100;
vars = [0, 0, 0, 0, 0];
vars_red = [0, 0, 0];
for ensemble = 1:num_ensembles
    time_max = 60; 
    dt = 0.04; 
    [y, t, w] = time_domain_sim(system_aug_ay, dt, time_max);
    [y_simple, t, w] = time_domain_sim(system_simple_ay, dt, time_max);
    vars = vars + var(y);
    vars_red = vars_red + var(y_simple);
end
vars = vars/num_ensembles;
vars_red = vars_red/num_ensembles;


