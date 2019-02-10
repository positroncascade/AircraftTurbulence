aircraft = AircraftModel();

%% full state space analysis
system = aircraft.state_space();
% damp(system)
K = [0, -0.05, 0, 0; 0, 0, 0, 0];
% damp(system.A - system.B * K)
system.A = system.A - system.B * K;


%% Reduced state space analysis phi = 0 and beta = 0;
system_red = aircraft.simple_state_space();
damp(system_red)
K = [0, -1; 0, 0];
damp(system_red.A - system_red.B * K)

%% augmented
system_aug = aircraft.augmented_state_space();
K = zeros(5, 10);
K(1, 2) = -0.05;
system_aug.A = system_aug.A - system_aug.B * K;
damp(system_aug)

