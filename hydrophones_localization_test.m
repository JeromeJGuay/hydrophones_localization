clear all
close all
clc

%% Known parameters

% sound-speed [m/s]
sound_speed = 1500;

sources = [
    477.613 164.939 5;
    414.672 283.547 5;
    340.623 202.005 5;
    550.427 255.749 5;
    494.891 346.558 5;
    381.35  138.994 5;
    ];


pyramid = [444 233 30];

%% Synthetic Data
receivers = [
    -1.0  -1.0   0.0; % receiver 1
     0.0   0.5   0.0; % receiver 2
     1.0  -1.0   0.0; % receiver 3
     0.0  -0.5  -1.0; % receiver 4
    ] + pyramid;


taus_meas = tools.compute_taus_receivers_to_sources(receivers,sources,sound_speed);
delays_meas = tools.compute_receivers_delays(taus_meas);

%% add noise to the delays
delays_meas = delays_meas + 0.0001 * (2*rand(size(delays_meas))-1);


%% Variables Initialization
receivers_0 = [
    -1.0  -1.0   0.0; % receiver 1
     0.0   0.5   0.0; % receiver 2
     1.0  -1.0   0.0; % receiver 3
     0.0  -0.5  -1.0; % receiver 4
    ] + pyramid;

%% Add noise to initial guess geometry
% receivers_0 = receivers_0 + 0.1 * (2*rand(size(receivers_0))-1);

%% add position x,y offset
% receivers_0 = receivers_0 + [5, 0, 0];

%% Solving
% None Linear constrainte (distance between hydrophones)
%max_distance = 3;
%nonlcon = @(X) nonlinear_constraint.max_distance(X, max_distance);

cfun = @(X) lowcost_functions.R(X, sources, sound_speed, delays_meas);

delta_zr0 = -1:.25:1;

xsols = zeros(horzcat(length(delta_zr0), size(receivers)));


for i=1:length(delta_zr0)
    r0 = receivers_0 + [0 0 1] * delta_zr0(i);
    xsol = fmincon(cfun, r0);
    xsols(i,:,:) = xsol;
end



%% Figures
receivers_expanded = permute( ...
    repmat(receivers, [1, 1, length(xsols)]), ...
    [3 1 2]); 
delta_r = sqrt(sum((xsols - receivers_expanded).^2, 3));

figure()
for i_r=1:4
    subplot(2,2,i_r)
    plot(delta_zr0, delta_r(:, i_r))
    title(sprintf('Receiver %d', i_r))
    xlabel("\Delta Z r0")
    ylabel("Position Error")
end

%%
figures.pyramid_3d_solution(receivers, receivers_0, xsol)
