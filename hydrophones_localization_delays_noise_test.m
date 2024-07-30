clear all
close all
clc
%% Switch

ADD_DELAYS_NOISE = false;
ADD_GEOMETRY_NOISE = false;
ADD_SOURCES_NOISE = false;
SOUND_SPEED_0_OFFSET = 0; % m/s

RECEIVERS_0_OFFSET = [0 0 0]; % m

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


pyramid = [444 233 100];

%% Synthetic Data
receivers = [
    -1.0  -1.0   0.0; % receiver 1
     0.0   0.5   0.0; % receiver 2
     1.0  -1.0   0.0; % receiver 3
     0.0  -0.5  -1.0; % receiver 4
    ] + pyramid;


taus_meas = pdist2(receivers, sources) / sound_speed;
delays_meas = tools.compute_receivers_delays(taus_meas);

%% add noise to the delays
if ADD_DELAYS_NOISE == true
    delays_meas = delays_meas + 0.0001 * (2*rand(size(delays_meas))-1);
end

%% Variables Initialization
receivers_0 = [
    -1.0  -1.0   0.0; % receiver 1
     0.0   0.5   0.0; % receiver 2
     1.0  -1.0   0.0; % receiver 3
     0.0  -0.5  -1.0; % receiver 4
    ] + pyramid;


sound_speed_0 = sound_speed + SOUND_SPEED_0_OFFSET;

%% Add noise to initial guess geometry
if ADD_GEOMETRY_NOISE == true
    receivers_0 = receivers_0 + 0.05 * (2*rand(size(receivers_0))-1);
end
%% add position x,y offset

receivers_0 = receivers_0 + RECEIVERS_0_OFFSET;

%% Solving for a range of vertical r0 position.

fmincon_options = optimoptions('fmincon','display','off');
delta_zr0 = 0:1:10;
    
xsols = zeros(horzcat(length(delta_zr0), size(receivers)));

cfun = @(X) lowcost_functions.R(X, sources, sound_speed_0, delays_meas);
for i=1:length(delta_zr0)
    r0 = receivers_0 + [0 0 1] * delta_zr0(i);
    
    X0 = r0(:); % vector form
    
    X1 = fmincon(cfun, X0,  [],[],[],[],[],[],[], fmincon_options); % vector form
    
    xsol = reshape(X1, length(X0)/3, 3);
    
    xsols(i,:,:) = xsol;
end




%% Figures Receivers position error
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
    ylabel("Position Error (m)")
    xlim([delta_zr0(1) delta_zr0(end)])
end


%% Figures distance receivers-receivers error
ddr = pdist2(receivers, receivers);
ddr_expanded = permute( ...
    repmat(ddr, [1, 1, length(xsols)]), ...
    [3 1 2]); 

ddr0 = zeros(length(delta_zr0), length(receivers), length(receivers));
for i=1:length(delta_zr0)
    ddr0(i,:,:) = pdist2(squeeze(xsols(i,:,:)), squeeze(xsols(i,:,:)));
end

delta_ddr = sqrt((ddr_expanded - ddr0).^2);
index_ij = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
subtitles = ['1-2'; '1-3'; '1-4'; '2-3'; '2-4'; '3-4'];
figure()
for i_rr=1:length(index_ij)
    i=index_ij(i_rr,1);
    j=index_ij(i_rr,2);
    subplot(2,3,i_rr)
    plot(delta_zr0, delta_ddr(:, i, j))
    title(subtitles(i_rr,:))
    xlabel("\Delta Z r0")
    ylabel("Distance Error (m)")
    xlim([delta_zr0(1) delta_zr0(end)])
end

%% Figure 3D solutions
figures.pyramid_3d_solution(receivers, receivers_0, xsol)
