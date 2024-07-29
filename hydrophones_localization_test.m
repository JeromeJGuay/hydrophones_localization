clear all
close all
clc

%% Known parameters

% sound-speed [m/s]
c0 = 1500;

% Sources Positions
xs = [477.613 414.672 340.623 550.427 494.891 381.35 ];
ys = [164.939 283.547 202.005 255.749 346.558 138.994];
zs = [5 5 5 5 5 5];

sources = [xs(:).';ys(:).';zs(:).'];
Ns = size(sources, 2);

% pyramid position
pyramid_x = 444;
pyramid_y = 233;
pyramid_z =  30;

% receiver locations
r_1_xyz = [-1  -1    0];
r_2_xyz = [ 0  .5    0];
r_3_xyz = [1   -1   0];
r_4_xyz = [0  -.5  -1];

receivers = ([r_1_xyz; r_2_xyz; r_3_xyz; r_4_xyz] + [pyramid_x pyramid_y pyramid_z])';
Nr = size(receivers, 2);

% synthetic travel time and delays
xr = receivers(1,:);
yr = receivers(2,:);
zr = receivers(3,:);

tau_meas = sqrt( (xr(:)-xs(:).').^2  + (yr(:)-ys(:).').^2 + (zr(:)-zs(:).').^2 ) / c0;

delays_meas = NaN(Nr,Nr,Ns);
for i_e = 1:Ns
    delays_meas(:,:,i_e) = tau_meas(:,i_e) - tau_meas(:,i_e).';
end

%% Variables Initialization
% receiver guess location

r0_1_xyz = [-1  -1    0];
r0_2_xyz = [ 0  .5    0];
r0_3_xyz = [1   -1   0];
r0_4_xyz = [0  -.5  -1];

r0 = ([r0_1_xyz; r0_2_xyz; r0_3_xyz; r0_4_xyz] + [pyramid_x pyramid_y pyramid_z])';

%% Solving
% None Linear constrainte (distance between hydrophones)
%max_distance = 3;
%nonlcon = @(X) nonlinear_constraint.max_distance(X, max_distance);

cfun = @(X) lowcost_functions.R(X, sources, c0, delays_meas);

% z_range = -20:20;
delta_zr0 = -10:0.5:10;

mean_positions = zeros(size(delta_zr0, 2),3);

for i=1:length(delta_zr0)
    r00 = r0 + [1;1;1] * delta_zr0(i);
    xsol = fmincon(cfun,r00);
    mean_positions(i,:) = mean(xsol, 2);
end

%% Figures
figure
plot(delta_zr0, mean_positions(:,1) - pyramid_x)
xlabel('\Delta Zr_0')
ylabel('solution x diff')

figure
plot(delta_zr0, mean_positions(:,2) - pyramid_y)
xlabel('\Delta Zr_0')
ylabel('solution y diff')

figure
plot(delta_zr0, mean_positions(:,3) - pyramid_z)
xlabel('\Delta Zr_0')
ylabel('\Delta Z')
%figures.pyramid_3d(r0, receivers, xsol)


%%
tools.compute_distances_receivers_to_sources