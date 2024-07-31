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

receivers = [
    -1.0  -1.0   0.0; % receiver 1
     0.0   0.5   0.0; % receiver 2
     1.0  -1.0   0.0; % receiver 3
     0.0  -0.5  -1.0; % receiver 4
    ] + pyramid;

%% Synthetic Data
taus_measured = pdist2(receivers, sources) / sound_speed;
delays_measured = tools.compute_receivers_delays(taus_measured);

%% Variables Initialization
receivers_0 = [
    -1.0  -1.0   0.0; % receiver 1
     0.0   0.5   0.0; % receiver 2
     1.0  -1.0   0.0; % receiver 3
     0.0  -0.5  -1.0; % receiver 4
    ] + pyramid;


%% Solving for a range of horizontal position.

DELTA_THETA = linspace(0, 2*pi, 12);
DELTA_RADIUS = 1e-2*10.^(0:.1:4);

rt_sols = zeros( ...
    horzcat( ...
    length(DELTA_RADIUS), length(DELTA_THETA), size(receivers) ...
    ) ...
    );


r0 = receivers_0;
dtaus0 = delays_measured;
sources_0 = sources;
sound_speed_0 = sound_speed;

for ri=1:length(DELTA_RADIUS)
    fprintf('Iteration Radial: %d\n', ri)
    for ti=1:length(DELTA_THETA)
    
        [dx, dy] = pol2cart(DELTA_THETA(ti), DELTA_RADIUS(ri));
            
        r0 = receivers_0 + [dx dy 0];
      
        rsol = solver.R(r0, sources_0, sound_speed_0, dtaus0);
       
        rt_sols(ri,ti,:,:) = rsol;
    end
end


%%  Receivers position error

%permute(permute(xsols, [2 3 1]) - receivers, [3 1 2]);

delta_r = sqrt(sum(( ...
    permute(permute(rt_sols, [3 4 1 2]) - receivers, [3 4 1 2]) ...
    ).^2, 4));

mean_delta_r = squeeze(mean(delta_r, [3]));

%% Figures
figure()
boxplot(mean_delta_r.', DELTA_RADIUS)

ax = gca;
ax.YAxis.Scale ="log";

xlabel("\Delta Radius (m) ")
ylabel("Position Error (m)")
title("Error on guess horizontal position")
grid on

%% Figures
[TT,RR] = meshgrid(DELTA_THETA, DELTA_RADIUS);
[X,Y,Z] = pol2cart(TT,RR, mean_delta_r);
figure()
surf(X,Y,Z)
ax = gca;
ax.ZAxis.Scale ="log";

colorbar
set(gca,'ColorScale','log')
%%
% figure()
% for i_r=1:4
%     subplot(2,2,i_r)
% %     plot(delta_dtaus0, delta_r(:, i_r))
%     loglog(DELTA_RADIUS, delta_r(:, i_r), '.-')
%     title(sprintf('Receiver %d', i_r))
%     xlabel(strcat("\Delta ", PARAM0))
%     ylabel("Position Error (m)")
%     xlim([DELTA_RADIUS(1) DELTA_RADIUS(end)])
%     grid on
% end


%% Figures distance receivers-receivers error
% drr = pdist2(receivers, receivers);
% 
% drr_expanded = permute( ...
%     repmat(drr, [1, 1, length(rt_sols)]), ...
%     [3 1 2]); 
% 
% drr0 = zeros(length(DELTA_RADIUS), length(DELTA_THETA), length(receivers), length(receivers));
% 
% for i=1:length(DELTA_RADIUS)
%     drr0(i,:,:) = pdist2(squeeze(rt_sols(i,:,:)), squeeze(rt_sols(i,:,:)));
% end
% 
% delta_drr = sqrt((drr_expanded - drr0).^2);
% index_ij = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
% subtitles = ['1-2'; '1-3'; '1-4'; '2-3'; '2-4'; '3-4'];
% figure()
% for i_rr=1:length(index_ij)
%     i=index_ij(i_rr,1);
%     j=index_ij(i_rr,2);
%     subplot(2,3,i_rr)
%     loglog(DELTA_RADIUS, delta_drr(:, i, j), '.-')
%     title(subtitles(i_rr,:))
%     xlabel(strcat("\Delta ", PARAM0))
%     ylabel("Distance Error (m)")
%     xlim([DELTA_RADIUS(1) DELTA_RADIUS(end)])
%     grid on
% end

%% Figure 3D solutions
%figures.pyramid_3d_solution(receivers, receivers_0, squeeze(xsols(end,:,:)))
