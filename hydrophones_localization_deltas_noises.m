clear all
close all
clc

%% Deltas Noise

% PARAM0 = 'delays'; % delays, geometry, pyramid
% DELTAS_PARAM0 = 1e-6*(10.^(0:.1:4)); % noise amplitude vector
% DELTAS_PARAM0 = 0*(10.^(0:1)); % noise amplitude vector

% PARAM0 = 'geometry'; % delays, geometry, pyramid
% DELTAS_PARAM0 = 1e-3*(10.^(0:.25:5)); % noise amplitude vector; % noise amplitude vector

% PARAM0 = 'sources'; % delays, geometry, pyramid
% DELTAS_PARAM0 = 1e-3*(10.^(0:.5:5)); % noise amplitude vector; % noise amplitude vector

PARAM0 = 'pyramid'; % delays, geometry, pyramid
DELTAS_PARAM0 = 1e-2*(10.^(0:.25:4)); % noise amplitude vector; % noise amplitude vector

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


%% Solving for a range of vertical r0 position.
    
xsols = zeros( ...
    horzcat( ...
    length(DELTAS_PARAM0), size(receivers) ...
    ) ...
    );

for i=1:length(DELTAS_PARAM0)
    
    fprintf('iteration %d\n', i)

    dtaus0 = delays_measured;
    r0 = receivers_0;
    sources_0 = sources;
    
    switch PARAM0 
        % Incertitude on measurements
        case 'delays'
            % noise on delays
            dtaus0 = delays_measured + tools.rand0(size(delays_measured), DELTAS_PARAM0(i));
        
        case 'sources'
             % noise on sources position
             sources_0 = sources_0 + tools.rand0(size(sources_0), DELTAS_PARAM0(i)).*[1 1 0.005];
        
        % Incertitude on initial Guess.
        case 'geometry'
             % noise on receiver geometry
             r0 = receivers_0 + tools.rand0(size(receivers_0), DELTAS_PARAM0(i));

        case 'pyramid'
            % noise on pyramid position
            r0 = receivers_0 + DELTAS_PARAM0(i) .*[1 1 0.1];
    end
   
    xsol = solver.R(r0, sources_0, sound_speed_0, dtaus0);
   
    xsols(i,:,:) = xsol;

end


%% Figures Receivers position error

%permute(permute(xsols, [2 3 1]) - receivers, [3 1 2]);

delta_r = sqrt(sum(( ...
    permute(permute(xsols, [2 3 1]) - receivers, [3 1 2]) ...
    ).^2, 3));

figure()
for i_r=1:4
    subplot(2,2,i_r)
%     plot(delta_dtaus0, delta_r(:, i_r))
    loglog(DELTAS_PARAM0, delta_r(:, i_r), '.-')
    title(sprintf('Receiver %d', i_r))
    xlabel(strcat("\Delta ", PARAM0))
    ylabel("Position Error (m)")
    xlim([DELTAS_PARAM0(1) DELTAS_PARAM0(end)])
    grid on
end


%% Figures distance receivers-receivers error
drr = pdist2(receivers, receivers);

drr_expanded = permute( ...
    repmat(drr, [1, 1, length(xsols)]), ...
    [3 1 2]); 

drr0 = zeros(length(DELTAS_PARAM0), length(receivers), length(receivers));

for i=1:length(DELTAS_PARAM0)
    drr0(i,:,:) = pdist2(squeeze(xsols(i,:,:)), squeeze(xsols(i,:,:)));
end

delta_drr = sqrt((drr_expanded - drr0).^2);
index_ij = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
subtitles = ['1-2'; '1-3'; '1-4'; '2-3'; '2-4'; '3-4'];
figure()
for i_rr=1:length(index_ij)
    i=index_ij(i_rr,1);
    j=index_ij(i_rr,2);
    subplot(2,3,i_rr)
    loglog(DELTAS_PARAM0, delta_drr(:, i, j), '.-')
    title(subtitles(i_rr,:))
    xlabel(strcat("\Delta ", PARAM0))
    ylabel("Distance Error (m)")
    xlim([DELTAS_PARAM0(1) DELTAS_PARAM0(end)])
    grid on
end

%% Figure 3D solutions
%figures.pyramid_3d_solution(receivers, receivers_0, squeeze(xsols(end,:,:)))
