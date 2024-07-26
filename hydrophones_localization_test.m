%% Step 4: 2 receivers N-pings locations
clear all
close all
clc

% 639
% xs = [1258.825 1385.942 1308.191 1214.396 1439.01 1383.474];
% ys = [748.713 739.447 867.322 822.844 854.348 919.213]; 

%637
xs = [477.613 414.672 340.623 550.427 494.891 381.35 ];
ys = [164.939 283.547 202.005 255.749 346.558 138.994];

zs = [5 5 5 5 5 5];
ns = length(xs);

%% Known parameters
% set source location
s = [xs(:).';ys(:).';zs(:).'];

% sound-speed [m/s]
c0 = 1500;

%% Test Data

% pyramid position
pyramid_x = 444.291;
pyramid_y = 233.51;
pyramid_z =  30.21;

% receiver locations
xr = [-1 -0 1 0] + pyramid_x;
yr = [-1 .5 -1 -0.5] + pyramid_y;
zr = [0 0 0 -1] + pyramid_z;

% yr = [-1 .5 -1 -0.5] + pyramid_y + rand(1,4)*1-.5;
% zr = [0 0 0 -1] + pyramid_z + rand(1,4)*1-.5;
% xr = [-1 -0 1 0] + pyramid_x + rand(1,4)*1-.5;

nr = length(xr);

r = [xr; yr; zr];

% synthetic travel time and delays
tau_mes = sqrt( (xr(:)-xs(:).').^2  + (yr(:)-ys(:).').^2 + (zr(:)-zs(:).').^2 ) / c0;

dtau_mes = NaN(nr,nr,ns);

for i_e = 1:ns
    dtau_mes(:,:,i_e) = tau_mes(:,i_e)-tau_mes(:,i_e).';
end



%% Variables Initialization
% receiver locations

r1_xyz = [-1  -1    0];
r2_xyz = [ 0  .5    0];
r3_xyz = [1   -1   0];
r4_xyz = [0  -.5  -1];

r0 = ([r1_xyz; r2_xyz; r3_xyz; r4_xyz] + [pyramid_x pyramid_y pyramid_z])';

%% Initial guess is real position with noise

%r0 = r + (rand(3,4)*2 -1) ;


%% swap initial receivers position
r0 = r0(:,[1 3 4 2]);


%% None Linear constrainte (distance between hydrophones)
max_distance = 3;
nonlcon = @(X) nonlinear_constraint.max_distance(X, max_distance);


%% optimize

% from percent of the period of the max freq
fmax = 256e3/2;
Tmin = 1/fmax;
tolfun = Tmin/100;

%cfun = @(X) loccostfun2(s,X,c0,dtau_mes);

% cfun = @(X) lowcost_R(X, s, c0, dtau_mes);
cfun = @(X) lowcost_functions.R(X, s, c0, dtau_mes);


%xsol = fmincon(cfun,r0,[],[],DD,dd);
%xsol = fmincon(cfun,r0,[],[],[],[],[], [], nonlcon);
xsol = fmincon(cfun,r0);

%% PLot
receivers_label = [1 2 3 4];
figure
hold on

% Ping Location
% plot3(xs,ys,zs,'kp','MarkerFaceColor','k')

% Initial Guess
P = r0;
plot3(P(1,:),P(2,:),P(3,:),'k.','MarkerSize',10)
trisurf(boundary(P',1),P(1,:)',P(2,:)',P(3,:)','FaceColor','black','FaceAlpha',0.1)
for i=1:nr
   text(P(1,i),P(2,i),P(3,i),['   ' ...
    num2str(receivers_label(i))],'HorizontalAlignment','left','FontSize',8);
end

% real position
P = r;
plot3(P(1,:),P(2,:),P(3,:),'blue.','MarkerSize',10)
trisurf(boundary(P',1),P(1,:)',P(2,:)',P(3,:)','FaceColor','blue','FaceAlpha',0.1)
for i=1:nr
   text(P(1,i),P(2,i),P(3,i),['   ' ...
    num2str(receivers_label(i))],'HorizontalAlignment','left','FontSize',8);
end


% Solutions
P = xsol;
plot3(P(1,:),P(2,:),P(3,:),'r.','MarkerSize',10)
trisurf(boundary(P',1),P(1,:)',P(2,:)',P(3,:)','FaceColor','red','FaceAlpha',0.1)
for i=1:nr
   text(P(1,i),P(2,i),P(3,i),['   ' ...
    num2str(receivers_label(i))],'HorizontalAlignment','left','FontSize',8);
end
hold off
set(gca,'ZDir','reverse')
grid on
% legend('Sources', 'Guess', 'Solution')
view(3)
daspect([1 1 1])
shg

%% test
dr1 = sqrt(sum((xsol(:,1) - xsol(:,[1 2 3 4])).^2));
dr2 = sqrt(sum((xsol(:,2) - xsol(:,[1 2 3 4])).^2));
dr3 = sqrt(sum((xsol(:,3) - xsol(:,[1 2 3 4])).^2));
dr4 = sqrt(sum((xsol(:,4) - xsol(:,[1 2 3 4])).^2));

dr = [dr1; dr2; dr3; dr4];

disp('Distance between hydrophone')
disp('      1     2     3     4')
fprintf('%d: %.2f  %.2f  %.2f  %.2f\n', [(1:4)', dr]')

% distance between hydrophone r < 2.5, c0 = 1465;
%       1     2     3     4
% 1: 0.00  1.93  2.03  1.53
% 2: 1.93  0.00  2.03  1.81
% 3: 2.03  2.03  0.00  1.76
% 4: 1.53  1.81  1.76  0.00

% distance between hydrophone Aucune contrainte, c0 = 1465;
%       1     2     3     4
% 1: 0.00  1.93  2.03  1.54
% 2: 1.93  0.00  2.03  1.82
% 3: 2.03  2.03  0.00  1.79
% 4: 1.54  1.82  1.79  0.00

% Distance between hydrophone c0=1500
%       1     2     3     4
% 1: 0.00  1.98  2.08  1.56
% 2: 1.98  0.00  2.08  1.85
% 3: 2.08  2.08  0.00  1.80
% 4: 1.56  1.85  1.80  0.00

% xsol(:, 1) c0=1500
%   436.4800
%   213.5600
%    34.3800

% % Distance between hydrophone c0=1450
%       1     2     3     4
% 1: 0.00  1.91  2.01  1.51
% 2: 1.91  0.00  2.01  1.79
% 3: 2.01  2.01  0.00  1.75
% 4: 1.51  1.79  1.75  0.00

% xsol(:, 1) c0=1450
%   436.5200
%   213.5500
%    34.3200

