function rockyriver_trajectories()

clear all
close all

set(0,'DefaultLineLineWidth',1.5, ...
      'DefaultAxesFontSize',14, ...
      'DefaultTextFontSize',14, ...
      'DefaultAxesLineWidth',1.5);

addpath('../MatlabFire')

% Parameters

tend = datenum('2007-12-09_07:02:00','yyyy-mm-dd_HH:MM:SS')
tstart = datenum('2007-12-09_07:12:00','yyyy-mm-dd_HH:MM:SS')
cloudlocation = [136.86 -35.83]; % long lat
cloudheight = 1000;            % m
cloudthickness = 0;          % m
nparticles = 1;

odeset('AbsTol',1e-2,'RelTol',1e-2);

% Data files

pathname = '/home/tmattner/Research/FireWX/Data/Mika/KIruns/RR19marFRST/';
filename = 'wrfout_d04_2007-12-09_06:00:00';
datafile = strcat(pathname,filename)

% Open the datafile.

data = wrffire(datafile)

% Local time (days)

t = datenum(data.getTimes,'yyyy-mm-dd_HH:MM:SS');

% Find time indices

[~,kstart] = min(abs(t - tstart));
[~,kend] = min(abs(t - tend));
kstep = sign(kend - kstart);
kfile = kstart:kstep:kend;
nk = length(kfile)

% Cell-centred coordinates in m. 

Lx = (data.nx - 1)*data.dx;
Ly = (data.ny - 1)*data.dy;
[x3,y3,~] = ndgrid(0.5*data.dx + linspace(0,Lx,data.nx), ...
                   0.5*data.dy + linspace(0,Ly,data.ny), ...
                   zeros(1,data.nz));
z3 = data.interpn(data.getVar('altitude',kstart));   

% Find indices corresponding to centre of particle cloud

x2 = data.getVar('XLONG',kstart);
y2 = data.getVar('XLAT',kstart);
[i,j] = data.coord2sub(x2,y2,cloudlocation(1),cloudlocation(2));
z2 = z3(i,j,:);
[~,k] = min(abs(z2(:) - cloudheight));

% Initial conditions. Seed particles randomly in a cube

x0 = x3(i,j,k) + 2*cloudthickness*(rand(1,nparticles) - 0.5);
y0 = y3(i,j,k) + 2*cloudthickness*(rand(1,nparticles) - 0.5);
z0 = z3(i,j,k) + 2*cloudthickness*(rand(1,nparticles) - 0.5);

% Preallocate storage

x = zeros(nk,nparticles);
y = zeros(nk,nparticles);
z = zeros(nk,nparticles);

% Store initial conditions

x(1,:) = x0;
y(1,:) = y0;
z(1,:) = z0;

% Convert time to sec

t = 86400*t;

% Interpolate velocity and height to centered grid at
% k-th time.

uk = data.interpx(data.getVar('U',kfile(1)));
vk = data.interpy(data.getVar('V',kfile(1)));
wk = data.interpn(data.getVar('W',kfile(1)));
z3k = data.interpn(data.getVar('altitude',kfile(1)));
tk = t(kfile(1));

% Prescribe velocity field for testing

theta = atan2(y3-y0(1)-100,x3-x0(1)-100);
r = sqrt((x3-x0(1)-100).^2 + (y3-y0(1)-100).^2)/100*pi/3;

uk = -r.*sin(theta);
vk = r.*cos(theta);
wk(:,:,:) = 1;

% Loop through all times

for k = 1:nk-1
    tic;

    % Interpolate velocity and height to centered grid at
    % (k+1)-th time.

    ukp1 = data.interpx(data.getVar('U',kfile(k+1)));
    vkp1 = data.interpy(data.getVar('V',kfile(k+1)));
    wkp1 = data.interpn(data.getVar('W',kfile(k+1)));
    z3kp1 = data.interpn(data.getVar('altitude',kfile(k+1)));
    tkp1 = t(kfile(k+1));

    % Prescribe velocity field for testing

    ukp1 = -r.*sin(theta);
    vkp1 = r.*cos(theta);
    wkp1(:,:,:) = 1;

    % Integrate particle positions from k-th to (k+1)-th time

    [~,tmp] = ode23(@velocity,[tk tkp1],[x(k,:)'; y(k,:)'; z(k,:)']);

    % Store positions at (k+1)-th time

    x(k+1,:) = tmp(end,1:nparticles);
    y(k+1,:) = tmp(end,nparticles+1:2*nparticles);
    z(k+1,:) = tmp(end,2*nparticles+1:3*nparticles);

    % Move data at (k+1)-th time to k-th time in preparations for
    % next step

    uk = ukp1;
    vk = vkp1;
    wk = wkp1;
    z3k = z3kp1;
    tk = tkp1;

    toc
end

% Modify to plot or output as required

plot3(x,y,z,'bo-')
hold on
plot3(x0,y0,z0,'ro')
hold off
xlabel('Easting (m)')
ylabel('Northing (m)')
print -djpeg tmp.jpg
saveas(gcf,'tmp.fig','fig')

% Close the file

data.delete;

  function dqdt = velocity(t,q)

    % VELOCITY interpolates gridded velocity data to particle position
    % fpor any time between the k-th and (k+1)-th.

    % Extract position data from column vector q

    xq = q(1:nparticles);
    yq = q(nparticles+1:2*nparticles);
    zq = q(2*nparticles+1:3*nparticles);

    % Interpolate velocity at particle positions xq, yq, zq

    dxdtk   = data.interp3d(x3,y3,z3k,  uk,  xq,yq,zq,'spline');
    dydtk   = data.interp3d(x3,y3,z3k,  vk,  xq,yq,zq,'spline');
    dzdtk   = data.interp3d(x3,y3,z3k,  wk,  xq,yq,zq,'spline');
    dxdtkp1 = data.interp3d(x3,y3,z3kp1,ukp1,xq,yq,zq,'spline');
    dydtkp1 = data.interp3d(x3,y3,z3kp1,vkp1,xq,yq,zq,'spline');
    dzdtkp1 = data.interp3d(x3,y3,z3kp1,wkp1,xq,yq,zq,'spline');

    % Interpolate particle velocity to time t

    dxdt = (t - tkp1)/(tk - tkp1)*dxdtk + (t - tk)/(tkp1 - tk)*dxdtkp1;
    dydt = (t - tkp1)/(tk - tkp1)*dydtk + (t - tk)/(tkp1 - tk)*dydtkp1;
    dzdt = (t - tkp1)/(tk - tkp1)*dzdtk + (t - tk)/(tkp1 - tk)*dzdtkp1;

    % Store partile velocities in output column vector dqdt

    dqdt = [dxdt; dydt; dzdt];

  end

end


