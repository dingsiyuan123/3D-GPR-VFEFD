%           Vector Finite Frequency Domain Simulation for GPR              %
%                   Analytical Solution of dipole                          %
% ------------------------ 2024.10 Siyuan Ding --------------------------- %
% Frequency-domain vector finite element forward modeling of 3D GPR data
% using exact PML absorbing boundary conditions %
clc;clear;tic
%% Define Model
Nx = 80; 
Ny = 80; 
Nz = 40;
ds = 0.05;
dx = ds; 
dy = ds;
dz = ds;
miu = 4*pi*1e-7;
ep0 = 8.85*1.0e-12;
ep=4;                % relative permittivity
sig=0.003;           % conductivity
%% Antenna
freq =1e8;
%% Analytical Solution
eta=120*pi/2;
Io=1;                % current in dipole
l=1;                 % length of dipole
% c = 1e8;
% lambda=c/freq;       % wavelength
% k=(2*pi)/lambda;     % propagation constant
omega = 2 * pi *freq;
k = sqrt(omega^2.*ep0*miu*(ep-1i*sig/omega/ep0));
x = Nx * dx/2;
y = Ny * dy/2;
z = Nz * dz/2;
xx = [-x:dx:-(0+dx) dx/50 (0+dx):dx:x];    % avoid zero
yy = [-y:dy:-(0+dy) dy/50 (0+dy):dy:y];
zz = [-z:dz:-(0+dz) dz/50 (0+dz):dz:z];
[X,Y,Z] = meshgrid(xx,yy,zz);

r = sqrt( X.^2 + Y.^2 + Z.^2 );            % distance from centre in grid
theta = atan( sqrt(X.^2 + Y.^2)./Z);       % theta at each grid points
phi = atan2(Y,X);

Er = ((eta*Io*l)./(2*pi*r.^2)).*(1+1./(1j*k*r)).*cos(theta).*exp(-1j*k*r);
Etheta = ((1j*eta*k*Io*1)./(4*pi*r)).*(1+1./(1j*k*r)-1./(k*r).^2).*sin(theta).*exp(-1j*k*r);
Ephi = 0;

Ex = Er.*sin(theta).*cos(phi) + Etheta.*cos(theta).*cos(phi) - Ephi*sin(phi);
Ey = Er.*sin(theta).*sin(phi) + Etheta.*cos(theta).*sin(phi) + Ephi*cos(phi);
Ez = Er.*cos(theta) - Etheta.*sin(theta);
%% Plot
size_t = ceil(size(Ex)/2);
fig1=figure(1);
axes1 = axes('Parent',fig1);
image1=slice(real(Ex),size_t(1),size_t(2),size_t(3));
set(image1,'EdgeColor','none','FaceColor','interp');
colormap jet
colorbar
caxis([-100, 100]);
xlabel('x');ylabel('y');zlabel('z');
set(axes1,'CLim',[-50 50],'ZDir','reverse');

fig2=figure(2);
axes2 = axes('Parent',fig2);
image2=slice(real(Ey),size_t(1),size_t(2),size_t(3));
set(image2,'EdgeColor','none','FaceColor','interp');
colormap jet
colorbar
caxis([-100, 100]);
xlabel('x');ylabel('y');zlabel('z');
set(axes1,'CLim',[-50 50],'ZDir','reverse');

fig3=figure(3);
axes3 = axes('Parent',fig3);
image3=slice(real(Ez),size_t(1),size_t(2),size_t(3));
set(image3,'EdgeColor','none','FaceColor','interp');
colormap jet
colorbar
caxis([-100, 100]);
xlabel('x');ylabel('y');zlabel('z');
set(axes1,'CLim',[-50 50],'ZDir','reverse');