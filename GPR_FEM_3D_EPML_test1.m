%           Vector Finite Frequency Domain Simulation for GPR              %
%                   Using EPML and Hexahedral mesh                                  %
% Governing Equations for 3D waves:
%             curl(1/t*curl(E))-k^2*t*E = -1i*omega*miu*J;
%             k^2 = omega^2*miu*(eps-1i*sigma/omega);
% ------------------------ 2024.10 Siyuan Ding --------------------------- %
% Frequency-domain vector finite element forward modeling of 3D GPR data
% using exact PML absorbing boundary conditions %
clc;clear;tic
addpath Solvers   % Solve the system of complex sparse linear equations using MUMPS. Please download from the website.
%% Define Model
Nx = 80; 
Ny = 80; 
Nz = 40;
npml = 1;           % number of layers in PML
Nx = Nx+2*npml;
Ny = Ny+2*npml;
Nz = Nz+2*npml;
ds = 0.05;          % discrete mesh spacing
dx = ds; 
dy = ds;
dz = ds;
epr(1:Ny,1:Nx,1:Nz)=4;                 % relative permittivity
sig(1:Ny,1:Nx,1:Nz)=0.003;             % conductivity
miu = 4*pi*1e-7;
ep0 = 8.85*1.0e-12;
%% Antenna and Location
freq =0.5e8;   
omega = 2*pi*freq;
nsitex = Nx/2;
nsitey = Ny/2;
nsitez = Nz/2;
%% Define EPML
wavenumber = sqrt(omega^2.*ep0*miu.*(epr-1i.*sig/omega/ep0));
% x direction
sigx = zeros(Ny,Nx,Nz);
for i=1:npml
    sigx(:,i,:)=1/(i-1/2)/dx;
end
sigx(:,Nx-npml+1:Nx,:)=sigx(:,npml:-1:1,:);

% y direction
sigy = zeros(Ny,Nx,Nz);
for i=1:npml
    sigy(i,:,:)=1/(i-1/2)/dy;
end
sigy(Ny-npml+1:Ny,:,:)=sigy(npml:-1:1,:,:);

% z direction
sigz = zeros(Ny,Nx,Nz);
for i=1:npml
    sigz(:,:,i)=1/(i-1/2)/dz;
end
sigz(:,:,Nz-npml+1:Nz)=sigz(:,:,npml:-1:1);

Sx =  1-(sigx*1i./wavenumber);
Sy =  1-(sigy*1i./wavenumber);
Sz =  1-(sigz*1i./wavenumber);
%%  Mesh 
Ne = Nx*Ny*Nz; % number of elements
NP = Nx*(Ny+1)*(Nz+1)+Ny*(Nx+1)*(Nz+1)+Nz*(Nx+1)*(Ny+1);% number of edges
%% Sort Elements and Edges
alpha_x = zeros(Ne,1);
alpha_y = zeros(Ne,1);
alpha_z = zeros(Ne,1);
beta_x = zeros(Ne,1);
beta_y = zeros(Ne,1);
beta_z = zeros(Ne,1);
k2 = zeros(Ne,1);

I12=zeros(Ne,12);
for l = 1:Nz
    for m = 1:Ny
        for n = 1:Nx
            % Sort the Elements
            NumNE=(l-1)*Ny*Nx+(m-1)*Nx+n;
            % Position of Source
            if l==nsitez && m==nsitey && n==nsitex
                site = Nx*(Nz+1)*(Ny+1)+Ny*(Nx+1)*(Nz+1)+(l-1)*(Nx+1)*(Ny+1)+(m-1)*(Nx+1)+n;
            end
            %%%%%%%%%%%%%%%%% x direction %%%%%%%%%%%%%%%%
            I12(NumNE,1) = (n-1)*(Ny+1)*(Nz+1)+(l-1)*(Ny+1)+m;
            I12(NumNE,2) = I12(NumNE,1)+1;
            I12(NumNE,3) = I12(NumNE,1)+Ny+1;
            I12(NumNE,4) = I12(NumNE,1)+Ny+2;
            %%%%%%%%%%%%%%%%% y direction %%%%%%%%%%%%%%%%
            I12(NumNE,5) = Nx*(Nz+1)*(Ny+1)+(m-1)*(Nx+1)*(Nz+1)+(n-1)*(Nz+1)+l;
            I12(NumNE,6) = I12(NumNE,5)+1;
            I12(NumNE,7) = I12(NumNE,5)+Nz+1;
            I12(NumNE,8) = I12(NumNE,5)+Nz+2;
            %%%%%%%%%%%%%%%%% z direction %%%%%%%%%%%%%%%%
            I12(NumNE,9)  = Nx*(Nz+1)*(Ny+1)+Ny*(Nx+1)*(Nz+1)+(l-1)*(Nx+1)*(Ny+1)+(m-1)*(Nx+1)+n;
            I12(NumNE,10) = I12(NumNE,9)+1;
            I12(NumNE,11) = I12(NumNE,9)+Nx+1;
            I12(NumNE,12) = I12(NumNE,9)+Nx+2;
            %%%%%%%%%%%% Sort the Parameters %%%%%%%%%%%%%
            k2(NumNE) = omega^2.*ep0*miu*(epr(m,n,l)-1i*sig(m,n,l)/omega/ep0);
            alpha_x(NumNE) = Sx(m,n,l)./(Sy(m,n,l).*Sz(m,n,l));
            alpha_y(NumNE) = Sy(m,n,l)./(Sx(m,n,l).*Sz(m,n,l));
            alpha_z(NumNE) = Sz(m,n,l)./(Sx(m,n,l).*Sy(m,n,l));
            beta_x(NumNE) = Sy(m,n,l).*Sz(m,n,l)./Sx(m,n,l);
            beta_y(NumNE) = Sx(m,n,l).*Sz(m,n,l)./Sy(m,n,l);
            beta_z(NumNE) = Sx(m,n,l).*Sy(m,n,l)./Sz(m,n,l);
        end
    end
end

%% Matrix Generation
K2=sparse(NP,NP);
Ke=[4 2 2 1;2 4 1 2;2 1 4 2;1 2 2 4]/36;
for j=1:4
    for k=1:4
        %%%%%%%%%%%% K11 %%%%%%%%%%%%%%%%%%%%%
        NK = I12(:,k); 
        NJ = I12(:,j);     
        ktemp = k2.*Ke(j,k)*dx.*dy.*dz.*beta_x;
        K2 = K2+sparse(NJ,NK,ktemp,NP,NP);
        %%%%%%%%%%%% K22 %%%%%%%%%%%%%%%%%%%%%
        NK = I12(:,k+4); 
        NJ = I12(:,j+4);   
        ktemp = k2.*Ke(j,k)*dx.*dy.*dz.*beta_y;
        K2 = K2+sparse(NJ,NK,ktemp,NP,NP);
        %%%%%%%%%%%% K33 %%%%%%%%%%%%%%%%%%%%%
        NK = I12(:,k+8);
        NJ = I12(:,j+8);   
        ktemp = k2.*Ke(j,k)*dx.*dy.*dz.*beta_z;
        K2 = K2+sparse(NJ,NK,ktemp,NP,NP); 
    end
end
% toc
% tic
K1=sparse(NP,NP);
k1=[2 -2 1 -1;-2 2 -1 1;1 -1 2 -2;-1 1 -2 2]/6;
k2=[2 1 -2 -1;1 2 -1 -2;-2 -1 2 1;-1 -2 1 2]/6;
k3=[2 1 -2 -1;-2 -1 2 1;1 2 -1 -2;-1 -2 1 2]/6;
for j=1:4
    for k=1:4        
        %%%%%%%%%    K11  %%%%%%
        NJ = I12(:,j);
        NK = I12(:,k);
        ktemp = alpha_z.*k1(j,k).*dx.*dz./dy +...
            alpha_y.*k2(j,k).*(dx.*dy./dz);
        K1 = K1+sparse(NJ,NK,ktemp,NP,NP);
        %%%%%%%%%    K22  %%%%%%
        NJ = I12(:,j+4);
        NK = I12(:,k+4);
        ktemp = alpha_x.*k1(j,k).*(dx.*dy./dz)+...
            alpha_z.*k2(j,k).*(dy.*dz./dx);
        K1 = K1+sparse(NJ,NK,ktemp,NP,NP);
        %%%%%%%%%    K33 %%%%%%
        NJ = I12(:,j+8);
        NK = I12(:,k+8);
        ktemp = alpha_y.*k1(j,k).*(dy.*dz./dx)+...
            alpha_x.* k2(j,k).*(dx.*dz./dy);
        K1 = K1+sparse(NJ,NK,ktemp,NP,NP);
        %%%%%%%%%%    K12 %%%%%%
        NJ = I12(:,j);  
        NK = I12(:,k+4);
        ktemp = -alpha_z.*k3(j,k).*dz;
        K1 = K1+sparse(NJ,NK,ktemp,NP,NP);
        %%%%%%%%%%    K21 %%%%%%
        NJ = I12(:,j+4);  
        NK = I12(:,k);
        ktemp = -alpha_z.*k3(k,j).*dz;
        K1 = K1+sparse(NJ,NK,ktemp,NP,NP);
        %%%%%%%%    K31 %%%%%%
        NJ = I12(:,j+8);  
        NK = I12(:,k);
        ktemp = -alpha_y.*k3(j,k).*dy;
        K1 = K1+sparse(NJ,NK,ktemp,NP,NP);
        %%%%%%%%    K13 %%%%%%
        NJ = I12(:,j);  
        NK = I12(:,k+8);
        ktemp = -alpha_y.*k3(k,j).*dy;
        K1 = K1+sparse(NJ,NK,ktemp,NP,NP);
        %%%%%%%    K23 %%%%%%
        NJ = I12(:,j+4);
        NK = I12(:,k+8);
        ktemp = -alpha_x.*k3(j,k).*dx;
        K1 = K1+sparse(NJ,NK,ktemp,NP,NP);
        %%%%%%%    K32 %%%%%%
        NJ = I12(:,j+8);  
        NK = I12(:,k+4);
        ktemp = -alpha_x.*k3(k,j).*dx;
        K1 = K1+sparse(NJ,NK,ktemp,NP,NP);
    end
end
A = K1-K2;
%% Source
b = sparse(NP,1);           
b(site) = -1i*omega*miu;           
time1 = toc
tic
%% Solve Equation
[x,info] = opti_zmumps(A,b); % Solve the system of complex sparse linear equations using MUMPS. Please download from the website.
% x = A\b;                   % Left division also works. While it will be difficult and take time if the size of the model is big.
time2 = toc
%% Post-processing
u = x;
%%%%%%%%%%%%%%   Ex  %%%%%%%%%%%%%%%
for k=1:Nz+1
    for i=1:Nx
        for j=1:Ny+1
            MEij=(i-1)*(Ny+1)*(Nz+1)+(k-1)*(Ny+1)+j;% edges in x direction
            Ex(j,i,k)=x(MEij);
        end
    end
end

%%%%%%%%% Ey %%%%%%%%%%%%%
for k=1:Nz+1
    for i=1:Nx+1
        for j=1:Ny
            MEij=Nx*(Nz+1)*(Ny+1)+(j-1)*(Nx+1)*(Nz+1)+(i-1)*(Nz+1)+k;% edges in y direction
            Ey(j,i,k)=x(MEij);
        end
    end
end
%%%%%%%%% Ez %%%%%%%%%%%%%
for k=1:Nz
    for i=1:Nx+1
        for j=1:Ny+1
            MEij=Nx*(Nz+1)*(Ny+1)+Ny*(Nx+1)*(Nz+1)+(k-1)*(Nx+1)*(Ny+1)+(j-1)*(Nx+1)+i;% edges in z direction
            Ez(j,i,k)=x(MEij);
        end
    end
end
%% Plot
fig1=figure(1);
axes1 = axes('Parent',fig1);
image1=slice(real(Ex),Nx/2,Ny/2,Nz/2);
set(image1,'EdgeColor','none','FaceColor','interp');
colormap jet
colorbar
caxis([-100, 100]);
xlabel('x');ylabel('y');zlabel('z');
set(axes1,'CLim',[-50 50],'ZDir','reverse');

fig2=figure(2);
axes2 = axes('Parent',fig2);
image2=slice(real(Ey),Nx/2,Ny/2,Nz/2);
set(image2,'EdgeColor','none','FaceColor','interp');
colormap jet
colorbar
caxis([-100, 100]);
xlabel('x');ylabel('y');zlabel('z');
set(axes1,'CLim',[-50 50],'ZDir','reverse');

fig3=figure(3);
axes3 = axes('Parent',fig3);
image3=slice(real(Ez),Nx/2,Ny/2,Nz/2);
set(image3,'EdgeColor','none','FaceColor','interp');
colormap jet
colorbar
caxis([-100, 100]);
xlabel('x');ylabel('y');zlabel('z');
set(axes1,'CLim',[-50 50],'ZDir','reverse');