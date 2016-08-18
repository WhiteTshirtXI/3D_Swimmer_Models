function [params,Shat,newRHS]=get_4roll_inputs(Ny,lam,xi,diffconst)

% This function lamll output the needed parameters for the update
% command form is [grid,params,Shat,newRHS]=get_4roll_inputs(Ny,lam,xi,nu)


% define the domain and the grid spacing
%
Lx = 2;
Ly = 2;

xmin=-.5;
ymin=-.5;

K  = Lx/Ly;  % for non-square domains

Nx = K*Ny;
dx = Ly/Ny; % dx==dy no matter if it is rectangular or square

% diffusion based on gs and grid spacing

nu = (diffconst*dx).^2;


% time stepping information
%
dt =(.01/(2^(log2(Ny)-6)));  % This gives dt = 0.01 for Ny =  64 and divides that by 2 as you refine (or multiplies as you coarsen



%%% --- problem parameters --- %%%
params.lam = lam;   % relaxation time
params.xi = xi;     % polymer to solvent viscosity ratio
params.nu = nu;     % diffusion
params.dt = dt;     % time step
params.diffconst = diffconst; % 

%%% --- grid parameters --- %%%%
%
grid.Lx   = Lx;
grid.Ly   = Ly;
grid.xmin = xmin;
grid.ymin = ymin;
grid.Nx   = Nx;
grid.Ny   = Ny;
grid.dx   = dx;

%%%---- initial data for stress ---%%%%

S(:,:,1)=zeros(Nx,Ny);
S(:,:,2)=zeros(Nx,Ny);
S(:,:,3)=zeros(Nx,Ny);
Shat = fft2(S);  % in fourier space

newRHS = zeros(Nx,Ny,3);  % initial Right hand side for time-step
