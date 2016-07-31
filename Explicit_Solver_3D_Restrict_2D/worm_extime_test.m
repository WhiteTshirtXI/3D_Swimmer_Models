clear
%  3D Stokes code with an 
%  with immersed boundary worm
% add path for the source code
%
addpath('./src/');
% info for a restart
%
restart = false;
lastsaved = .5;
% define the domain and the grid spacing
%
Lx = 2;
Ly = 2;
Lz = 2;
% center the domain at the origin
%
xmin=-Lx/2;
ymin=-Ly/2;
zmin=-Lz/2;
% number of grid points and the grid spacing
%
Nx = 64;
Ny = 64;
Nz = 64;
dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;
% Viscoelastic Fluid Parameters
%
lam = 0;        % relaxation time
diffconst = 0;   % diffusion =  (diffconst*dx)^2
xi=0.5;           % polymer viscosity/solvent viscosity
% worm paramters
%
L  = 1.0;                      % length of the worm [mm]
Tper = 1.0;                    % period of the gait [s]
Nprefact=1.0;                  % ratio of ds to dx
N  = round( L/(Nprefact*dx) ); % number of points to discretize the worm
ds = L/(N-1);                  % length of each segment
s  = ds*(0:(N-1))';            % arclength coordinate of the worm
kb = 2;                        % bending stiffness
ks = 2500;                     % stretching stiffness
% curvature function that defines the gait
%
k0 = 2.0;
kappa_fun = @(s,t)(k0*sin(2*pi/Tper*t + pi*s));
% time stepping information
%
dt     = 2.5e-5;           % time step [s]
Tend   = 1;              % end time [s]
Nt     = round(Tend/dt);  % number of time steps to take
saveit = round(0.01/dt);  % frequency of output swimmer positions
saveall = 10*saveit;      % frequency of output all data
% output locations
%
datadir    = './data';  
runname    = 'exworm_3D';
fileprefix = sprintf('%s_n%03d',runname,Ny);
paramfile  = sprintf('%s/PARAMS_%s.txt',datadir,fileprefix);
%%%%%%%%%%-----END INPUT PARAMETERS-----%%%%%%%%%%
% pack up the grid parameters into a data structure
%
grid.Lx   = Lx;
grid.Ly   = Ly;
grid.Lz   = Lz;
grid.xmin = xmin;
grid.ymin = ymin;
grid.zmin = zmin;
grid.Nx   = Nx;
grid.Ny   = Ny;
grid.Nz   = Nz;
grid.N    = N;
grid.ds   = ds;
grid.dx   = dx;
grid.dy   = dy;
grid.dz   = dz;
% preallocate for speed
%
outcount_tot = round(Tend/(saveit*dt))+1;
XTworm=zeros(N,3,outcount_tot);
% initialize swimmer body position
%
kappa0 = kappa_fun(s,0);
X = initialize_worm(kappa0, ds);
% write parameters to file
%
fileID = fopen(paramfile,'w');
fprintf(fileID,'Lx = %d\n',Lx);
fprintf(fileID,'Ly = %d\n',Ly);
fprintf(fileID,'Lz = %d\n',Lz);
fprintf(fileID,'Nx = %d\n',Nx); 
fprintf(fileID,'Ny = %d\n',Ny);
fprintf(fileID,'Nz = %d\n',Nz);
fprintf(fileID,'dt = %0.8f\n',dt);
fprintf(fileID,'N (worm points) = %f\n',N);
fprintf(fileID, ' bending stiffness, kb = %4.4f\n',kb);
fprintf(fileID, ' stretching stiffness, ks =  %4.4f\n',ks);
fprintf(fileID, ' End time = %f\n',Tend);
fclose(fileID);
% initialize output counter and time -- overwritten if restart
%
outcount=1;
t = 0.0;
% if this is a restart, reload from the appropriate data file
%
if(restart)
    t=lastsaved;
    outcount=round(t*round(1/(saveit*dt))+1)
    fout1 = sprintf('%s/%s_t%f.mat',datadir,fileprefix,t);
    load(fout1);
    Nt   = round((Tend-t)/dt);
    X = XTworm(:,:,outcount);
end
% begin main loop in time
%
tic
for tint=0:Nt
    fprintf('Time step %i of %i, time=%f \n',tint,Nt,t);
    % record the worm positons
    %
    if(mod(tint,saveit)==0)
        XTworm(:,:,outcount)=X;
        outcount=outcount+1;
	fprintf('  recording output number %g \n',outcount-1);
    end
    % save all the data 
    % modify save need shat
    if(mod(tint,saveall) == 0 && t~=0)
        foutw = sprintf('%s/%s_t%f.mat',datadir,fileprefix,t);
        save(foutw,'U','Uw','XTworm');
    end
    % compute the curvature at the current time
    kappa0 = kappa_fun(s,t);
    % set the external body forces to zero for stokes solve
    fbhat = zeros(Nx, Ny, Nz, 3);
    % advance swimmer in time using forward euler
    [X,Uw,U,Eb,uhat] = EXstep_stokes(X,dt,fbhat,ks,kb,kappa0,grid);
    % fprintf('  Elastic energy = %g \n \n',Eb);
    % update time
    t = t+dt; 
    %toc
end
comp_time = toc;
fprintf('total computation time  = %g \n',comp_time);
