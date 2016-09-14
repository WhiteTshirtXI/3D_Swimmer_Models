clear
%  2D Stokes code with an 
%  with immersed boundary worm using implicit-time stepping

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
Nx = 16;
Ny = 16;
Nz = 16;
dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;

%fluid parameters
lam = 0.3;
xi = 0.5;
diffconst = 8;
%tend = 10.0;
%t0 = 0;                %<--- i think this is unnecessary
%savetime = 1;

% worm paramters
%
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
%kappa_fun = @(s,t)(k0*sin(2*pi/Tper*(Tper-t) + pi*(L-s)));

% time stepping information
%
dt     = 1e-3;           % time step [s]
Tend   = 3;              % end time [s]

Nt     = round(Tend/dt);  % number of time steps to take
saveit = round(0.01/dt);  % frequency of output swimmer positions
saveall = 10*saveit;      % frequency of output all data

% solver tolerances
%
rtol     = 1e-3;    % relative tolerance for objective function
rXtol    = 1e-6;   % relative tolerance for function values 
gmrestol = 5e-3;    % gmres tolerance for jacobian solve
  

% output locations
%
datadir    = './data';  
runname    = 'imworm_3D_VE';
fileprefix = sprintf('%s_n%03d',runname,Ny);
paramfile  = sprintf('%s/PARAMS_%s.txt',datadir,fileprefix);

% preallocate for speed
%
outcount_tot=round(Tend/(saveit*dt))+1;
XTworm=zeros(N,3,outcount_tot);

% initialize swimmer body position
%
kappa0 = kappa_fun(s,0);
X = initialize_worm(kappa0,ds);
%X = initialize_worm2(X,kappa0,ds);
   
% Obtain 4-roll Parameters, Stress Tensor, and Right-Hand Side Equation
[grid,params,Shat,newRHS] = get_4roll_inputs_3d(Ny,lam,xi,diffconst);   
% pack up the grid parameters into a data structure
%
grid.N    = N;
grid.ds   = ds;
grid.dy   = dy;
grid.dz   = dz;
% write parameters to file
%
fileID = fopen(paramfile,'w');
fprintf(fileID,'Lx = %f\n',Lx);
fprintf(fileID,'Ly = %f\n',Ly);
fprintf(fileID,'Lz = %f\n',Lz);
fprintf(fileID,'Nx = %f\n',Nx);
fprintf(fileID,'Ny = %d\n',Ny);
fprintf(fileID,'Nz = %f\n',Nz);
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
    if(isnan(Shat(1,1,1,1)))  % if the stress blows up stop running
        sprintf('warning, nans');
        break
    end
    fprintf('Time step %i of %i, time=%f \n',tint,Nt,t);

    % record the worm positons
    %
    if(mod(tint,saveit)==0)
        XTworm(:,:,outcount)=X;
        outcount=outcount+1;
	fprintf('  recording output number %g \n',outcount-1);
    end
    
    % save all the data 
    %
    if( mod(tint,saveall)==0 && t~=0)
        foutw = sprintf('%s/%s_t%f.mat',datadir,fileprefix,t);
	    save(foutw,'U','Uw','XTworm','Shat');
    end
    
    % compute the curvature at the current time
    %
    kappa0 = kappa_fun(s,t);
    if(lam == 0)  % Stokes solve
        fbhat = zeros(Nx, Ny, Nz, 3);   
    else  % Compute VE force
        fbhat = get_veforcehat_3d(Shat,xi,grid);
    end
    % advance swimmer in time using backward euler
    [X,Uw,U,Uhat,output] = IMstep_stokes_newton(X,dt,fbhat,ks,kb,kappa0,grid,rtol,rXtol,gmrestol,lam,xi);
    % Update the stress tensor if not solving Stokes
    if(lam~=0)
        [Shat, newRHS] = update_Shat_3d(Uhat,grid,Shat,params.nu,dt,lam,newRHS);
    end    
    % update time
    %
    t = t+dt;
        
    %toc
end

comp_time = toc;

fprintf('total computation time  = %g \n',comp_time);
