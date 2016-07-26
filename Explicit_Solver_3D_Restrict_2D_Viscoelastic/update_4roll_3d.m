function [Uhat,Shat,newRHS]=update_4roll_3d(Ny,lam,xi,diffconst,t0,tend,savepersec)

[grid,params,Shat,newRHS]=get_4roll_inputs_3d(Ny,lam,xi,diffconst);

% unpack parameters

lam = params.lam;
xi = params.xi;
nu = params.nu;
dc = params.diffconst;
dt = params.dt;

Lx = grid.Lx;
Ly = grid.Ly;
Lz = grid.Lz;
xmin = grid.xmin;
ymin = grid.ymin;
zmin = grid.zmin;
Nx = grid.Nx;
Ny = grid.Ny;
Nz = grid.Nz;
dx = grid.dx;


datadir=sprintf('./SOB_4roll_3d/lam%1.1f',lam);

runname    = '4roll_3d';

fileprefix = sprintf('%s_n%03d_lam%1.2f_dc%d',runname,Ny,lam,dc);

[status,message,messageid] = mkdir(datadir);  % make directory if it is not already present

paramfile  = sprintf('%s/PARAMS_%s.txt',datadir,fileprefix)  % parameter file


% write parameters to file
%
fileID = fopen(paramfile,'w');
fprintf(fileID,'lam = %f\n',lam);
fprintf(fileID,'xi = %f\n',xi);
fprintf(fileID,'nu = %f\n',nu);
fprintf(fileID,'Lx = %f\n',Lx);
fprintf(fileID,'Ly = %f\n',Ly);
fprintf(fileID,'Lz = %f\n',Lz);
fprintf(fileID, ' xmin = %f\n',xmin);
fprintf(fileID, ' ymin = %f\n',ymin);
fprintf(fileID, ' zmin = %f\n',zmin);
fprintf(fileID,'dx = %f\n',dx);
fprintf(fileID, 'Ny = %d\n',Ny);
fprintf(fileID,'dt = %0.8f\n',dt);
fprintf(fileID,'save per sec. = %1.4\n',savepersec);
fprintf(fileID, ' End time = %f\n',tend);

fclose(fileID);

Nt     = round((tend-t0)/dt);  % number of time steps to take
persecs  = round(1/dt);  % frequency of output swimmer positions
saveall = persecs/savepersec;


if t0==0
    t=0;
else   % load in old data for restart
    fin = sprintf('%s/%s_t%1.2f.mat',datadir_in,fileprefix_in,t0);
    load(fin);
    t=t0;
end

frollFhat_3d = forcing_4roll_3d(grid);  % compute the constant (in time) forcing

% Start the time stepping

for tint=0:Nt    
    if(isnan(Shat(1,1,1,1)))  % if the stress blows up stop running
        sprintf('warning, nans');
        break
    end
     
    % save all the data 
    %
    if( mod(tint,saveall)==0 && t~=0)
        foutw = sprintf('%s/%s_t%f.mat',datadir,fileprefix,t);
	    save(foutw,'Uhat','Shat','newRHS');
        fprintf('Time step %i of %i, time=%f \n',tint,Nt,t);
  
        
    end
    
    
    if(lam==0)  % Stokes solve
        
        fbhat = frollFhat_3d;   
        
    else  % Compute VE force
        
        veforce_hat = get_veforcehat_3d(Shat,xi,grid);
        
        fbhat = veforce_hat+frollFhat_3d;
        
    end
    
    % Find Velocity
    % Solve Stokes equation with VE stress
    
    Uhat = stokes_solve_fourier_3d(fbhat,grid.Lx,grid.Ly,grid.Lz);
    
    % Update the stress tensor if not solving Stokes
    if(lam~=0)
        
        [Shat, newRHS] = update_Shat_3d(Uhat,grid,Shat,nu,dt,lam,newRHS);
       
    end
    
     % update time
    %
    t = t+dt;
        
    
end



