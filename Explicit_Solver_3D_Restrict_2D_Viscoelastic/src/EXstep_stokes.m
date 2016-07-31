function [X,Uw,U,E,Uhat] = EXstep_stokes(Xn,dt,fbhat_ext,ks,kb,kappa,grid)
     
  % form the spreading operator
  spfactor = grid.ds/(grid.dx*grid.dx*grid.dx);
  Sm = spreadmatrix3_vc_vec(Xn,grid.dx,grid.Nx,grid.Ny,grid.Nz,grid.xmin,grid.ymin,grid.zmin);
    
  % Evaluate the forces at the current position
  [Fb,Kx] = bending_force_vec3(Xn,kappa,kb,grid.ds);  
  [Fs,St] = stretch_force_vec3(Xn,ks,grid.ds);
  F  = Fb + Fs;
  
  % compute the elastic energy of the worm
  E = kb*sum((Kx(2:end-1)-kappa(2:end-1)).^2) + ks*sum( (St(1:end-1)-1).^2);
  E = grid.ds*0.5*E;
   
  % spread forces
  fb = spfactor * reshape(Sm*F,grid.Nx,grid.Ny,grid.Nz,3);
  % Preallocate space for speed.
  fbhat = zeros(grid.Nx, grid.Ny, grid.Nz, 3);
  
  % compute total force in fourier space
  for d = 1:3
      fbhat(:,:,:,d) = fftn(fb(:,:,:,d));
  end
  fbhat = fbhat + fbhat_ext;
  
  % solve stokes in fourier space
  Uhat = stokes_solve_fourier_3d(fbhat,grid.Lx,grid.Ly,grid.Lz);
  % Preallocate space for speed.
  U = zeros(grid.Nz, grid.Ny, grid.Nz, 3);
  
  for i = 1:3
      U(:,:,:,i) = real(ifftn(Uhat(:,:,:,i)));
  end
  
  % interpolate velocity back to the IB points
  Uw = Sm'*reshape(U,grid.Nx*grid.Ny*grid.Nz,3);
  
  % update the point loctions
  X = Xn + dt*Uw;
  
