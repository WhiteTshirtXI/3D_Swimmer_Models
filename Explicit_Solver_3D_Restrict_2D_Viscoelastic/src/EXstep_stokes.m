function [X,Uw,u,E,uhat] = EXstep_stokes(Xn,dt,fbhat_ext,ks,kb,kappa,grid)
     
  % form the spreading operator
  spfactor = grid.ds/(grid.dx*grid.dx*grid.dx);
  Sm = spreadmatrix3_vc(Xn,grid.Nx,grid.Ny,grid.Nz,grid.dx);
    
  % Evaluate the forces at the current position
  [Fb,Kx] = bending_force3(Xn,kappa,kb,grid.ds);  
  [Fs,St] = stretch_force3(Xn,ks,grid.ds);
  F  = Fb + Fs;
  
  % compute the elastic energy of the worm
  E = kb*sum((Kx(2:end-1)-kappa(2:end-1)).^2) + ks*sum( (St(1:end-1)-1).^2);
  E = grid.ds*0.5*E;
   
  % spread forces
  fb = spfactor * reshape(Sm*F,grid.Nx,grid.Ny,grid.Nz,3);
  
  % compute total force in fourier space
  for d = 1:3
      fbhat(:,:,:,d) = fftn(fb(:,:,:,d));
  end
  fbhat = fbhat + fbhat_ext;
  
  % solve stokes in fourier space
  [u,p,uhat,phat] = stokes_solve3(fbhat,grid.Lx,grid.Ly,grid.Lz);
    
  % interpolate velocity back to the IB points
  Uw = Sm'*reshape(u,grid.Nx*grid.Ny*grid.Nz,3);
  
  % update the point loctions
  X = Xn + dt*Uw;
  
