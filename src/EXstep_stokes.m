function [X,Uw,u,E] = EXstep_stokes(Xn,dt,fbhat_ext,ks,kb,kappa,grid);
     
  % form the spreading operator
  %
  spfactor = grid.ds/(grid.dx*grid.dx);
  Sm = spreadmatrix_vc_vec(Xn,grid.dx,grid.Nx,grid.Ny,grid.xmin,grid.ymin);
    
  % Evaluate the forces at the current position
  %
  [Fb,Kx] = bending_force_vec(Xn,kappa,kb,grid.ds);  
  [Fs,St] = stretch_force_vec(Xn,ks,grid.ds);
  F  = Fb + Fs;
  
  % compute the elastic energy of the worm
  %
  E = kb*sum((Kx(2:end-1)-kappa(2:end-1)).^2) + ks*sum( (St(1:end-1)-1).^2);
  E = grid.ds*0.5*E;
   
  % spread forces
  %
  fb = spfactor * reshape(Sm*F,grid.Nx,grid.Ny,2);
  
  % compute total force in fourier space
  %
  fbhat = fft2( fb ) + fbhat_ext;
  
  % solve stokes in fourier space
  %
  [uhat,phat]=stokes_solve_fourier(fbhat,grid.Lx,grid.Ly);
    
  % transform velocity back to real space
  %
  u = real( ifft2(uhat) );
    
  % interpolate velocity back to the IB points
  %
  Uw = Sm'*reshape(u,grid.Nx*grid.Ny,2);
  
  % update the point loctions
  %
  X = Xn + dt*Uw;
  
