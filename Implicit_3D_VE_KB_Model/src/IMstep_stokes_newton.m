function [X,Uw,U,Uhat,output] = IMstep_stokes_newton(Xn,dt,fbhat,ks,kb,kappa,grid,rtol,rXtol,gmrestol,lam,xi)
   
    gmresiter   = 20;
    gmresrstart = 20;
    maxiter     = 75;
    
    % form the spreading operator
    %
    spfactor = grid.ds/(grid.dx*grid.dx*grid.dx);
    Sm = spreadmatrix3_vc_vec(Xn,grid.dx,grid.Nx,grid.Ny,grid.Nz,grid.xmin,grid.ymin,grid.zmin);
    
    % record the number of IB points
    %
    N = size(Xn,1);

    % initialize and begin solver loop
    %
    X = Xn;
    jac_calls = 0;
    output.fncalls = 1;
    output.iter    = 0;
    output.gmiter  = 0;
    
    % make an initial funciton call
    %
    [G,Uw,U] = IMfun(X,Xn,dt,Sm,spfactor,ks,kb,kappa,fbhat,grid,lam,xi);

    
    % for a Stokelets matrix for preconditioning
    %
    epsilon = 1.5*grid.ds;
    mu = 1.0;
    M = form_reg_stokes_matrix_3D(Xn,epsilon,mu);
    M = grid.ds*M;
 
    
    % begin main loop
    %
    for k=1:maxiter
            
        % compute the force Jacobians
        %
        Jb = bend_force_jac3(X,kappa,kb,grid.ds);
        Js = stretch_force_jac3(X,ks,grid.ds);
        JF = Jb + Js;

        % make approximate J for preconditioning
        %
        %J = eye(2*N,2*N) - dt*SSblock*JF;
        J = speye(3*N,3*N) - dt*M*JF;
        
        % invert Jacobian using J as preconditioner
        %
        A = @(Y)JMfun(Y,JF,Sm,spfactor,dt,grid);
        [dX,flag,relres,iter] = gmres(A,-G(:),20,gmrestol,20,J);
        dX = reshape(dX,N,3);
        
        % record the number of function calls
        %
        output.iter    = output.iter + 1;
        output.gmiter  = output.gmiter + (iter(1)-1)*gmresrstart + iter(2);
        output.fncalls = output.fncalls + 1 + (iter(1)-1)*gmresrstart + iter(2);
        
        % update X
        %
        X  = X + dX;
         
        % eval function
        %
        [G,Uw,U,Uhat] = IMfun(X,Xn,dt,Sm,spfactor,ks,kb,kappa,fbhat,grid,lam,xi);
 
        % record the size of G 
        %
        G2 = sqrt(mean( sum(G.^2,2)));
        Gmax = max(abs(G(:)));
        G2r = G2./sqrt( mean( sum(Xn.^2,2)));
        
        % record the size of the updates
        %
        dX2 = norm(dX(:));
        dX2r = dX2/norm(Xn(:));
        
        fprintf('iter = %i;  abs norm X = %8.2e; rel norm X = %8.2e; max(G)=%8.2e; G2=%8.2e\n',...
                k,dX2,dX2r,Gmax,G2);

        % check stopping condition 
        %
        if( (G2r < rtol) |  (dX2r < rXtol) )
            fprintf('JFNK converged in %i steps \n',k);
            fprintf('  total gmres iterations %i \n',output.gmiter);
            fprintf('  total function calls   %i \n',output.fncalls);
            fprintf('  rel(dX) = %8.2e, abs(dX) = %8.2e \n',dX2r,dX2);
            fprintf('  rel(G2) = %8.2e; abs(G2)  = %8.2e \n',G2r,G2);
            fprintf('  max(G2) = %8.2e \n\n',Gmax);
            break
        end
          
    end
        
    
    % did not converge, print warning
    %
    if( k==maxiter )
        fprintf('WARNING: JFNK iteration did not converge \n');
        fprintf('  max iterations = %g \n',maxiter);
        fprintf('  Xtol = %8.2e, rel(dX) = %8.2e \n',rXtol,dX2r);
        fprintf('  Gtol = %8.2e, G2 = %8.2e, Gmax = %8.2e \n',rtol,G2,Gmax);
        fprintf('  attempting to continue the simulation \n');
    end
        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
function W = JMfun(Y,JF,Sm,spfactor,dt,grid);
%
% function to apply the Jacobian of IMfun to the vector Y
%


  % multiply by force Jacobian
  % 
  Z = JF*Y;
   
  % reshape before spreading
  %
  F  = reshape(Z,grid.N,3);
   
  % spread and reshape to Eulerian grid 
  %
  fb = spfactor * reshape(Sm*F,grid.Nx,grid.Ny,grid.Nz,3); 
   
  % solve Stokes equations
  %
  fbhat = zeros(grid.Nx, grid.Ny, grid.Nz, 3);
  for d = 1:3
      fbhat(:,:,:,d) = fftn(fb(:,:,:,d));
  end
  uhat = stokes_solve_fourier_3d(fbhat,grid.Lx,grid.Ly,grid.Lz);
  
  U = zeros(grid.Nx, grid.Ny, grid.Nz, 3);
  for i = 1:3
      U(:,:,:,i) = real ( ifftn ( uhat(:,:,:,i)));
  end
  
  % interpolate back to the IB points
  %
  Uw = Sm'*reshape(U,grid.Nx*grid.Ny*grid.Nz,3);   
  
  % finish application of Jacobian
  %
  W = Y(:) - dt*Uw(:);

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
function [G,Uw,U,Uhat] = IMfun(X,Xn,dt,Sm,spfactor,ks,kb,kappa,fbhat_ext,grid,lam,xi)
%
% Use this function to solve G= X(n+1)-X(n)-dt*U = 0
%

   % Evaluate the forces at the current position
   %
   [Fb,Kx] = bending_force_vec3(X,kappa,kb,grid.ds);  
   [Fs,St] = stretch_force_vec3(X,ks,grid.ds);
   F  = Fb + Fs;
    
   % spread forces
   %
   fb = spfactor * reshape(Sm*F,grid.Nx,grid.Ny,grid.Nz,3);
    
   % solve stokes
   %
   fbhat = zeros(grid.Nx, grid.Ny, grid.Nz, 3);
   for d = 1:3
       fbhat(:,:,:,d) = fftn(fb(:,:,:,d));
   end
   
   fbhat = fbhat + fbhat_ext;
   Uhat = stokes_solve_fourier_3d(fbhat,grid.Lx,grid.Ly,grid.Lz);
   
   U = zeros(grid.Nx, grid.Ny, grid.Nz, 3);
   for i = 1:3
       U(:,:,:,i) = real( ifftn( Uhat(:,:,:,i)));
   end
   % Newtonian Swimmer
   if(lam == 0)
       U = U*(1/(1+xi));
   end
   
   % interpolate back to the IB points
   %
   Uw = Sm'*reshape(U,grid.Nx*grid.Ny*grid.Nz,3);
   
   % eval function 
   %
   G = (X - dt*Uw) - Xn;
