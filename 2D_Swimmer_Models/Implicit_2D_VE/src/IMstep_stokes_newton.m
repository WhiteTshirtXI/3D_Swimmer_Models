function [X,Uw,U,output] = IMstep_stokes_newton(Xn,dt,fbhat,ks,kb,kappa,grid,rtol,rXtol,gmrestol)
   
    gmresiter   = 20;
    gmresrstart = 20;
    maxiter     = 75;
    
    % form the spreading operator
    %
    spfactor = grid.ds/(grid.dx*grid.dx);
    Sm = spreadmatrix_vc_vec(Xn,grid.dx,grid.Nx,grid.Ny,grid.xmin,grid.ymin);
    
    % record the number of IB points
    %
    N = size(Xn,1);

    % form block S'*S operator for preconditioning
    %
    SS = spfactor*Sm'*Sm;
    Zm = 0*SS;
    SSblock = [[SS, Zm];[Zm,SS]];

    % initialize and begin solver loop
    %
    X = Xn;
    jac_calls = 0;
    output.fncalls = 1;
    output.iter    = 0;
    output.gmiter  = 0;
    
    % make an initial funciton call
    %
    [G,Uw,U] = IMfun(X,Xn,dt,Sm,spfactor,ks,kb,kappa,fbhat,grid);

    
    % for a Stokelets matrix for preconditioning
    %
    epsilon = 1.5*grid.ds;
    mu = 1.0;
    M = form_reg_stokes_matrix(Xn,epsilon,mu);
    M = grid.ds*M;
 
    
    % begin main loop
    %
    for k=1:maxiter
            
        % compute the force Jacobians
        %
        Jb = bend_force_jac(X,kappa,kb,grid.ds);
        Js = stretch_force_jac(X,ks,grid.ds);
        JF = Jb + Js;

        % make approximate J for preconditioning
        %
        %J = eye(2*N,2*N) - dt*SSblock*JF;
        J = speye(2*N,2*N) - dt*M*JF;
        
        % invert Jacobian using J as preconditioner
        %
        A = @(Y)JMfun(Y,JF,Sm,spfactor,dt,grid);
        [dX,flag,relres,iter] = gmres(A,-G(:),20,gmrestol,20,J);
        dX = reshape(dX,N,2);
        
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
        [G,Uw,U] = IMfun(X,Xn,dt,Sm,spfactor,ks,kb,kappa,fbhat,grid);
 
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
    
function W = JMfun(Y,JF,Sm,spfactor,dt,grid)
%
% function to apply the Jacobian of IMfun to the vector Y
%


  % multiply by force Jacobian
  % 
  Z = JF*Y;
   
  % reshape before spreading
  %
  F  = reshape(Z,grid.N,2);
   
  % spread and reshape to Eulerian grid 
  %
  fb = spfactor * reshape(Sm*F,grid.Nx,grid.Ny,2); 
   
  % solve Stokes equations
  %
  fbhat = fft2( fb );
  [uhat,phat]=stokes_solve_fourier(fbhat,grid.Lx,grid.Ly);
  U = real( ifft2( uhat ) );
  
  
    
  % interpolate back to the IB points
  %
  Uw = Sm'*reshape(U,grid.Nx*grid.Ny,2);   
  
  % finish application of Jacobian
  %
  W = Y(:) - dt*Uw(:);

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
function [G,Uw,U] = IMfun(X,Xn,dt,Sm,spfactor,ks,kb,kappa,fbhat_ext,grid)
%
% Use this function to solve G= X(n+1)-X(n)-dt*U = 0
%

   % Evaluate the forces at the current position
   %
   [Fb,Kx] = bending_force_vec(X,kappa,kb,grid.ds);  
   [Fs,St] = stretch_force_vec(X,ks,grid.ds);
   F  = Fb + Fs;
    
   % spread forces
   %
   fb = spfactor * reshape(Sm*F,grid.Nx,grid.Ny,2);
    
   % solve stokes
   %
   fbhat = fft2( fb ) + fbhat_ext;
   [uhat,phat]=stokes_solve_fourier(fbhat,grid.Lx,grid.Ly);
   U = real( ifft2( uhat ) );
   
   %update Shat
   %[Shat,newRHS] = update_Shat(uhat,grid,Shat,nu,dt,lam,newRHS);
   %^ update stress in wrapper
    
   % interpolate back to the IB points
   %
   Uw = Sm'*reshape(U,grid.Nx*grid.Ny,2);
   
   % eval function 
   %
   G = (X - dt*Uw) - Xn;
