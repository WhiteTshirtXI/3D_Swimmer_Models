function X = initialize_worm2(X0,kappa,ds);
    
  % initialize location
  %
  X = X0;
      
      
  tol  = 1e-4;
  options = optimset('Jacobian','on','Display','off');

  % some stiffness parameters -- just used to find shape
  %
  kb=1;
  ks=2500;
  
  % damping parameter
  %
  dt = 1;    
  
  fprintf('Initializing worm location \n');
  for j=1:1000
      [Z,fval,flag,output] = fsolve(@(Y)IMstep(Y,X,dt,ds,ks,kb,kappa),X, ...
                                    options);
      
      
      
      if( flag <= 0 )
          fprintf('  initialization solver failed, reducing damping from %e to %e \n',dt,dt/10);
          dt = dt/2;
      else 
          X = Z;
          [Fb,Kx] = bending_force_vec(X,kappa,kb,ds);  
          err = max(abs(Kx(2:end-1)-kappa(2:end-1))) /max(abs(kappa));
          fprintf('  iter = %i; err = %e; flag=%g \n',j,err,flag);
          if( err < tol | err < 1e-10 )
              break
          end
      end
          
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
function [G,J] = IMstep(X,Xn,dt,ds,ks,kb,kappa);
    [Fb,Kx] = bending_force_vec(X,kappa,kb,ds);  
    [Fs,St] = stretch_force_vec(X,ks,ds);
    F  = Fb + Fs;
    
    G = (X - dt*F) - Xn;
    
    N = size(X,1);
    Jb = bend_force_jac(X,kappa,kb,ds);
    Js = stretch_force_jac(X,ks,ds);
    J = eye(2*N,2*N) - dt*(Js+Jb);
    
    
