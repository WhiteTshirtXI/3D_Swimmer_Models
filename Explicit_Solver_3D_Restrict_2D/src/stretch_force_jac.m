%
% stretch_force_jac -- compute the Jacobian of the stretching force
%
function J = stretch_force_jac(X,ks,ds)
    
  % record the number of points
  %
  N = size(X,1);

 
  % compute the differences of the point location
  %  note that D(i) is the forward difference for point i
  %
  D = X(2:N,:) - X(1:N-1,:);
  Dp = [D; [0 0]];
  Dm = [[0 0]; D];

  % compute the current lengths
  %
  L = sqrt( sum(Dp.^2,2));
  A = 1 - ds./L;
  
  
  % initialize the sub plots of the Jacobians to be zero
  %
  J11 = zeros(N,N);
  J12 = J11;
  J21 = J11;
  J22 = J11;
  
  % loop over the springs
  %
  for i=1:N-1
      
      L3 = (L(i)).^3;
      
      J11(i,i  ) = J11(i,i  ) - A(i) - ds*Dp(i,1)^2/L3;
      J11(i,i+1) = J11(i,i+1) + A(i) + ds*Dp(i,1)^2/L3;
      
      J12(i,i  ) = J12(i,i  ) - ds*Dp(i,1)*Dp(i,2)/L3;
      J12(i,i+1) = J12(i,i+1) + ds*Dp(i,1)*Dp(i,2)/L3;
      
      J21(i,i  ) = J21(i,i  ) - ds*Dp(i,1)*Dp(i,2)/L3;
      J21(i,i+1) = J21(i,i+1) + ds*Dp(i,1)*Dp(i,2)/L3;
      
      J22(i,i  ) = J22(i,i  ) - A(i) - ds*Dp(i,2)^2/L3;
      J22(i,i+1) = J22(i,i+1) + A(i) + ds*Dp(i,2)^2/L3;
  

      J11(i+1,i  ) = J11(i+1,i  ) + A(i) + ds*Dp(i,1)^2/L3;
      J11(i+1,i+1) = J11(i+1,i+1) - A(i) - ds*Dp(i,1)^2/L3;
      
      J12(i+1,i  ) = J12(i+1,i  ) + ds*Dp(i,1)*Dp(i,2)/L3;
      J12(i+1,i+1) = J12(i+1,i+1) - ds*Dp(i,1)*Dp(i,2)/L3;
      
      J21(i+1,i  ) = J21(i+1,i  ) + ds*Dp(i,1)*Dp(i,2)/L3;
      J21(i+1,i+1) = J21(i+1,i+1) - ds*Dp(i,1)*Dp(i,2)/L3;
      
      J22(i+1,i  ) = J22(i+1,i  ) + A(i) + ds*Dp(i,2)^2/L3;
      J22(i+1,i+1) = J22(i+1,i+1) - A(i) - ds*Dp(i,2)^2/L3;
       
  end
  
  % put the blocks into one big matrix and rescale
  %
  J = [ [J11, J12]; [J21, J22]];
  J = ks*J/ds^2;