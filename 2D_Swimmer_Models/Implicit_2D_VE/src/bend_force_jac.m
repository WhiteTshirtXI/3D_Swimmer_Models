%
% bend_force_jac -- compute the Jacobian of the bending
%
function J = bend_force_jac(X,kappa,kb,ds)
    
  % record the number of points
  %
  N = size(X,1);

  % initialize the sub plots of the Jacobians to be zero
  %
  J11 = zeros(N,N);
  J12 = J11;
  J21 = J11;
  J22 = J11;
  
 
  % compute the differences of the point location
  %  note that D(i) is the forward difference for point i
  %
  D = X(2:N,:) - X(1:N-1,:);
  Dp = [D; [0 0]];
  Dm = [[0 0]; D];

  % compute the energy density at each point
  %
  Kx = zeros(N,1);
  W = zeros(N,1);
  K = 2:N-1;
  Kx(K) =  (Dp(K,1).*Dm(K,2) - Dp(K,2).*Dm(K,1))/ds^3;
  W(K) = Kx(K) - kappa(K);
  
  % derivatives of W
  %
  dWdXim1 =  Dp(:,2)/ds^3;
  dWdXi   = -(Dp(:,2)+ Dm(:,2))/ds^3;
  dWdXip1 =  Dm(:,2)/ds^3;

  dWdYim1 = -Dp(:,1)/ds^3;
  dWdYi   =  (Dp(:,1)+ Dm(:,1))/ds^3;
  dWdYip1 = -Dm(:,1)/ds^3;
 

  
  
  % loop over the interior points
  %
  for i=2:N-1
      
      J11(i-1,i-1) = J11(i-1,i-1) - dWdXim1(i)*dWdXim1(i); %
      J11(i-1,i  ) = J11(i-1,i  ) - dWdXi(i)  *dWdXim1(i); %
      J11(i-1,i+1) = J11(i-1,i+1) - dWdXip1(i)*dWdXim1(i); %
      
      J21(i-1,i-1) = J21(i-1,i-1) - dWdXim1(i)*dWdYim1(i); %
      J21(i-1,i  ) = J21(i-1,i  ) - dWdXi(i)  *dWdYim1(i) + W(i)/ds^3; %%????
      J21(i-1,i+1) = J21(i-1,i+1) - dWdXip1(i)*dWdYim1(i) - W(i)/ds^3;
      
      J12(i-1,i-1) = J12(i-1,i-1) - dWdYim1(i)*dWdXim1(i); %
      J12(i-1,i  ) = J12(i-1,i  ) - dWdYi(i)  *dWdXim1(i) - W(i)/ds^3;
      J12(i-1,i+1) = J12(i-1,i+1) - dWdYip1(i)*dWdXim1(i) + W(i)/ds^3;

      J22(i-1,i-1) = J22(i-1,i-1) - dWdYim1(i)*dWdYim1(i); %
      J22(i-1,i  ) = J22(i-1,i  ) - dWdYi(i)  *dWdYim1(i); %
      J22(i-1,i+1) = J22(i-1,i+1) - dWdYip1(i)*dWdYim1(i); %
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
      J11(i,i-1) = J11(i,i-1) - dWdXim1(i)*dWdXi(i); %
      J11(i,i  ) = J11(i,i  ) - dWdXi(i)  *dWdXi(i); %
      J11(i,i+1) = J11(i,i+1) - dWdXip1(i)*dWdXi(i); %
      
      J21(i,i-1) = J21(i,i-1) - dWdXim1(i)*dWdYi(i) - W(i)/ds^3;
      J21(i,i  ) = J21(i,i  ) - dWdXi(i)  *dWdYi(i); %
      J21(i,i+1) = J21(i,i+1) - dWdXip1(i)*dWdYi(i) + W(i)/ds^3;
      
      J12(i,i-1) = J12(i,i-1) - dWdYim1(i)*dWdXi(i) + W(i)/ds^3; 
      J12(i,i  ) = J12(i,i  ) - dWdYi(i)  *dWdXi(i); %
      J12(i,i+1) = J12(i,i+1) - dWdYip1(i)*dWdXi(i) - W(i)/ds^3;

      J22(i,i-1) = J22(i,i-1) - dWdYim1(i)*dWdYi(i); %
      J22(i,i  ) = J22(i,i  ) - dWdYi(i)  *dWdYi(i); %
      J22(i,i+1) = J22(i,i+1) - dWdYip1(i)*dWdYi(i); %

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      J11(i+1,i-1) = J11(i+1,i-1) - dWdXim1(i)*dWdXip1(i); %
      J11(i+1,i  ) = J11(i+1,i  ) - dWdXi(i)  *dWdXip1(i); %
      J11(i+1,i+1) = J11(i+1,i+1) - dWdXip1(i)*dWdXip1(i); %
      
      J21(i+1,i-1) = J21(i+1,i-1) - dWdXim1(i)*dWdYip1(i) + W(i)/ds^3; 
      J21(i+1,i  ) = J21(i+1,i  ) - dWdXi(i)  *dWdYip1(i) - W(i)/ds^3; 
      J21(i+1,i+1) = J21(i+1,i+1) - dWdXip1(i)*dWdYip1(i);   %          
      
      J12(i+1,i-1) = J12(i+1,i-1) - dWdYim1(i)*dWdXip1(i) - W(i)/ds^3;
      J12(i+1,i  ) = J12(i+1,i  ) - dWdYi(i)  *dWdXip1(i) + W(i)/ds^3;
      J12(i+1,i+1) = J12(i+1,i+1) - dWdYip1(i)*dWdXip1(i); %

      J22(i+1,i-1) = J22(i+1,i-1) - dWdYim1(i)*dWdYip1(i); %
      J22(i+1,i  ) = J22(i+1,i  ) - dWdYi(i)  *dWdYip1(i); %
      J22(i+1,i+1) = J22(i+1,i+1) - dWdYip1(i)*dWdYip1(i); %
       
  end
  
  % put the blocks into one big matrix and rescale
  %
  J = kb*[ [J11, J12]; [J21, J22]];
