function [F,Kx] = bending_force(X,kappa,kb,ds);
  
  % record the number of points
  %
  N = size(X,1);
  
  % initialize the forces
  %
  F = zeros(N,2);
  
  
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
  
  
  % loop over the interior points and update forces
  %
  for j=2:N-1

   F(j-1,:) = F(j-1,:) - W(j)*[  Dp(j,2)        , -Dp(j,1)        ];
   F(j  ,:) = F(j  ,:) - W(j)*[  -Dm(j,2)-Dp(j,2),  Dp(j,1)+Dm(j,1)];
   F(j+1,:) = F(j+1,:) - W(j)*[  Dm(j,2)        ,  -Dm(j,1)        ];
       
  end
    
  % rescale the forces
  %
  F = kb * F/ds^3;
