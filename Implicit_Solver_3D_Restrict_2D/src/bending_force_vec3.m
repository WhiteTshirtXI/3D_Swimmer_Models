function [F,Kx] = bending_force_vec3(X,kappa,kb,ds);
  
  % record the number of points
  %
  N = size(X,1);
  
  % initialize the forces
  %
  F = zeros(N,3);
  
  
  % compute the differences of the point location
  %  note that D(i) is the forward difference for point i
  %
  D = X(2:N,:) - X(1:N-1,:);
  Dp = [D; [0 0 0]];
  Dm = [[0 0 0]; D];
  
  % compute the energy density at each point
  %
  Kx = zeros(N,1);
  W = zeros(N,1);
  K = 2:N-1;
  Kx(K) =  (Dp(K,1).*Dm(K,2) - Dp(K,2).*Dm(K,1))/ds^3;
  W(K) = Kx(K) - kappa(K);
  
  % update the bending forces
  %
  J = 2:N-1;
  WW = repmat(W,1,3);
  A(J,1) = 0; 
  F(J-1,:) = F(J-1,:) - WW(J,:).*[  Dp(J,2)        , -Dp(J,1)    , A(J)    ];
  F(J  ,:) = F(J  ,:) - WW(J,:).*[ -Dm(J,2)-Dp(J,2),  Dp(J,1)+Dm(J,1), A(J)];
  F(J+1,:) = F(J+1,:) - WW(J,:).*[  Dm(J,2)        ,  -Dm(J,1)   , A(J)   ];
 
  
% $$$   
% $$$   
% $$$   % loop over the interior points and update forces
% $$$   %
% $$$   G = zeros(N,2);
% $$$   for j=2:N-1
% $$$    G(j-1,:) = G(j-1,:) - W(j)*[  Dp(j,2)        , -Dp(j,1)        ];
% $$$    G(j  ,:) = G(j  ,:) - W(j)*[  -Dm(j,2)-Dp(j,2),  Dp(j,1)+Dm(j,1)];
% $$$    G(j+1,:) = G(j+1,:) - W(j)*[  Dm(j,2)        ,  -Dm(j,1)        ];
% $$$   end
% $$$   fprintf('max diff of F and G %g \n',max( abs(F(:)-G(:))));  
% $$$   
  
  % rescale the forces
  %
  F = kb * F/ds^3;
