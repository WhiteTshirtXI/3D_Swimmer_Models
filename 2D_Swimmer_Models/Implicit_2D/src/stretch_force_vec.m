function [F,St] = stretch_force_vec(X,ks,ds);
  
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

  % compute the current lengths
  %
  L = sqrt( sum(Dp.^2,2));
  
  % update the forces spring-by-spring
  %
  J=1:N-1;
  T = repmat( (L(J)-ds)./L(J),1,2).*Dp(J,:);
  F(J  ,:) = F(J  ,:) + T;
  F(J+1,:) = F(J+1,:) - T;
   
  % size of |X_{s}|
  %
  St = L/ds;
  
  
  % rescale forces
  %
  F = ks/ds^2 * F;
