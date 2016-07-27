function [F,St] = stretch_force3(X,ks,ds);
  
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

  % compute the current lengths
  %
  L = sqrt( sum(Dp.^2,2));
  
  
  % loop over the springs and update the forces
  %
  for j=1:N-1
    
    F(j  ,:) = F(j,:)   +  (L(j)-ds)*Dp(j,:)/L(j);
    F(j+1,:) = F(j+1,:) -  (L(j)-ds)*Dp(j,:)/L(j);
    
  end

  % size of |X_{s}|
  %
  St = L/ds;
  
  
  % rescale forces
  %
  F = ks/ds^2 * F;
