function M = form_reg_stokes_matrix(X,epsilon,mu);
  
  Nb = size(X,1);
  Xm = X(:,ones(1,Nb));
  Ym = X(:,2*ones(1,Nb));
  XX = (Xm-Xm').^2;
  YY = (Ym-Ym').^2;
  XY = (Xm-Xm').*(Ym-Ym');
  
  R    = sqrt( XX + YY );
  Re   = sqrt( R.^2 + epsilon.^2);
  P1 = log(Re+epsilon) - epsilon*(Re+2*epsilon)./(Re.*(Re+epsilon));
  P2 =                           (Re+2*epsilon)./(Re.*(Re+epsilon).^2);
  
  
  M = [ [-P1 + P2.*XX,       P2.*XY]; 
        [      P2.*XY, -P1 + P2.*YY]
      ];
  M = M/(4*pi*mu);
