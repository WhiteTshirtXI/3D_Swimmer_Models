function M = form_reg_stokes_matrix_3D(X,epsilon,mu);
  
  Nb = size(X,1);
  Xm = X(:,ones(1,Nb));
  Ym = X(:,2*ones(1,Nb));
  Zm = X(:,3*ones(1,Nb));
  
  XX = (Xm-Xm').^2;
  YY = (Ym-Ym').^2;
  ZZ = (Zm-Zm').^2;
 
  XY = (Xm-Xm').*(Ym-Ym');
  XZ = (Xm-Xm').*(Zm-Zm');
  YZ = (Ym-Ym').*(Zm-Zm');
  
  
  R    = sqrt( XX + YY + ZZ);
  Re   = sqrt( R.^2 + epsilon.^2);
  
  P2 = 1./Re.^3;
  P1 = (R.*R + 2*epsilon^2).*P2;
 
  
  
  M = [ [P1 + P2.*XX,      P2.*XY,      P2.*XZ]; 
        [     P2.*XY, P1 + P2.*YY,      P2.*YZ];
	[     P2.*XZ,      P2.*YZ, P1 + P2.*ZZ];
      ];
  M = M/(8*pi*mu);
