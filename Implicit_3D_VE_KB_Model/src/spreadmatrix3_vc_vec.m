%
% spreadmatrix_vc.m
%
% compute the matrix for the spreading operator
%   discretized by a vertex centered grid with N points in each
%   direction 
%
% this is coded for Periodic boundary conditions 
%
% input,  X -- matrix with ib point locations (size Nib x 3 )
%         N -- number of grid points in each direction (does not include
%              the boundary points) 
%
% output, S -- scaled spreading operator of size N*N*N x Nib
%
%
function S = spreadmatrix3_vc_vec(X,dx,Nx,Ny,Nz,xmin,ymin,zmin)
  
  % record the number of unknowns
  %
  Nib = size(X,1);

  % allocate space for S
  %
  N3 = Nx*Ny*Nz;
  %S = spalloc( N3, Nib, 64*Nib);       %this was giving me problems
  
  % convert X to grid coordinates (xg,yg,zg)
  %
  xg = (X(:,1)-xmin)/dx  + 1;
  yg = (X(:,2)-ymin)/dx  + 1;
  zg = (X(:,3)-zmin)/dx  + 1;
  
  % indices of grid point down and to the left
  %
  I0 = floor( xg );
  J0 = floor( yg );
  K0 = floor( zg );

  % compute shifts of the indices
  %
  Im = I0-1;
  I1 = I0+1;
  I2 = I0+2;
  
  Jm = J0-1;
  J1 = J0+1;
  J2 = J0+2;

  Km = K0-1;
  K1 = K0+1;
  K2 = K0+2;

  % compute the weights
  %
  Wxm = delta(Im - xg);
  Wx0 = delta(I0 - xg);
  Wx1 = delta(I1 - xg);
  Wx2 = delta(I2 - xg);
  
  Wym = delta(Jm - yg);
  Wy0 = delta(J0 - yg);
  Wy1 = delta(J1 - yg);
  Wy2 = delta(J2 - yg);
  
  Wzm = delta(Km - zg);
  Wz0 = delta(K0 - zg);
  Wz1 = delta(K1 - zg);
  Wz2 = delta(K2 - zg);

  % assemble these into single arrrays with four columns
  %
  Wx = [Wxm Wx0 Wx1 Wx2];
  Wy = [Wym Wy0 Wy1 Wy2];
  Wz = [Wzm Wz0 Wz1 Wz2];

  
  % done computing weights, make I, J, K periodic
  %
  Im = mod(Im-1,Nx) + 1;
  I0 = mod(I0-1,Nx) + 1;
  I1 = mod(I1-1,Nx) + 1;
  I2 = mod(I2-1,Nx) + 1;
 
  Jm = mod(Jm-1,Ny) + 1;
  J0 = mod(J0-1,Ny) + 1;
  J1 = mod(J1-1,Ny) + 1;
  J2 = mod(J2-1,Ny) + 1;
   
  Km = mod(Km-1,Nz) + 1;
  K0 = mod(K0-1,Nz) + 1;
  K1 = mod(K1-1,Nz) + 1;
  K2 = mod(K2-1,Nz) + 1;
  
  % assemble into single arrays
  %
  II = [Im I0 I1 I2];
  JJ = [Jm J0 J1 J2];
  KK = [Km K0 K1 K2];

   
  % make Iv, Jv, Kv 4d arrays
  %
  Iv = reshape(II,[Nib 4 1 1]);
  Jv = reshape(JJ,[Nib 1 4 1]);
  Kv = reshape(KK,[Nib 1 1 4]);

  v = [1 1 1 1];
  Iv = Iv(:,:,v,v);
  Jv = Jv(:,v,:,v);
  Kv = Kv(:,v,v,:);

  % same trick for the W's
  %
  Wx = reshape(Wx,[Nib 4 1 1]);
  Wy = reshape(Wy,[Nib 1 4 1]);
  Wz = reshape(Wz,[Nib 1 1 4]);
  
  Wx = Wx(:,:,v,v);
  Wy = Wy(:,v,:,v);
  Wz = Wz(:,v,v,:);
  W = Wx.*Wy.*Wz;
  
  % column numbers 
  %
  Kc = repmat( (1:Nib)',64,1);
  Kc = Kc(:);
  
  % row numbers
  %
  Kr = sub2ind([Nx, Ny, Nz],Iv(:),Jv(:),Kv(:));
  
  
  % make a sparse matrix of the weights
  %
  S = sparse(Kr(:),Kc(:),W(:),N3,Nib);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% delta -- form of the 4-point discrete delta function
%
function phi = delta(r)
%  phi = 0.25*( 1.0 + cos(0.5*pi*r));

    ra = abs(r);
    phi1 = 0.125*( 3 - 2*ra + sqrt(1 + 4*ra- 4*r.^2));
    phi2 = 0.125*( 5 - 2*ra - sqrt(-7+12*ra- 4*r.^2));
  
    phi = phi1.*double( ra < 1 ) + phi2.*double( ra>=1 );
    phi = phi.*double(ra<2);
