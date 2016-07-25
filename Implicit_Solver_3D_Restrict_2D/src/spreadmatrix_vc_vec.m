%
% spreadmatrix_vc_vec.m
%
% compute the matrix for the spreading operator
%  assumes a discretization of x = xmin + dx*(i-1); i=1..Nx
%                              y = ymin + dx*(j-1); j=1..Ny
%  it also assumes a periodic domain so that xmin == xmin + dx*Nx
%                                            ymin == ymin + dx*Ny
%
% input,  X  --   matrix with ib point locations (size Nib x 2 )
%         dx --   grid spacing
%         Nx --   number of grid points in the x-direction
%         Ny --   number of grid points in the y-direction
%         xmin -- left edge of the domain
%         ymin -- bottom edge of the domain
%
% output, S -- scaled spreading operator of size Nx*Ny x Nib
%
function S = spreadmatrix_vc_vec(X,dx,Nx,Ny,xmin,ymin);
     
  % record the number of unknowns
  %
  Nib = size(X,1);

  % record the number of grid points
  %
  Nsq = Nx*Ny;
 
  % convert X to grid coordinates (xg,yg)
  %
  xg = (X(:,1)-xmin)/dx  + 1;
  yg = (X(:,2)-ymin)/dx  + 1;
  
  % indices of grid point down and to the left
  %
  I0 = floor( xg );
  J0 = floor( yg );

  % compute shifts of the indices
  %
  Im = I0-1;
  I1 = I0+1;
  I2 = I0+2;
  
  Jm = J0-1;
  J1 = J0+1;
  J2 = J0+2;
  
  % compute the wights
  %
  Wxm = delta(Im - xg);
  Wx0 = delta(I0 - xg);
  Wx1 = delta(I1 - xg);
  Wx2 = delta(I2 - xg);

  Wym = delta(Jm - yg);
  Wy0 = delta(J0 - yg);
  Wy1 = delta(J1 - yg);
  Wy2 = delta(J2 - yg);

  
  % done computing weights, make I and J's periodic
  %
  Im = mod(Im-1,Nx) + 1;
  I0 = mod(I0-1,Nx) + 1;
  I1 = mod(I1-1,Nx) + 1;
  I2 = mod(I2-1,Nx) + 1;
 
  Jm = mod(Jm-1,Ny) + 1;
  J0 = mod(J0-1,Ny) + 1;
  J1 = mod(J1-1,Ny) + 1;
  J2 = mod(J2-1,Ny) + 1;
   

  % make four copies of each I,J corresponding to the 16 point stencil of
  % the delta function
  %
  Iv = [repmat(Im,1,4), repmat(I0,1,4), repmat(I1,1,4), repmat(I2,1,4)];
  Jv = repmat([Jm, J0, J1, J2], 1,4);
  
  
  % compute the elements of the matrix
  %
  Wx = [repmat(Wxm,1,4), repmat(Wx0,1,4), repmat(Wx1,1,4), repmat(Wx2,1,4)];
  Wy = repmat([Wym, Wy0, Wy1, Wy2], 1,4);
  W  = Wx(:).*Wy(:);
  
  % column numbers
  %
  Kc = repmat( (1:Nib)',16,1);
  Kc = Kc(:);

  % row numbers
  %
  Kr = sub2ind([Nx,Ny],Iv(:),Jv(:));
  
  % make the (scaled) spreading matrix
  %
  S = sparse(Kr(:),Kc(:),W(:), Nsq, Nib);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% delta -- form of the 4-point discrete delta function
%
function phi = delta(r);
%  phi = 0.25*( 1.0 + cos(0.5*pi*r));

    ra = abs(r);
    phi1 = 0.125*( 3 - 2*ra + sqrt(1 + 4*ra- 4*r.^2));
    phi2 = 0.125*( 5 - 2*ra - sqrt(-7+12*ra- 4*r.^2));
  
    phi = phi1.*double( ra < 1 ) + phi2.*double( ra>=1 );
    phi = phi.*double(ra<2);

  
  
  
  
  
