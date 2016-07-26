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
%  this code assumes for a periodic [0,Lx]^3 grid
%
function S = spreadmatrix3_vc(X,Nx,Ny,Nz,dx)
  
  % record the number of unknowns
  %
  Nib = size(X,1);

  % allocate space for S
  %
  N3 = Nx*Ny*Nz;
  S = spalloc( N3, Nib, 64*Nib);

  
  % loop over the IB points
  %
  for n=1:Nib
    n
    % (xg,yg,zg) is the scaled IB loction 
    %
    xg = X(n,1)/dx + 1;
    yg = X(n,2)/dx + 1;
    zg = X(n,3)/dx + 1;
    
    % grid point down and to the left
    %
    i0 = floor( xg );
    j0 = floor( yg );
    k0 = floor( zg );
    
    
    % adjust for peridicity
    %
    i0 = mod(i0-1,Nx) + 1;
    j0 = mod(j0-1,Ny) + 1;
    k0 = mod(k0-1,Nz) + 1;
    
    
    % compute the weights
    %
    for i = i0-1:i0+2
      
      % evalute x-weight
      %
      wx = delta(i-xg);
      
      % adjust for periodic
      %
      ii = mod(i-1,Nx) + 1;
      
      for j=j0-1:j0+2
        
        % evaluate y-weight
        %
        wy = delta(j-yg);
    
        % adjust for periodic
        %
        jj = mod(j-1,Ny) + 1;
             
        for k=k0-1:k0+2

            % adjust for periodic
            %
            wz = delta(k-zg);
            
            % adjust for periodic
            %
            kk = mod(k-1,Nz) + 1;
                        
            % row of S -- assuming standard row ordering
            %
            kr = sub2ind([Nx,Ny,Nz],ii,jj,kk);
                        
            % update S
            %
            S(kr,n) = S(kr,n) + wx*wy*wz;
	  
	  end % end loop over k
      end   % end loop over j 
    end     % end loop over i

  end       % end loop over IB points
 
  
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


  
  
  
  
  
  
