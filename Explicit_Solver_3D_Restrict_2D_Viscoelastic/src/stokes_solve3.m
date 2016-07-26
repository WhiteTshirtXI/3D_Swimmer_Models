%
% solve stokes equations on periodic domain (0,Lx)x(0,Ly)x(0,Lz)
%
function [u,p,uhat,phat] = stokes_solve3(f,Lx,Ly,Lz);
  
   % record the number of grid points
   Nxyz = size(f);
   Lxyz=[Lx,Ly,Lz];
   
   for d=1:3
     N  = Nxyz(d);
     N1 = floor( (N-1)/2 );
     N2 = N/2*ones(rem(N+1,2));
     freq(d,:) = 2*pi/Lxyz(d) *[(0:N1)  N2 (-N1:-1)]';
   end
   [kx ky kz]=ndgrid(freq(1,:),freq(2,:),freq(3,:));   
   ksq = kx.^2 + ky.^2 + kz.^2;
   i=sqrt(-1);
   
   % avoid division by zero
   ksq(1,1,1) = 1.0;
   
   % take fourier transform of the force
   fhat = zeros(size(f));
   for d=1:3
     fhat(:,:,:,d) = fftn(f(:,:,:,d));
   end
         
   % take the divergence of the force
   divf = i*(kx.*fhat(:,:,:,1) + ky.*fhat(:,:,:,2) + kz.*fhat(:,:,:,3));
   
   % solve for the pressure
   phat = -divf./ksq;

   % solve for the velocity
   uhat(:,:,:,1) = (fhat(:,:,:,1) - i*kx.*phat)./ksq;
   uhat(:,:,:,2) = (fhat(:,:,:,2) - i*ky.*phat)./ksq;
   uhat(:,:,:,3) = (fhat(:,:,:,3) - i*kz.*phat)./ksq;
      
   % force mean zero
   uhat(1,1,1,:) = 0.0;
   
   % transform back
   for d=1:3
     u(:,:,:,d) = real(ifftn(uhat(:,:,:,d)));
   end
   
   p = real(ifftn(phat));