function [uhat,phat]=stokes_solve_fourier(fhat,Lx,Ly)

% Solve \lap u -\grad p + f = 0

% record the number of grid points in each direction
%
szf=size(fhat);
nx=szf(1);
ny=szf(2);

% compute the wave numbers
%
N1x =  floor((nx-1)/2);
N2x = (nx/2)*ones(rem(nx+1,2));
freqx =(2*pi/Lx)* [(0:N1x)  N2x (-N1x:-1)]';

N1y =  floor((ny-1)/2);
N2y = (ny/2)*ones(rem(ny+1,2));
freqy = (2*pi/Ly)*[(0:N1y)  N2y (-N1y:-1)]';

[k1, k2]=ndgrid(freqx,freqy);


% compute k.k and set the zero mode to 1 to avoid division by zero
%
ksq=k1.^2+k2.^2;
ksq(1,1)=1;

% solve for pressure
%
i=sqrt(-1);
divf_hat=i*k1.*fhat(:,:,1)+i*k2.*fhat(:,:,2);
phat=(-1./ksq).*(divf_hat);


% solve for velocity
%
uhat = zeros(nx,ny,2);
uhat(:,:,1)=(fhat(:,:,1)-i*k1.*phat)./ksq;
uhat(:,:,2)=(fhat(:,:,2)-i*k2.*phat)./ksq;

