function uhat=stokes_solve_fourier_3d(fhat,Lx,Ly,Lz)

%3d version
%Using the fourier transform compute the Fourier transform 
%of u from
%the equation  \lap u -grad p + f = 0


szf=size(fhat);
nx=szf(1);
ny=szf(2);
nz=szf(3);

% compute the wave numbers
%
N1x =  floor((nx-1)/2);
N2x = (nx/2)*ones(rem(nx+1,2));
freqx =(2*pi/Lx)* [(0:N1x)  N2x (-N1x:-1)]';

N1y =  floor((ny-1)/2);
N2y = (ny/2)*ones(rem(ny+1,2));
freqy = (2*pi/Ly)*[(0:N1y)  N2y (-N1y:-1)]';

N1z =  floor((nz-1)/2);
N2z = (nz/2)*ones(rem(nz+1,2));
freqz = (2*pi/Lz)*[(0:N1z)  N2z (-N1z:-1)]';

[k1, k2, k3]=ndgrid(freqx,freqy,freqz);

i=sqrt(-1);


ksq=k1.^2+k2.^2+k3.^2;

%rename ksq(1)=1 to avoid dividing by zero
ksq(1,1,1)=1;

divf_hat=i*k1.*fhat(:,:,:,1)+i*k2.*fhat(:,:,:,2)+i*k3.*fhat(:,:,:,3);
phat=(-1./ksq).*(divf_hat);


% solve for velocity
%
uhat = zeros(nx,ny,nz,3);
uhat(:,:,:,1)=(fhat(:,:,:,1)-i*k1.*phat)./ksq;
uhat(:,:,:,2)=(fhat(:,:,:,2)-i*k2.*phat)./ksq;
uhat(:,:,:,3)=(fhat(:,:,:,3)-i*k3.*phat)./ksq;


                