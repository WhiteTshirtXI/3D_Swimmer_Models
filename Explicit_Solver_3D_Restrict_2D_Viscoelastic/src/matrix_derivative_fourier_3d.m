function fphat=matrix_derivative_fourier_3d(Fhat,Lx,Ly,Lz)

% record the number of grid points in each direction
%
szf=size(Fhat);
nx=szf(1);
ny=szf(2);
nz=szf(3);
if ndims(Fhat)<4
    nd=1;
else
    nd=szf(4);
end

fphat = zeros(nx,ny,nz,3*nd);
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

tic=1;
for ind=1:3:3*nd
    fphat(:,:,:,ind+0)=1i.*k1.*Fhat(:,:,:,tic);
    fphat(:,:,:,ind+1)=1i.*k2.*Fhat(:,:,:,tic);
    fphat(:,:,:,ind+2)=1i.*k3.*Fhat(:,:,:,tic);
    tic=tic+1;
end

