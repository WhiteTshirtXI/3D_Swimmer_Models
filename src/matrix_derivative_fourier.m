function [fphat]=matrix_derivative_fourier(fhat,Lx,Ly)

% record the number of grid points in each direction
%
szf=size(fhat);
nx=szf(1);
ny=szf(2);
if ndims(fhat)<3
    nd=1;
else
    nd=szf(3);
end

fphat = zeros(nx,ny,2*nd);
% compute the wave numbers
%
N1x =  floor((nx-1)/2);
N2x = (nx/2)*ones(rem(nx+1,2));
freqx =(2*pi/Lx)* [(0:N1x)  N2x (-N1x:-1)]';

N1y =  floor((ny-1)/2);
N2y = (ny/2)*ones(rem(ny+1,2));
freqy = (2*pi/Ly)*[(0:N1y)  N2y (-N1y:-1)]';

[k1, k2]=ndgrid(freqx,freqy);

i=sqrt(-1);

cd=1;
for ll=1:2:2*nd
    fphat(:,:,ll)=i*k1.*fhat(:,:,cd);
    fphat(:,:,ll+1)=i*k2.*fhat(:,:,cd);
    cd=cd+1;
end

