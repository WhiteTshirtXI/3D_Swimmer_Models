function Ffil=hfil(F,Lx,Ly)

% Compute the Hou-filter to cut-off high fourier modes

szf=size(F);
nx=szf(1);
ny=szf(2);
if ndims(F)<3
    nd=1;
else
    nd=szf(3);
end

% compute the wave numbers
%
N1x =  floor((nx-1)/2);
N2x = (nx/2)*ones(rem(nx+1,2));
freqx =(2*pi/Lx)* [(0:N1x)  N2x (-N1x:-1)]';

N1y =  floor((ny-1)/2);
N2y = (ny/2)*ones(rem(ny+1,2));
freqy = (2*pi/Ly)*[(0:N1y)  N2y (-N1y:-1)]';

[k1, k2]=ndgrid(freqx,freqy);

cutoff=exp(-36*(sqrt((k1/(nx/2)).^2+(k2/(ny/2)).^2)).^36);

for i=1:nd
    fftcutoff(:,:,i)= cutoff;
end

Ffil=fftcutoff.*F;

