function  [Shat, newRHS] = update_Shat_3d(Uhat,grid,Shat,nu,dt,lam,newRHS)
% update the Stress advection equation

% Note to translate into something that reads better use the following:
% 
% u = U(:,:,1);
% v = U(:,:,2);
% w = U(:,:,3);
% 
% S11 = S(:,:,:,1);
% S12 = S(:,:,:,2);
% S13 = S(:,:,:,3);
% S22 = S(:,:,:,4);
% S23 = S(:,:,:,5);
% S33 = S(:,:,:,6);
% 
% dxu = gU(:,:,1);
% dyu = gU(:,:,2);
% dzu = gU(:,:,3);
% dxv = gU(:,:,4);
% dyv = gU(:,:,5);
% dzv = gU(:,:,6);
% dxw = gU(:,:,7);
% dyw = gU(:,:,8);
% dzw = gU(:,:,9);

% dxS11 = S(:,:,:,1);
% dyS11 = S(:,:,:,2);
% dzS11 = S(:,:,:,3);

% dxS12 = S(:,:,:,4);
% dyS12 = S(:,:,:,5);
% dzS12 = S(:,:,:,6);

% dxS13 = S(:,:,:,7);
% dyS13 = S(:,:,:,8);
% dzS13 = S(:,:,:,9);

% dxS22 = S(:,:,:,10);
% dyS22 = S(:,:,:,11);
% dzS22 = S(:,:,:,12);

% dxS23 = S(:,:,:,13);
% dyS23 = S(:,:,:,14);
% dzS23 = S(:,:,:,15);

% dxS33 = S(:,:,:,16);
% dyS33 = S(:,:,:,17);
% dzS33 = S(:,:,:,18);


% record the number of grid points in each direction
%
szf=size(Uhat);
nx=szf(1);
ny=szf(2);
nz=szf(3);
Lx=grid.Lx;
Ly=grid.Ly;
Lz=grid.Lz;
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

%update stress

gradUh=matrix_derivative_fourier_3d(Uhat,Lx,Ly,Lz);

gradSh=matrix_derivative_fourier_3d(Shat,Lx,Ly,Lz);

% Filter the high frequencies to avoid aliasing errors 

Uhatz=hfil_3d(Uhat,Lx,Ly,Lz);
Shatz=hfil_3d(Shat,Lx,Ly,Lz);
gradUhz=hfil_3d(gradUh,Lx,Ly,Lz);
gradShz=hfil_3d(gradSh,Lx,Ly,Lz);

UdgradShat=nludgradshat_3d(Uhatz,gradShz);

SgUthat=nlsguthat_3d(gradUhz,Shatz);

ksq=k1.^2+k2.^2+k3.^2;

g(:,:,:,1)=nu*ksq*dt/2;
g(:,:,:,2)=nu*ksq*dt/2;
g(:,:,:,3)=nu*ksq*dt/2;
g(:,:,:,4)=nu*ksq*dt/2;
g(:,:,:,5)=nu*ksq*dt/2;
g(:,:,:,6)=nu*ksq*dt/2;

oldRHS=newRHS;

newRHS=shat_update_3d(UdgradShat,SgUthat,gradUh,lam,Shat);

Shat=((ones(grid.Nx,grid.Ny,grid.Nz,6)-g)./(ones(grid.Nx,grid.Ny,grid.Nz,6)+g)).*Shat+...         % ABCN update lamth visc. nu
    (1./(ones(grid.Nx,grid.Ny,grid.Nz,6)+g)).*dt/2.*(3*newRHS-oldRHS);   % if g=0 get regular AB


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function udgradshat=nludgradshat_3d(Uhate,gradShe)

% nludgradsa stands for nonlinear (U)dot(gradSigma) (and a for de-aliased)
% this function computes the above from Fourier data
% by filtering the spectrum and ifft and multiplying
% and fft and truncating

% Preallocate
sizeU = size(Uhate);
sizeS = size(gradShe);
U = zeros(sizeU(1), sizeU(2), sizeU(3), 3);
gS = zeros(sizeS(1), sizeS(2), sizeS(3), 18);
udgradsext = zeros(sizeU(1), sizeU(2), sizeU(3), 6);
udgradshat = zeros(sizeU(1), sizeU(2), sizeU(3), 6);
for ind=1:3
    U(:,:,:,ind)=real(ifftn(Uhate(:,:,:,ind)));
end
for ind = 1:18
    gS(:,:,:,ind)=real(ifftn(gradShe(:,:,:,ind)));
end
udgradsext(:,:,:,1)=U(:,:,:,1).*gS(:,:,:, 1)+U(:,:,:,2).*gS(:,:,:, 2)+U(:,:,:,3).*gS(:,:,:, 3);
udgradsext(:,:,:,2)=U(:,:,:,1).*gS(:,:,:, 4)+U(:,:,:,2).*gS(:,:,:, 5)+U(:,:,:,3).*gS(:,:,:, 6);
udgradsext(:,:,:,3)=U(:,:,:,1).*gS(:,:,:, 7)+U(:,:,:,2).*gS(:,:,:, 8)+U(:,:,:,3).*gS(:,:,:, 9);
udgradsext(:,:,:,4)=U(:,:,:,1).*gS(:,:,:,10)+U(:,:,:,2).*gS(:,:,:,11)+U(:,:,:,3).*gS(:,:,:,12);
udgradsext(:,:,:,5)=U(:,:,:,1).*gS(:,:,:,13)+U(:,:,:,2).*gS(:,:,:,14)+U(:,:,:,3).*gS(:,:,:,15);
udgradsext(:,:,:,6)=U(:,:,:,1).*gS(:,:,:,16)+U(:,:,:,2).*gS(:,:,:,17)+U(:,:,:,3).*gS(:,:,:,18);
for ind = 1:6
    udgradshat(:,:,:,ind)=fftn(udgradsext(:,:,:,ind)); 
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sigout=shat_update_3d(UdgradShat,SgUthat,gradUh,lam,Shat)


% This is the function associated lamth the AB time step
% it is in fourier space
% it looks like F=-U(dotgrad)S+gradU*S+S*(gradU)^T-(1/lam)*S+(1/lam)*2D;
% where 2D=gradU+gradU^T

% in fourier space

szf=size(Shat);
nx=szf(1);
ny=szf(2);
nz=szf(3);

%Preallocate
sigout = zeros(nx, ny, nz, 6);

sigout(:,:,:,1) =    -UdgradShat(:,:,:,1)+   2*(SgUthat(:,:,:,1)) ...
                    -(1/lam)*Shat(:,:,:,1)+(1/lam)*2*gradUh(:,:,:,1);

sigout(:,:,:,2) =    -UdgradShat(:,:,:,2)+   (SgUthat(:,:,:,2)+SgUthat(:,:,:,4))...
                    -(1/lam)*Shat(:,:,:,2) +(1/lam)*(gradUh(:,:,:,2)+gradUh(:,:,:,4));

sigout(:,:,:,3) =    -UdgradShat(:,:,:,3)+   (SgUthat(:,:,:,3)+SgUthat(:,:,:,7))...
                    -(1/lam)*Shat(:,:,:,3)+(1/lam)*(gradUh(:,:,:,3)+gradUh(:,:,:,7));

sigout(:,:,:,4) =    -UdgradShat(:,:,:,4)+   2*(SgUthat(:,:,:,5))...
                    -(1/lam)*Shat(:,:,:,4)+(1/lam)*2*gradUh(:,:,:,5);

sigout(:,:,:,5) =    -UdgradShat(:,:,:,5)+   (SgUthat(:,:,:,6)+SgUthat(:,:,:,8))...
                    -(1/lam)*Shat(:,:,:,5)+(1/lam)*(gradUh(:,:,:,6)+gradUh(:,:,:,8));

sigout(:,:,:,6) =    -UdgradShat(:,:,:,6)+   2*(SgUthat(:,:,:,9))  ...
                    -(1/lam)*Shat(:,:,:,6)+(1/lam)*2*gradUh(:,:,:,9);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sguthat=nlsguthat_3d(gradUhe,Shate)

% nlsguthat_3d stands for nonlinear (sigma)(gradU)^T (and a for de-aliased)
% this function computes the above from Fourier data
% by extending the spectrum and ifft and multiplying
% and fft and truncating and ifft again

% Preallocate
gU = zeros(size(gradUhe));
S = zeros(size(Shate));
for ind=1:9
    gU(:,:,:,ind)=real(ifftn(gradUhe(:,:,:,ind)));
end

for ind=1:6
    S(:,:,:,ind)=real(ifftn(Shate(:,:,:,ind)));
end

sizeS = size(S);
sgutext = zeros(sizeS(1), sizeS(2), sizeS(3), 9);
sguthat = zeros(size(sgutext));

sgutext(:,:,:,1)=S(:,:,:,1).*gU(:,:,:,1)+S(:,:,:,2).*gU(:,:,:,2)+S(:,:,:,3).*gU(:,:,:,3);
sgutext(:,:,:,2)=S(:,:,:,1).*gU(:,:,:,4)+S(:,:,:,2).*gU(:,:,:,5)+S(:,:,:,3).*gU(:,:,:,6);
sgutext(:,:,:,3)=S(:,:,:,1).*gU(:,:,:,7)+S(:,:,:,2).*gU(:,:,:,8)+S(:,:,:,3).*gU(:,:,:,9);
sgutext(:,:,:,4)=S(:,:,:,2).*gU(:,:,:,1)+S(:,:,:,4).*gU(:,:,:,2)+S(:,:,:,5).*gU(:,:,:,3);
sgutext(:,:,:,5)=S(:,:,:,2).*gU(:,:,:,4)+S(:,:,:,4).*gU(:,:,:,5)+S(:,:,:,5).*gU(:,:,:,6);
sgutext(:,:,:,6)=S(:,:,:,2).*gU(:,:,:,7)+S(:,:,:,4).*gU(:,:,:,8)+S(:,:,:,5).*gU(:,:,:,9);
sgutext(:,:,:,7)=S(:,:,:,3).*gU(:,:,:,1)+S(:,:,:,5).*gU(:,:,:,2)+S(:,:,:,6).*gU(:,:,:,3);
sgutext(:,:,:,8)=S(:,:,:,3).*gU(:,:,:,4)+S(:,:,:,5).*gU(:,:,:,5)+S(:,:,:,6).*gU(:,:,:,6);
sgutext(:,:,:,9)=S(:,:,:,3).*gU(:,:,:,7)+S(:,:,:,5).*gU(:,:,:,8)+S(:,:,:,6).*gU(:,:,:,9);

for ind = 1:9
    sguthat(:,:,:,ind)=fftn(sgutext(:,:,:,ind));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






