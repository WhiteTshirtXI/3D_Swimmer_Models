function  [Shat, newRHS] = update_Shat(Uhat,grid,Shat,nu,dt,lam,newRHS)

% update the Stress advection equation

% Note to translate into something that reads better use the following:
% 
% 
% S11 = S(:,:,1);
% S12 = S(:,:,2);
% S22 = S(:,:,3);
% 
% dxu = gU(:,:,1);
% dyu = gU(:,:,2);
% dxv = gU(:,:,3);
% dyv = gU(:,:,4);

% record the number of grid points in each direction
%
szf=size(Uhat);
nx=szf(1);
ny=szf(2);
Lx=grid.Lx;
Ly=grid.Ly;
% compute the wave numbers
%
N1x =  floor((nx-1)/2);
N2x = (nx/2)*ones(rem(nx+1,2));
freqx =(2*pi/Lx)* [(0:N1x)  N2x (-N1x:-1)]';


N1y =  floor((ny-1)/2);
N2y = (ny/2)*ones(rem(ny+1,2));
freqy = (2*pi/Ly)*[(0:N1y)  N2y (-N1y:-1)]';


[k1, k2]=ndgrid(freqx,freqy);

%update stress

gradUh=matrix_derivative_fourier(Uhat,grid.Lx,grid.Ly);

gradSh=matrix_derivative_fourier(Shat,grid.Lx,grid.Ly);

% Filter the high frequencies to avoid aliasing errors 

Uhatz=hfil(Uhat,grid.Lx,grid.Ly);
Shatz=hfil(Shat,grid.Lx,grid.Ly);
gradUhz=hfil(gradUh,grid.Lx,grid.Ly);
gradShz=hfil(gradSh,grid.Lx,grid.Ly);

UdgradShat=nludgradshat(Uhatz,gradShz);

SgUthat=nlsguthat(gradUhz,Shatz);

ksq=k1.^2+k2.^2;

g(:,:,1)=nu*ksq*dt/2;
g(:,:,2)=nu*ksq*dt/2;
g(:,:,3)=nu*ksq*dt/2;

oldRHS=newRHS;

newRHS=shat_update(UdgradShat,SgUthat,gradUh,lam,Shat);

Shat=((ones(grid.Nx,grid.Ny,3)-g)./(ones(grid.Nx,grid.Ny,3)+g)).*Shat+...         % ABCN update lamth visc. nu
    (1./(ones(grid.Nx,grid.Ny,3)+g)).*dt/2.*(3*newRHS-oldRHS);   % if g=0 get regular AB


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function udgradshat=nludgradshat(Uhatz,gradShz)

% nludgradsa stands for nonlinear (U)dot(gradSigma) (and a for de-aliased)
% this function computes the above from Fourier data
% by filtering the spectrum and ifft and multiplying
% and fft and truncating

U=real(ifft2(Uhatz));
gS=real(ifft2(gradShz));

udgrads(:,:,1)=U(:,:,1).*gS(:,:,1)+U(:,:,2).*gS(:,:,2);    % u dx S11 + v dy S11
udgrads(:,:,2)=U(:,:,1).*gS(:,:,3)+U(:,:,2).*gS(:,:,4);    % u dx S12 + v dy S12
udgrads(:,:,3)=U(:,:,1).*gS(:,:,5)+U(:,:,2).*gS(:,:,6);    % u dx S22 + v dy S22

udgradshat=fft2(udgrads);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sig2=shat_update(UdgradShat,SgUthat,gradUh,lam,Shat)


% This is the function associated lamth the AB time step
% it is in fourier space
% it looks like F=-U(dotgrad)S+gradU*S+S*(gradU)^T-(1/lam)*S+(1/lam)*2D;
% where 2D=gradU+gradU^T

% in fourier space



sig2(:,:,1)=-UdgradShat(:,:,1)+2*SgUthat(:,:,1)-(1/lam)*Shat(:,:,1)+(1/lam)*2*gradUh(:,:,1);
sig2(:,:,2)=-UdgradShat(:,:,2)+(SgUthat(:,:,2)+SgUthat(:,:,3))-(1/lam)*Shat(:,:,2)+(1/lam)*(gradUh(:,:,2)+gradUh(:,:,3));
sig2(:,:,3)=-UdgradShat(:,:,3)+2*SgUthat(:,:,4)-(1/lam)*Shat(:,:,3)+(1/lam)*2*gradUh(:,:,4);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sguthat=nlsguthat(gradUhz,Shatz)

% nlsguta stands for nonlinear (sigma)(gradU)^T (and a for de-aliased)
% this function computes the above from Fourier data
% by extending the spectrum and ifft and multiplying
% and fft and truncating and ifft again

gU=real(ifft2(gradUhz));
S=real(ifft2(Shatz));

sgut(:,:,1)=S(:,:,1).*gU(:,:,1)+S(:,:,2).*gU(:,:,2);
sgut(:,:,2)=S(:,:,1).*gU(:,:,3)+S(:,:,2).*gU(:,:,4);
sgut(:,:,3)=S(:,:,2).*gU(:,:,1)+S(:,:,3).*gU(:,:,2);
sgut(:,:,4)=S(:,:,2).*gU(:,:,3)+S(:,:,3).*gU(:,:,4);

sguthat=fft2(sgut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






