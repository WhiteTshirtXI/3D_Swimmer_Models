%
% animate the swimmer in lab frame 
%
load('./data/imworm_3D_R_2D_VE_t3.000000.mat');


% this info should really be read in
%
Lx = 2;
Ly = 2;
Lz = 2;
xmin=-Lx/2;
ymin=-Ly/2;
zmin=-Lz/2;
Nx = 64;
Ny = 64;
Nz = 64;
dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nx;

% record the number of outputs of position
%
Nt = size(XTworm,4);
for k=1:Nt
  plot3(XTworm(:,1,k),XTworm(:,2,k), XTworm(:,3,k),'bo');
  axis([xmin xmin+Lx ymin ymin+Ly zmin zmin+Lz]);
  set(gca,'plotboxaspectratio',[Lx Ly Lz]);
  pause(0.01);
end



