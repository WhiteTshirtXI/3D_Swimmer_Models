%
% animate the swimmer in lab frame 
%


% this info should really be read in
%
Lx = 2;
Ly = 2;
Lz = 2;
xmin=-Lx/2;
ymin=-Ly/2;
zmin=-Lz/2;
Ny = 16;
Nx = 16;
Nz = 16;
dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;

% time stepping info
%
t0    = 0.1;
dtout = 0.1;
Tend  = 1.0;


% grid point positions
%
x = xmin + dx*(0:Nx-1)';
y = ymin + dy*(0:Ny-1)';
z = zmin + dz*(0:Nz-1)';
[x,y,z] = ndgrid(x,y,z);


% record the number of outputs of position
%
k=1;
for t = t0:dtout:Tend
  filename = sprintf('./data/imworm_3D_VE_n016_t%f.mat',t);
  load(filename);
  quiver3(x,y,z,U(:,:,:,1),U(:,:,:,2),U(:,:,:,3));
  hold on;
  plot3(XTworm(:,1,k),XTworm(:,2,k),XTworm(:,3,k),'bo');
  axis([xmin xmin+Lx ymin ymin+Ly zmin zmin+Lz]);
  set(gca,'plotboxaspectratio',[Lx Ly Lz]);
  pause(0.2);
  hold off;
  k = k+10;
end



