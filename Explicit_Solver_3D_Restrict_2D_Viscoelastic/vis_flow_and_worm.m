%
% animate the swimmer in lab frame 
%


% this info should really be read in
%
Lx = 2;
Ly = 2;
Lz = 3;
xmin=-Lx/2;
ymin=-Ly/2;
zmin=-Lz/2;
Ny = 64;
Nx = 64;
Nz = 64;
dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;

% time stepping info
%
t0    = 0.1;
dtout = 0.1;
Tend  = 10.0;


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
  filename = sprintf('./data/exworm_n128_t%f.mat',t);
  load(filename);
  quiver(x,y,z,U(:,:,:,1),U(:,:,:,2),U(:,:,:,3));
  hold on;
  plot(XTworm(:,1,k),XTworm(:,2,k),'bo');
  axis([xmin xmin+Lx ymin ymin+Ly zmin zmin+Lz]);
  set(gca,'plotboxaspectratio',[Lx Ly Lz 1]);
  pause(0.05);
  hold off;
  k = k+10;
end


