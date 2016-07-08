%
% animate the swimmer in lab frame 
%


% this info should really be read in
%
Lx = 2;
Ly = 2;
xmin=-Lx/2;
ymin=-Ly/2;
K  = Lx/Ly;
Ny = 64;
Nx = K*Ny;
dx = Lx/Nx;

% time stepping info
%
t0    = 0.1;
dtout = 0.1;
Tend  = 2.0;


% grid point positions
%
x = xmin + dx*(0:Nx-1)';
y = ymin + dx*(0:Ny-1)';
[x,y]=ndgrid(x,y);


% record the number of outputs of position
%
k=1;
for t = t0:dtout:Tend
  filename = sprintf('./diffdata/fluidVelocity_t%f.mat',t);
  load(filename);
  quiver(x,y,fluidForce(:,:,1),fluidForce(:,:,2));
  hold on;
  axis([xmin xmin+Lx ymin ymin+Ly]);
  set(gca,'plotboxaspectratio',[Lx Ly 1]);
  pause(0.05);
  hold off;
  k = k+10;
end



