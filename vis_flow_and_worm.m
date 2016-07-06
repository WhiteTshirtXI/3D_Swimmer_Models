%
% animate the swimmer in lab frame 
%
load('./data/imworm_n032_t4.000000.mat');

% this info should really be read in
%
Lx = 2;
Ly = 1;
xmin=-Lx/2;
ymin=-Ly/2;
K  = Lx/Ly;
Ny = 32;
Nx = K*Ny;
dx = Lx/Nx;

% time stepping info
%
t0    = 0.1;
dtout = 0.1;
Tend  = 4.0;


% grid point positions
%
x = xmin + dx*(0:Nx-1)';
y = ymin + dx*(0:Ny-1)';
[x,y]=ndgrid(x,y);


% record the number of outputs of position
%
k=1;
for t = t0:dtout:Tend
  filename = sprintf('./data/exworm_n032_t%f.mat',t);
  load(filename);
  quiver(x,y,U(:,:,1),U(:,:,2));
  hold on;
  plot(XTworm(:,1,k),XTworm(:,2,k),'bo');
  axis([xmin xmin+Lx ymin ymin+Ly]);
  set(gca,'plotboxaspectratio',[Lx Ly 1]);
  pause(0.05);
  hold off;
  k = k+10;
end



