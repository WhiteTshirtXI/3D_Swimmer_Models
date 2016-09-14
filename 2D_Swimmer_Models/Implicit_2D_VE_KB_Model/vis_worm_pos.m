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


% record the number of outputs of position
%
Nt = size(XTworm,3);
for k=1:Nt
  plot(XTworm(:,1,k),XTworm(:,2,k),'bo');
  axis([xmin xmin+Lx ymin ymin+Ly]);
  set(gca,'plotboxaspectratio',[Lx Ly 1]);
  pause(0.01);
end



