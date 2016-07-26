function Fhat=forcing_4roll_3d(grid)

% This creates the fourier transform of the force field which is periodic
% in x  y z directions


Lx = grid.Lx;
Ly = grid.Ly;
Lz = grid.Lz;

xmin = grid.xmin;
ymin = grid.ymin;
zmin = grid.zmin;

dx = grid.dx;


[x,y,z]=ndgrid(xmin:dx:xmin+Lx-dx,ymin:dx:ymin+Ly-dx,zmin:dx:zmin+Lz-dx);

F(:,:,:,1) =  2*sin(2*pi/Lx*x).*cos(2*pi/Ly*y);
F(:,:,:,2) =   -2*cos(2*pi/Lx*x).*sin(2*pi/Ly*y);
F(:,:,:,3) =  zeros(grid.Nx,grid.Ny,grid.Nz);


for ind = 1:3
    Fhat(:,:,:,ind) = fftn(F(:,:,:,ind));
end



