function [veforcehat] = get_veforcehat_3d(Shat,xi,grid)


gradShat = matrix_derivative_fourier_3d(Shat,grid.Lx,grid.Ly,grid.Lz);

% dxS11 = gs(1)   
% dyS11 = gs(2)
% dzS11 = gs(3)
% dxS12 = gs(4)
% dyS12 = gs(5)
% dzS12 = gs(6)
% dxS13 = gs(7)
% dyS13 = gs(8)
% dzS13 = gs(9)
% dxS22 = gs(10)
% dyS22 = gs(11)
% dzS22 = gs(12)
% dxS23 = gs(13)
% dyS23 = gs(14)
% dzS23 = gs(15)
% dxS33 = gs(16)
% dyS33 = gs(17)
% dzS33 = gs(18)

% divS_1 = dxS11+dyS12+dzS13 = 1 + 5 + 9
% divS_2 = dxS12+dyS22+dzS23 = 4 + 11 + 15
% divS_3 = dxS13+dyS23+dzS33 = 7 + 14 + 18

% Preallocate for speed.
sizeS = size(gradShat);
divShat = zeros(sizeS(1), sizeS(2), sizeS(3), 3);

divShat(:,:,:,1)=gradShat(:,:,:,1)+gradShat(:,:,:,5)+gradShat(:,:,:,9);
divShat(:,:,:,2)=gradShat(:,:,:,4)+gradShat(:,:,:,11)+gradShat(:,:,:,15);
divShat(:,:,:,3)=gradShat(:,:,:,7)+gradShat(:,:,:,14)+gradShat(:,:,:,18);

veforcehat = xi*divShat;