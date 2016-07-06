function [veforcehat] = get_veforcehat(Shat,xi,grid)


gradShat = matrix_derivative_fourier(Shat,grid.Lx,grid.Ly);

divShat(:,:,1)=gradShat(:,:,1)+gradShat(:,:,4);
divShat(:,:,2)=gradShat(:,:,3)+gradShat(:,:,6);

veforcehat = xi*divShat;