function [A] = kappa_burrower(S,t,Tper,L)
% Curvature function that defines the gait of the worm. 
% L - length of the worm
% Tper - gait time in seconds
% S - Coordinate matrix of worm
% t - time
% Burrower
en = size(S,1);
for i=1:en
    A(i) = (5.3 - 3.1*S(i))*cos(4*pi*(t - 0.25*S(i) + L/4));
end
A = A';
end
