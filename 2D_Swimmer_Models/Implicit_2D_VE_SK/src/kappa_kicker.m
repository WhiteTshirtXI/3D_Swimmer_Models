function [A] = kappa_kicker(S,t,Tper,L)
% Curvature function that defines the gait of the worm. 
% L - length of the worm
% Tper - gait time in seconds
% S - Coordinate matrix of worm
% t - time
% Kicker
en = size(S,1);
for i=1:en
    A(i) = (5.3 - 3.1*(L-S(i)))*cos(4*pi*((Tper-t) - 0.25*(L-S(i)) + L/4));
end
A = A';
end
