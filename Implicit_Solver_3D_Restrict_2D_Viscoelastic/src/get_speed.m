function [MAV, U, Xcm]=get_speed(XTworm)

% function [MAV, U]=get_speed(XTworm)
%
% Note that this function assumes that your data has
% size(XTworm) = (Ns,2,T) where Ns is the number of points along the
% swimmer, 2 for the x,y positions of the swimmer, and T the time - vector
% it is assumed that you have saved 100 time slices per period of the
% stroke


Nt_per_T= 100;  % here is where 100 time slices is assumed
Xcm=squeeze(mean(XTworm)); % get the center of mass in x,y

XXX = Xcm(:,(Nt_per_T+1):end) - Xcm(:,1:(end-Nt_per_T));  % take moving average over period
XXX = sqrt( sum( XXX.^2 ) );
a=1; b=ones(Nt_per_T,1)/Nt_per_T;
MAV = filter(b,a,XXX); 
U = MAV(end);
