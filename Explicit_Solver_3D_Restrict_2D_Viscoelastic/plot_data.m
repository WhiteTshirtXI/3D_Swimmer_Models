clear
% Load first data set
load('./data/imworm_3D_VE_n128_t3.000000.mat');
    X1 = XTworm;
    U1 = U;
    Sh1 = Shat;
% Load second data set
load('./data/imworm_3D_VE_n128_t3.000000.mat');
    X2 = XTworm;
    U2 = U;
    Sh2 = Shat;
% Get worm statistics
[MAV1, F1, XCM1] = get_speed(X1);
[MAV2, F2, XCM2] = get_speed(X2);
% Plot Moving Average
    figure
    % Plot Blue Dots
    plot(MAV1, 'b.');
    hold
    % Plot Red O's
    plot(MAV2, 'ro');
    % Labels
    title('Worm Moving Average');
    xlabel('Time in .01 secs');
    ylabel('Speed in X grid units per T secs');
% Plot Center Of Mass
    % X-Coordinate
    figure
    % Plot Blue Dots
    plot(XCM1(1,:)', 'b.');
    hold
    %Plot Red O's
    plot(XCM2(1,:)', 'ro');
    % Labels
    title('Center of Mass');
    xlabel('Time in .01 secs');
    ylabel('X-coordinate');
    % Y-Coordinate
    figure
    % Plot Blue Dots
    plot(XCM1(2,:)', 'b.');
    hold
    %Plot Red O's
    plot(XCM2(2,:)', 'ro');
    % Labels
    title('Center of Mass');
    xlabel('Time in .01 secs');
    ylabel('Y-coordinate');
    % Z-Coordinate
    figure
    % Plot Blue Dots
    plot(XCM1(3,:)', 'b.');
    hold
    %Plot Red O's
    plot(XCM2(3,:)', 'ro');
    % Labels
    title('Center of Mass');
    xlabel('Time in .01 secs');
    ylabel('Z-coordinate');
% To-Do : Bending Energy L2-Norm

for i = 1:6
    S1(:,:,:,i) = real(ifftn(Sh1(:,:,:,i)));
end

trs = S1(:,:,:,1) + S1(:,:,:,4) + S1(:,:,:,6);
% Note: Adjust increment value to increase plot slice density
Sx = -1:0.2:1;
Sy = -1:0.2:1;
Sz = -1:0.2:1;
cvals = linspace(-1,1,128);
[X,Y,Z] = meshgrid(linspace(-1,1,128),linspace(-1,1,128),linspace(-1,1,128));

% First 3-D Figure (Worm Focus)
figure
contourslice(X,Y,Z,trs,Sx,Sy,Sz,cvals);
axis([-1,1,-1,1,-1,1]);
daspect([1,1,1]);
campos([10,-20,10]);
title('Worm View');
box on
% Second 3-D Figure (Full Focus)
figure
contourslice(X,Y,Z,trs,Sx,Sy,Sz);
axis([-1,1,-1,1,-1,1]);
daspect([1,1,1]);
campos([10,-20,10]);
title('Full View');
box on