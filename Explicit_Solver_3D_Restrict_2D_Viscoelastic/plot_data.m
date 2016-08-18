clear
% Load first data set
load('./data/imworm_3D_VE_n256_t3.000000.mat');
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
% Calculate norm and rescale
norm = sqrt(sum(S1(:).^2))*1e-2;
% Note: Adjust increment value to increase plot slice density
Sx = [];
Sy = [];
Sz = [0 0];
cvals = linspace(0,norm+1,2^11);
[X,Y,Z] = meshgrid(linspace(-1,1,256),linspace(-1,1,256),linspace(-1,1,256));

% First 3-D Figure (Worm Focus)
figure
contourslice(X,Y,Z,trs,Sx,Sy,Sz,cvals);
axis([-1,1,-1,1,-1,1]);
daspect([1,1,1]);
campos([10,-20,10]);
title('Stress Plot');
box on

for i = 1:6
    S2(:,:,:,i) = real(ifftn(Sh2(:,:,:,i)));
end

trs2 = S2(:,:,:,1) + S2(:,:,:,4) + S2(:,:,:,6);
% Calculate norm and rescale
norm2 = sqrt(sum(S2(:).^2))*1e-2;
% Note: Adjust increment value to increase plot slice density
Sx = [];
Sy = [];
Sz = [0 0];
cvals = linspace(0,norm2+1,2^11);
[X,Y,Z] = meshgrid(linspace(-1,1,128),linspace(-1,1,128),linspace(-1,1,128));

% First 3-D Figure (Worm Focus)
figure
contourslice(X,Y,Z,trs2,Sx,Sy,Sz,cvals);
axis([-1,1,-1,1,-1,1]);
daspect([1,1,1]);
campos([10,-20,10]);
title('Stress Plot');
box on