clear
% Load first data set
load('./data/imworm_3D_VE_n128_t3.000000.mat');
    X1 = XTworm;
    U1 = U;
% Load second data set
load('./data/imworm_3D_VE_n128_t2.900000.mat');
    X2 = XTworm;
    U2 = U;
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