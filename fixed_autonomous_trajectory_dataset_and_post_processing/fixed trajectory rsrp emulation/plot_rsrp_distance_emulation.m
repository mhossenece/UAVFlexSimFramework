clc;
clear;
close all;

% Load the CSV files into tables
dataLW1 = readtable('lw1_merged_icc.csv'); % for LW1
dataLW2 = readtable('lw2_merged_icc.csv'); % for LW2
dataLW3 = readtable('lw3_merged_icc.csv'); % for LW3
dataLW4 = readtable('lw4_merged_icc.csv'); % for LW4

% Extract relevant columns for LW1
timeLW1 = dataLW1.time; % Assuming the time is in string format
latitudeLW1 = dataLW1.Latitude;
longitudeLW1 = dataLW1.Longitude;
altitudeLW1 = dataLW1.Altitude; % Assuming altitude is provided in meters
rsrpLW1 = dataLW1.rsrp;
snrLW1 = dataLW1.snr;
%datarateLW1 = dataLW1.datarate;

% Extract relevant columns for LW2
timeLW2 = dataLW2.time; % Assuming the time is in string format
latitudeLW2 = dataLW2.Latitude;
longitudeLW2 = dataLW2.Longitude;
altitudeLW2 = dataLW2.Altitude; % Assuming altitude is provided in meters
rsrpLW2 = dataLW2.rsrp;
snrLW2 = dataLW2.snr;
%datarateLW2 = dataLW2.datarate;


% Extract relevant columns for LW3
timeLW3 = dataLW3.time; % Assuming the time is in string format
latitudeLW3 = dataLW3.Latitude;
longitudeLW3 = dataLW3.Longitude;
altitudeLW3 = dataLW3.Altitude; % Assuming altitude is provided in meters
rsrpLW3 = dataLW3.rsrp;
snrLW3 = dataLW3.snr;
%datarateLW3 = dataLW3.datarate;


% Extract relevant columns for LW4
timeLW4 = dataLW4.time; % Assuming the time is in string format
latitudeLW4 = dataLW4.Latitude;
longitudeLW4 = dataLW4.Longitude;
altitudeLW4 = dataLW4.Altitude; % Assuming altitude is provided in meters
rsrpLW4 = dataLW4.rsrp;
snrLW4 = dataLW4.snr;
%datarateLW4 = dataLW4.datarate;


% Calculate 2D distances between successive points in meters for LW1
distances_2d_meters_LW1 = zeros(size(latitudeLW1));
for i = 2:length(latitudeLW1)
    distances_2d_meters_LW1(i) = distances_2d_meters_LW1(i-1) + haversine(latitudeLW1(i-1), longitudeLW1(i-1), latitudeLW1(i), longitudeLW1(i));
end

% Calculate 2D distances between successive points in meters for LW2
distances_2d_meters_LW2 = zeros(size(latitudeLW2));
for i = 2:length(latitudeLW2)
    distances_2d_meters_LW2(i) = distances_2d_meters_LW2(i-1) + haversine(latitudeLW2(i-1), longitudeLW2(i-1), latitudeLW2(i), longitudeLW2(i));
end

% Calculate 2D distances between successive points in meters for LW3
distances_2d_meters_LW3 = zeros(size(latitudeLW3));
for i = 2:length(latitudeLW3)
    distances_2d_meters_LW3(i) = distances_2d_meters_LW3(i-1) + haversine(latitudeLW3(i-1), longitudeLW3(i-1), latitudeLW3(i), longitudeLW3(i));
end

% Calculate 2D distances between successive points in meters for LW4
distances_2d_meters_LW4 = zeros(size(latitudeLW4));
for i = 2:length(latitudeLW4)
    distances_2d_meters_LW4(i) = distances_2d_meters_LW4(i-1) + haversine(latitudeLW4(i-1), longitudeLW4(i-1), latitudeLW4(i), longitudeLW4(i));
end

% Calculate 3D distances between successive points in meters for LW1
distances_3d_meters_LW1 = zeros(size(latitudeLW1));
for i = 2:length(latitudeLW1)
    distances_3d_meters_LW1(i) = distances_3d_meters_LW1(i-1) + haversine3d(latitudeLW1(i-1), longitudeLW1(i-1), altitudeLW1(i-1), ...
        latitudeLW1(i), longitudeLW1(i), altitudeLW1(i));
end

% Calculate 3D distances between successive points in meters for LW2
distances_3d_meters_LW2 = zeros(size(latitudeLW2));
for i = 2:length(latitudeLW2)
    distances_3d_meters_LW2(i) = distances_3d_meters_LW2(i-1) + haversine3d(latitudeLW2(i-1), longitudeLW2(i-1), altitudeLW2(i-1), ...
        latitudeLW2(i), longitudeLW2(i), altitudeLW2(i));
end

% Calculate 3D distances between successive points in meters for LW3
distances_3d_meters_LW3 = zeros(size(latitudeLW3));
for i = 2:length(latitudeLW3)
    distances_3d_meters_LW3(i) = distances_3d_meters_LW3(i-1) + haversine3d(latitudeLW3(i-1), longitudeLW3(i-1), altitudeLW3(i-1), ...
        latitudeLW3(i), longitudeLW3(i), altitudeLW3(i));
end

% Calculate 3D distances between successive points in meters for LW4
distances_3d_meters_LW4 = zeros(size(latitudeLW4));
for i = 2:length(latitudeLW4)
    distances_3d_meters_LW4(i) = distances_3d_meters_LW4(i-1) + haversine3d(latitudeLW4(i-1), longitudeLW4(i-1), altitudeLW4(i-1), ...
        latitudeLW4(i), longitudeLW4(i), altitudeLW4(i));
end

figure
% Plot 2D Distance vs RSRP for LW1, LW2, LW3, and LW4
%subplot(2, 1, 2);
plot(distances_2d_meters_LW1, rsrpLW1, 'black-', 'LineWidth', 2); %k --> black
hold on;
%plot(distances_3d_meters_LW1, rsrpLW1, 'c-', 'LineWidth', 1.5);
plot(distances_2d_meters_LW2, rsrpLW2, 'r.-', 'LineWidth', 2);
%plot(distances_3d_meters_LW2, rsrpLW2, 'b-', 'LineWidth', 1.5);
plot(distances_2d_meters_LW3, rsrpLW3, 'g-.', 'LineWidth', 2);

%plot(distances_3d_meters_LW3, rsrpLW3, 'Color',[0.4940 0.1840 0.5560], 'LineWidth', 1.5);
plot(distances_2d_meters_LW4, rsrpLW4, 'b:', 'LineWidth',2);
%plot(distances_3d_meters_LW4, rsrpLW4, 'orange-', 'LineWidth', 1.5);
%plot(distances_3d_meters_LW4, rsrpLW4, 'Color', [1, 0.5, 0], 'LineWidth', 1.5); % Custom orange color
%title('Distance vs RSRP  measured in Matlab Simulation');
xlabel('Distance (m)','FontSize', 18);
ylabel('RSRP (dBm)','FontSize', 18);
set(gca, 'FontSize', 18); 
xlim([0 1780]);
%legend('2D Distance LW1', '3D Distance LW1', '2D Distance LW2', '3D Distance LW2', ...
%    '2D Distance LW3', '3D Distance LW3','2D Distance LW4', '3D Distance LW4');

%legend('2D Distance LW1', '3D Distance LW1', '2D Distance LW2', '3D Distance LW2', ...
%    '2D Distance LW3', '3D Distance LW3','2D Distance LW4', '3D Distance LW4');

legend('LW1', 'LW2', 'LW3', 'LW4', 'Position', [0.2, 0.7, 0.2, 0.2]); 

grid on;
hold off

set(gcf, 'PaperUnits', 'inches'); % Set paper units to inches
set(gcf, 'PaperSize', [6.5 4.5]); % Set to standard letter size (width, height)
set(gcf, 'PaperPosition', [0 0 8.5 11]); % Ensure the figure fits within the page

% Save the figure as a PDF
print(gcf, 'fixed_distance_rsrp_emulation.pdf', '-dpdf', '-r300', '-bestfit');
%print(gcf, 'fixed_distance_rsrp_emulation.png', '-dpng', '-r300');
%saveas(gcf, 'fixed_distance_rsrp_emulation.fig');






% Function to calculate 2D distance using Haversine formula
function distance = haversine(lat1, lon1, lat2, lon2)
    % Earth's radius in meters
    R = 6371000;  % Radius of Earth in meters

    % Convert degrees to radians
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    % Differences in latitudes and longitudes
    dLat = lat2 - lat1;
    dLon = lon2 - lon1;

    % Haversine formula for 2D distance
    a = sin(dLat/2)^2 + cos(lat1) * cos(lat2) * sin(dLon/2)^2;
    c = 2 * asin(sqrt(a));

    % Calculate the 2D distance in meters
    distance = R * c;
end

% Function to calculate 3D distance using Haversine formula
function distance = haversine3d(lat1, lon1, alt1, lat2, lon2, alt2)
    % Earth's radius in meters
    R = 6371000;  % Radius of Earth in meters

    % Convert degrees to radians
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    % Differences in latitudes, longitudes, and altitudes
    dLat = lat2 - lat1;
    dLon = lon2 - lon1;
    dAlt = alt2 - alt1;

    % Haversine formula for 3D distance
    a = sin(dLat/2)^2 + cos(lat1) * cos(lat2) * sin(dLon/2)^2;
    c = 2 * asin(sqrt(a));

    % Calculate the 3D distance in meters
    distance = sqrt((R * c)^2 + dAlt^2);
end
