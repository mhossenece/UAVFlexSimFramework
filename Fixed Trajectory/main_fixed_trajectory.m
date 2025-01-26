clc;
clear;
close all;
global start_time 
start_time =[] % remove the previous value of start_time
figure;

ax = axes;

hold(ax, 'on');
grid on

xlabel(ax, 'Longitude');
ylabel(ax, 'Latitude');
zlabel(ax, 'Altitude (m)');

% Load the Geofence from KML File
geofence = kml2struct('AERPAW_UAV_Geofence_Phase_1.kml');

% Extract the geofence boundaries
geofence_lat = geofence.Lat;
geofence_lon = geofence.Lon;
geofence_alt = zeros(size(geofence_lat)); % Assuming geofence is at ground level

% Plot the geofence in 3D
%geof = plot3(ax, geofence_lon, geofence_lat, geofence_alt, 'k-', 'LineWidth', 2);

% eNodeB Positions
%            LW1, LW2, LW3, LW4
lat_eNBs = [35.7275, 35.728056, 35.725, 35.733056];
lon_eNBs = [-78.695833, -78.700833, -78.691667, -78.698333];
alt_eNBs = [10, 10, 10, 10]; % altitudes for eNodeBs
%data_sizes = [1, 4, 3, 2]; % [2 3 4 1] data sizes in Mbytes UDP max: 8kBytes=8192
%data_sizes = [100, 400, 200, 300]; % [2 4 3 1]
data_sizes = [10, 40, 20, 30]; % [2 4 3 1]
%data_sizes = [1, 2, 3, 4]; % [4 3 2 1]
%data_sizes = [3, 2, 1, 4]; % [4 1 2 3]
%data_sizes = [1, 2, 3, 4]; % [4 3 2 1]
is_data_collection = [];


% if sorted_indices (2) ==3 && sorted_indices (3) == 1 %[2 3 4 1] or [4 3 2 1] 
%     num_waypoints = 5
% elseif sorted_indices (3) == 3 && sorted_indices (4) == 1 % [2 4 3 1] or [4 2 3 1]
    
%data_sizes = [1024, 3072, 2048, 8192]; % data sizes in bytes UDP max: 8kBytes=8192
%#File size: 1KB=1024B, 2kB=2048B, 3kB=3072B, 4kB=4096B, 5kB=5120B, 6kB=6144, 7kB=7168B, 8kB=8192B, Xxx9kB=9216B

% Sorting eNodeBs by data size in descending order
[sorted_data_sizes, sorted_indices] = sort(data_sizes, 'descend');
sorted_lat_eNBs = lat_eNBs(sorted_indices);
sorted_lon_eNBs = lon_eNBs(sorted_indices);
adjust = 0
%colors = ["m", "g", "b", "black"];


%colors = ["k", "k", "k", "k"];
colors = ["m", "g", "b", "r"];
for i = 1:4
        if i==1
        sym = "s"; %square
    elseif i==2
        sym = "^"; %triangle
    elseif i==3
        sym = "d"; %diamond
    elseif i==4
        sym = "p"; %cross
        end
    %plot3(ax, lon_eNBs(i), lat_eNBs(i), alt_eNBs(i), colors(i) + sym, 'MarkerSize', 10, 'MarkerFaceColor', colors(i));
    plot3(ax, lon_eNBs(i), lat_eNBs(i), alt_eNBs(i), colors(i) + sym, 'MarkerSize', 10,'MarkerFaceColor', colors(i));
end


% Plot eNodeBs and their coverage ranges in 3D
%for i = 1:4
%    plot3(ax, lon_eNBs(i), lat_eNBs(i), alt_eNBs(i), colors(i) + "^", 'MarkerSize', 10, 'MarkerFaceColor', colors(i));
%    text(ax, lon_eNBs(i), lat_eNBs(i), alt_eNBs(i), ['LW' num2str(i)], 'Color', colors(i), 'FontSize', 10, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
%end


% Drone Position
lat_drone = 35.7274823;
lon_drone = -78.6962747;
alt_drone = -0.017; % Initial altitude

% Initial Plot of Drone Position in 3D
dronePlot = plot3(ax, lon_drone, lat_drone, alt_drone, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'black');
%droneOrigin = plot3(ax, lon_drone, lat_drone, alt_drone, 'co', 'MarkerSize', 8, 'MarkerFaceColor', 'cyan');

% Function to check if a point is inside the geofence
is_in_geofence = @(lat, lon) inpolygon(lon, lat, geofence_lon, geofence_lat);

% Define the bounding box for the geofence
min_lat = min(geofence_lat);
max_lat = max(geofence_lat);
min_lon = min(geofence_lon);
max_lon = max(geofence_lon);

% Initialize variables for the best trajectory
best_trajectory = [];
min_distance = Inf;

% Function to calculate 2D distance between two lat-lon points (Haversine formula)
haversine = @(lat1, lon1, lat2, lon2) 2 * 6371000 * ...
    asin(sqrt(sin(deg2rad(lat2 - lat1) / 2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(deg2rad(lon2 - lon1) / 2)^2));


% MGEN Traffic Setup
% Note: This is a simplified version and may require the actual MGEN command syntax and execution based on your environment and setup.

%Transmission
% 0.0 ON  1 UDP SRC 5001 DST 127.0.0.1/5000 POISSON [1 7187]
% 10.0 OFF 1

%#Receiver reception
% 0.0 LISTEN UDP 5000-5001,6000,6003
% 2.0 IGNORE UDP 5000
% 3.0 IGNORE UDP 5001,6000,6003
% 10.0 OFF 1
%PERIODIC
%0.0 ON %s UDP SRC %d DST %s/%s POISSON [1 %d]\n
%0.0 ON %s UDP SRC %d DST %s/%s PERIODIC [1 %d]\n

% Generate MGEN UDP traffic
generate_mgen_UDP = @(flow_id,src_port, dst_ip,dst_port, data_size) ...    
    sprintf(['#UDP Transmission\n' ...
    'TXBUFFER 10000\n' ...
    '0.0 ON %s UDP SRC %d DST %s/%s POISSON [1 512:%d]\n' ...
    '10 OFF %d\n' ...
    '#UDP Reception\n' ...
    '0.0 LISTEN UDP 5000-5001,6000,6003\n' ...
    '15.0 IGNORE UDP 5000\n' ...
    '20.0 IGNORE UDP 5001,6000,6003\n' ...
    '10.0 OFF 1'...
    ], ...
    flow_id, src_port, dst_ip, dst_port, data_size, ...
    flow_id);


%fprintf('%.1f ON %d TCP DST %.6f/%d PERIODIC [1 %d] COUNT 1\n', ...
 %           time, i, longitude, 7000 + i, 1048576);


% Function to generate MGEN traffic
generate_mgen_TCP1 = @(flow_id,src_port,dst_ip, dst_port, data_size) ...    
    sprintf(['#UDP Transmission\n' ...
    '0.0 ON 1 TCP SRC 5000 DST 192.168.1.100/7000 PERIODIC [1 5242880] COUNT 1\n' ...
    '10.0 OFF 1\n' ...
    '#TCP Reception\n' ...
    '0.0 LISTEN TCP 5000\n' ...
    '25.0 OFF 1'] ...
    );
    %src_port, src_port, dst_port, data_size);

% Function to generate MGEN traffic
generate_mgen_TCP1 = @(flow_id, dst_ip, dst_port, data_size) ...       
    sprintf(['#TCP Transmission\n' ...
    'TXBUFFER 10000\n' ...
    '21.0 ON %d TCP DST %s/%d PERIODIC [-1 %d] COUNT 1 \n' ...
    '10 OFF %d\n' ...    
    '10.0 IGNORE TCP %d' ...
    '#TCP Reception\n' ...
    '20.0 LISTEN TCP %d\n' ...
    '25.0 IGNORE TCP %d\n' ...    
    %'10.0 OFF 1'...
    ], ...
    flow_id, dst_ip, dst_port, data_size,flow_id,dst_port,...
    dst_port,dst_port);
                     %(sorted_indices_array(indexLW), dst_ip,dst_port, data_rate, data_size_LW * 1048576);    
generate_mgen_TCP = @(flow_id, dst_ip, dst_port, data_rate, data_size) ...       
    sprintf(['#TCP Transmission\n' ...
    'TXBUFFER 10000\n' ...
    '21.0 ON %d TCP DST %s/%d PERIODIC [%f %d] COUNT 1 \n' ...
    '10 OFF %d\n' ...    
    '10.0 IGNORE TCP %d' ...
    '#TCP Reception\n' ...
    '20.0 LISTEN TCP %d\n' ...
    '25.0 IGNORE TCP %d\n' ...    
    %'10.0 OFF 1'...
    ], ...
    flow_id, dst_ip, dst_port, data_rate, data_size,flow_id,dst_port,...
    dst_port,dst_port);


%rate = 100; % packets per second

% % Iterate through each base station to create MGEN traffic files
% dst_ip = '127.0.0.1';
% %flow_id = '1';
% for k = 1:4
%    % data_size = data_sizes(k);
%    % data_size = data_size * 1048576; % 1MBytes = 1048576 Bytes
%     mgen_command = generate_mgen_TCP(sorted_indices(k), dst_ip,dst_port, data_sizes(sorted_indices(k)) * 1048576);    
%     %mgen_command = generate_mgen_UDP(flow_id,src_port, dst_ip,dst_port, data_size);    
%     %mgen_filename = sprintf('mgen_traffic_LW%d_udp.mgn', k);
%     mgen_filename = sprintf('mgen_traffic_LW%d_tcp.mgn', sorted_indices(k));
%     fid = fopen(mgen_filename, 'w');
%     fprintf(fid, '%s\n', mgen_command);
%     fclose(fid);
% end

dst_ip = '127.0.0.1';




altitude = 30; %UAV altitude
waypoint_index = 3; %[4 1 2 3]

if sorted_indices (2) == 3 && sorted_indices (4) == 1 %[2 3 4 1] or [4 3 2 1] 
    num_waypoints = 5;

elseif sorted_indices (2) == 1 && sorted_indices (4) == 3 %[2 1 4 3] or [4 1 2 3] 
    num_waypoints = 5;

elseif sorted_indices (3) == 3 && sorted_indices (4) == 1 % [2 4 3 1] or [4 2 3 1]
    num_waypoints = 6;

else
    num_waypoints = 7;
end
%num_waypoints = 0; % 5 for 4 LWs 

% 7 if LW3-->LW1
% 6 if LW1-->LW3 | LW2-->LW3 for all
% and Vice versa
savDis = [];
savTrajWayps = [];    

maxTrajIter = 1;%20
best_index = 0;

% Optimum trajectory
for traj = 1:maxTrajIter
    waypoints = zeros(num_waypoints + 2, 3); % Including initial and final positions
    
    % Initial and final position (same)
    waypoints(1, :) = [lat_drone, lon_drone, -0.017];
    waypoints(2, :) = [lat_drone, lon_drone, 30];
    waypoints(end, :) = [lat_drone, lon_drone, -0.017];
    
    % Check for equal data sizes and sort by indices if equal
    for i = 1:length(sorted_data_sizes) - 1
        if sorted_data_sizes(i) == sorted_data_sizes(i + 1)
            % Sort by index to maintain a predictable order
            tie_indices = find(sorted_data_sizes == sorted_data_sizes(i));
            tie_sorted_indices = sort(sorted_indices(tie_indices));
            sorted_lat_eNBs(tie_indices) = lat_eNBs(tie_sorted_indices);
            sorted_lon_eNBs(tie_indices) = lon_eNBs(tie_sorted_indices);
        end
    end

    %Waypoints creation
    [waypoints, is_data_collection_ok] = create_waypoints(waypoints,sorted_indices,sorted_lat_eNBs,sorted_lon_eNBs,num_waypoints,geofence_lon,geofence_lat,min_lon,max_lon,min_lat,max_lat,lon_drone,lat_drone,altitude,waypoint_index);

    % Calculate the total distance of the trajectory
    total_distance = 0;
    for i = 1:size(waypoints, 1) - 1
        total_distance = total_distance + haversine(waypoints(i,1), waypoints(i,2), waypoints(i+1,1), waypoints(i+1,2));
    end
    
    savDis = [savDis; total_distance];
    savTrajWayps = [savTrajWayps; waypoints];
    
    % Update the best trajectory if this one is shorter
    if total_distance < min_distance
        min_distance = total_distance;
        best_trajectory = waypoints;        
        is_data_collection = is_data_collection_ok;
    end
end



% %Ensure again from beginning all waypoints are within the geofence 
for i = 1:size(waypoints, 1)
    if ~is_in_geofence(waypoints(i, 1), waypoints(i, 2))        
        error('Waypoint (%.6f, %.6f) is outside the geofence!', waypoints(i, 1), waypoints(i, 2));
    end
end


%% Write the best waypoints to a CSV file
filename = 'waypoints_best_trajectory.csv';

header = 'latitude,longitude,altitude';
fid = fopen(filename, 'w');

if fid == -1
    error('Cannot open file: %s', filename);
end

% Write the header to the file
fprintf(fid, '%s\n', header);

% Write each row to the file
for i = 1:size(best_trajectory, 1)
    fprintf(fid, '%.6f,%.6f,%.6f\n', best_trajectory(i, :));
end

fclose(fid);


%% Write all the trajectory waypoints in csv file
all_waypoints_in_csv(savTrajWayps, maxTrajIter, num_waypoints)

%% Keep distances in csv. At the end of csv file, you can find the min distance with index/ iteration number
% Write the waypoints to a CSV file
filename = 'distance_all.csv';
fid = fopen(filename, 'w');

if fid == -1
    error('Cannot open file: %s', filename);
end

for i = 1:traj % Write each row to the file
  %  i
    fprintf(fid, '%d\n', savDis(i));
end

fprintf(fid, 'Minimum Distance = %d at iteration %d', min(savDis(1:traj)),find(savDis(1:traj) == min_distance));
fclose(fid);


%% dlmwrite or csvwrite does not work in writing the full waypoint fraction
% part
% % Append the waypoints data to the file
% dlmwrite(filename, round(waypoints(:, 1:3), 4), '-append');
% fclose(fid)
% % Append the waypoints data to the file with specific formatting
% % for i = 1:size(waypoints, 1)
% %     fprintf(fid, '%.3f,%.3f,%\.4f\n', waypoints(i, 1), waypoints(i, 2), waypoints(i, 3));
% end
%csvwrite(filename, waypoints);
disp(['Waypoints have been written to ', filename]);

best_trajectory1 = [
    35.727482,-78.696275,0.000000
    35.727482,-78.696275,30.000000
    35.727973,-78.698737,30.000000
    35.729217,-78.698387,30.000000
    35.724449,-78.695997,30.000000
    35.726801,-78.696770,30.000000
    35.727550,-78.696224,30.000000
    35.727482,-78.696275,0.000000
];

best_trajectory = [
    35.7274823	-78.6962747	0
    35.7274823	-78.6962747	30
    35.728582	-78.699871	30
    35.730273	-78.698383	30
    35.724332	-78.69444	30
    35.727177	-78.696719	30
    35.727323	-78.696263	30
    35.7274823	-78.6962747	30
    35.7274823	-78.6962747	0
];

% Plot the path in 3D
%taking best waypoints
waypoints = best_trajectory;
geof = plot3(ax, geofence_lon, geofence_lat, geofence_alt, 'r-', 'LineWidth', 2);
trajectory = plot3(ax, waypoints(:, 2), waypoints(:, 1), waypoints(:, 3), 'k-.', 'LineWidth', 1.3);
%trajectory = plot3(ax, best_trajectory(:, 2), best_trajectory(:, 1), best_trajectory(:, 3), 'k-.', 'LineWidth', 1);

%legend([geof, trajectory, droneOrigin, dronePlot], {'Geofence', 'UAV Trajectory', 'UAV Initial Pos.', 'UAV Current Pos.'}, 'Location', 'northeast');
legend_entries = cell(1, 4); % Allocate 5 entries: 4 for eNodeBs + 1 for geofence


colors = ["m", "g", "b", "k"];

% Plot eNodeBs and their coverage ranges in 3D
for i = 1:4
%    plot3(ax, lon_eNBs(i), lat_eNBs(i), alt_eNBs(i), colors(i) + sym, 'MarkerSize', 10, 'MarkerFaceColor', colors(i));
    legend_entries{i} = ['LW' num2str(i)]; % Store the legend entry for this plot
    %text(ax, lon_eNBs(i), lat_eNBs(i), alt_eNBs(i), ['LW' num2str(i)], 'Color', colors(i), 'FontSize', 10, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end
legend_entries{5} = 'Drone'; 
legend_entries{6} = 'Geofence'; % First entry is for the geofence
legend_entries{7} = 'Trajectory';
legend(ax, legend_entries);

% Apply axis limits
xlim([-78.701 -78.6914]);
ylim([35.722 35.7339]);
xtickangle(45);


% Initialize arrays for SNR and RSRP measurements
num_waypoints = size(waypoints, 1);
SNR_values = zeros(num_waypoints, 4); % For four LWs
RSRP_values = zeros(num_waypoints, 4);
altitude_values = zeros(num_waypoints, 1); % To store altitude values

% Open files for writing
logFileIDs = cell(4, 1); %FOR rsrp, snr -- theoretical
logFileIDs_bias = cell(4, 1); %FOR rsrp, snr -- theoretical
for k = 1:4
    logFileIDs{k} = fopen(['LW' num2str(k) '_log.txt'], 'w');
    fprintf(logFileIDs{k}, 'time, Longitude, Latitude, Altitude, rsrp, snr, datarate\n');

    logFileIDs_bias{k} = fopen(['LW' num2str(k) '_log_bias.txt'], 'w');
    fprintf(logFileIDs_bias{k}, 'time, Longitude, Latitude, Altitude, rsrp\n');
end

%
% Open files for writing
% logTrafficIDs = cell(4, 1);
% for k = 1:4
%     logTrafficIDs{k} = fopen(['LW' num2str(k) '_logTraffic.txt'], 'w');
%     fprintf(logTrafficIDs{k}, 'time, Longitude, Latitude, Altitude, rsrp, snr, datarate\n');
% end


% AWGN Channel Simulation and SNR, RSRP Calculation
% Simulation parameters
transmit_power = 10; % dBm
noise_power = -90; % dBm

% Antenna gains 
antenna_gain_uav = 2; % dBi  omnidirectional antenna
antenna_gain_enb = 10; % dBi directional antenna


% Initialize start time
flight_start_time = tic;
pause_time_const = 1;

initial_flight_speed = 5; %m/s
bs_indx = 1;

%gps_data = cell(num_waypoints,4,1)

gps_data = cell(4);
total_time_TCP_starts = 0;
extra_waypoints_time = 0;
% Move the drone along the path
total_path_distance = 0;
for i=1:length(best_trajectory)-1
    total_path_distance = total_path_distance + haversine(best_trajectory(i,1), best_trajectory(i,2), best_trajectory(i+1,1), best_trajectory(i+1,2));
end

needed_time = total_path_distance / initial_flight_speed;
total_num_steps = ceil(needed_time / pause_time_const); 
signal_strength_for_handover = zeros(4,total_num_steps);
%traj_lat_lon_alt = zeros(1:3);

indx_step = 1;


for i = 1:size(waypoints, 1) - 1
    % Determine if current segment is ascent or descent
    
    if waypoints(i, 3) == 0 && waypoints(i + 1, 3) == 30
        %% Traffic
        %collect_files_from_bs(sorted_indices(i), logTrafficIDs{i}) 
        %sorten_indices(i) depending on data size, UAV will go to a BS

        initial_flight_speeds = .5;
        %initial_flight_speed = 10;
        % Ascending from 0m to 30m vertically first
        ascent_distance = abs(30 - 0);
        ascent_time = ascent_distance / initial_flight_speeds;
        ascending_time_vertically_30m = ascent_time;
        num_steps_ascent = ceil(ascent_time / pause_time_const)
       % pause_time_ascent = ascent_time / num_steps_ascent;
        pause_time_ascent = 1;

        lat_points_ascent = linspace(waypoints(i, 1), waypoints(i + 1, 1), num_steps_ascent);
        lon_points_ascent = linspace(waypoints(i, 2), waypoints(i + 1, 2), num_steps_ascent);
        alt_points_ascent = linspace(waypoints(i, 3), waypoints(i + 1, 3), num_steps_ascent);
        lat_points_ascent = linspace(waypoints(i, 1), waypoints(i + 1, 1), num_steps_ascent);

        % Display each value with 12 decimal places
        for j = 1:length(lat_points_ascent)
            fprintf('%.12f\n', lat_points_ascent(j));
        end

        
        for j = 1:num_steps_ascent
            set(dronePlot, 'XData', lon_points_ascent(j), 'YData', lat_points_ascent(j), 'ZData', alt_points_ascent(j));

            x = [lon_points_ascent(j), lat_points_ascent(j), alt_points_ascent(j)];   
            traj_lat_lon_alt(indx_step,1:3) = x;


            %Did not order sorted_lat_eNBs...
            for k = 1:4
                
                %alt = alt_points_ascent(j);
                %2D Distance
                %distance_to_enb = haversine(lat_points_ascent(j), lon_points_ascent(j), lat_eNBs(k), lon_eNBs(k));
                %3D Distance
                distance_to_enb = haversine3d(lat_points_ascent(j), lon_points_ascent(j),alt_points_ascent(j), lat_eNBs(k), lon_eNBs(k),alt_eNBs(k));
                %haversine3d(lat1, lon1, alt1, lat2, lon2, alt2)
                %path_loss = 20 * log10(distance_to_enb) + 40 + randn * 5;
                path_loss = 20 * log10(distance_to_enb) + 20 * log10(3410*10^6) - 147.55 ;
                RSRP = transmit_power + antenna_gain_enb - path_loss + antenna_gain_uav + adjust;                           
            
                %awgn_noise = randn * sqrt(0.5 * 10^(noise_power / 10));
                %constant_path_loss = 75.5 + randn;
                %RSRP2 = transmit_power + antenna_gain_enb - constant_path_loss + antenna_gain_uav + adjust ;                          
                %received_power = RSRP2 + awgn_noise; % 

                received_power = RSRP; %
                snrgp = snr_gap;
                SNR = received_power - noise_power ; %+ snrgp;

                RSRP_values(i, k) = RSRP;
                SNR_values(i, k) = SNR;
                
                signal_strength_for_handover(k, indx_step) = RSRP ; 

                %C = 0.5 * log 
                % Apply the Shannon-Hartley theorem to calculate channel capacity
                
                %C = 1.4 *log2(1 + 10^(SNR/10)); % Capacity in bits per second (bps)
                %ce = cod_eff();
                %C = C * ce; % 0.25;
                C2 = calculateDataRate(received_power, noise_power);
                C = C2; %(C2 + C)/2;

                timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
                %fprintf(logFileIDs{k}, '%s, %.6f, %.6f, %.2f, %.2f, %.2f\n', timestamp, lon_points_ascent(j), lat_points_ascent(j), alt_points_ascent(j), RSRP, SNR);
                fprintf(logFileIDs{k}, '%s, %.12f, %.12f, %.12f, %.12f, %.12f, %.12f\n', timestamp, lon_points_ascent(j), lat_points_ascent(j), alt_points_ascent(j), RSRP, SNR, C);
        

        %        fprintf('Vertially Ascent to 30m from 0m, Step %d, eNodeB %d: RSRP = %.2f dBm, SNR = %.2f dB,  C = %.2f Mbps\n', j, k, RSRP, SNR, C);

            end
            indx_step = indx_step + 1;
            pause(pause_time_ascent);
        end
        
        %collect_files_from_bs(sorted_indices(i), logTrafficIDs{i}) %collect data from LWs
    elseif waypoints(i, 3) == 30 && waypoints(i + 1, 3) == 0
        %collect_files_from_bs(sorted_indices(i), logTrafficIDs{i}) 
        initial_flight_speeds = 0.5
        % Descending from 30m to 0m vertically at the end
        descent_distance = abs(0 - 30)
        descent_time = descent_distance / initial_flight_speeds
        num_steps_descent = ceil(descent_time / pause_time_const)
        %pause_time_descent = descent_time / num_steps_descent;
        pause_time_descent = 1;

        lat_points_descent = linspace(waypoints(i, 1), waypoints(i + 1, 1), num_steps_descent);
        lon_points_descent = linspace(waypoints(i, 2), waypoints(i + 1, 2), num_steps_descent);
        alt_points_descent = linspace(waypoints(i, 3), waypoints(i + 1, 3), num_steps_descent);
        
        for j = 1:num_steps_descent
            set(dronePlot, 'XData', lon_points_descent(j), 'YData', lat_points_descent(j), 'ZData', alt_points_descent(j));
            x = [lon_points_descent(j), lat_points_descent(j), alt_points_descent(j)];   
            traj_lat_lon_alt(indx_step,1:3) = x;

            %Did not order sorted_lat_eNBs...
            for k = 1:4
                %2D Distance
                %distance_to_enb = haversine(lat_points_descent(j), lon_points_descent(j), lat_eNBs(k), lon_eNBs(k));
                %3D Distance
                distance_to_enb = haversine3d(lat_points_descent(j), lon_points_descent(j),alt_points_descent(j), lat_eNBs(k), lon_eNBs(k),alt_eNBs(k));
                
                %path_loss = 20 * log10(distance_to_enb) + 40 + randn * 5;
                path_loss = 20 * log10(distance_to_enb) + 20 * log10(3410*10^6) - 147.55; 
                RSRP = transmit_power + antenna_gain_enb - path_loss + antenna_gain_uav + adjust;
                %awgn_noise = randn * sqrt(0.5 * 10^(noise_power / 10));

                %constant_path_loss = 75.5 + randn;
                %RSRP2 = transmit_power + antenna_gain_enb - constant_path_loss + antenna_gain_uav + adjust ;                          
                %received_power = RSRP2 + awgn_noise; % 
                received_power = RSRP ; %+ awgn_noise; % 
                

                %received_power = RSRP + awgn_noise;
                %snrgp = snr_gap;
                SNR = received_power - noise_power ; %+ snrgp;

                RSRP_values(i, k) = RSRP;
                SNR_values(i, k) = SNR;
                %C = 1.4 *log2(1 + 10^(SNR/10));
                
                %C = C * cod_eff();
                
                C2 = calculateDataRate(received_power, noise_power);
                C = C2; %(C2 + C)/2;
                signal_strength_for_handover(k, indx_step) = RSRP;                     
                
                timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
                fprintf(logFileIDs{k}, '%s, %.12f, %.12f, %.12f, %.12f, %.12f, %.12f\n', timestamp, lon_points_descent(j), lat_points_descent(j), alt_points_descent(j), RSRP, SNR, C);
                                       
                %fprintf('Vertically Descent to 0m from 30m, Step %d, eNodeB %d: RSRP = %.2f dBm, SNR = %.2f dB, C = %.2f Mbps\n', j, k, RSRP, SNR, C);
            end
            indx_step = indx_step + 1;
            pause(pause_time_descent);
        end
        %collect_files_from_bs(sorted_indices(i-1), logTrafficIDs{i-1}) %collect data from LWs
    
    else
        %t_c(k) = t_c(k) + t_s(data_sizes(k)*8*1024*1024) / C + t_c_sening;  %data_colleciton_time
        t_c = zeros(1, 4);
        % Horizontal movement
        %distance = haversine(waypoints(i, 1), waypoints(i, 2), waypoints(i + 1, 1), waypoints(i + 1, 2));
        distance = haversine3d(waypoints(i, 1), waypoints(i, 2), waypoints(i, 3),waypoints(i + 1, 1), waypoints(i + 1, 2),waypoints(i+1, 3));
        time = distance / initial_flight_speed
        %total_time_TCP_starts = total_time_TCP_starts + time
        total_time_between_waypoints = time;
        num_steps = ceil(time / pause_time_const) 
        %pause_time = time / num_steps;
        pause_time = 1;

        lat_points = linspace(waypoints(i, 1), waypoints(i + 1, 1), num_steps);
        lon_points = linspace(waypoints(i, 2), waypoints(i + 1, 2), num_steps);
        alt_points = linspace(waypoints(i, 3), waypoints(i + 1, 3), num_steps);
        gps_data = [];
        data_rate = [];
        bias_val = 0
        for j = 1:num_steps
            set(dronePlot, 'XData', lon_points(j), 'YData', lat_points(j), 'ZData', alt_points(j));
            %haversine(waypoints(i, 1), waypoints(i, 2), waypoints(i + 1, 1), waypoints(i + 1, 2)); %each time measure distance when the UAV moves
            %Did not order sorted_lat_eNBs

            x = [lon_points(j), lat_points(j), alt_points(j)];   
            traj_lat_lon_alt(indx_step,1:3) = x;


            for k = 1:4
                %2D Distance
                %distance_to_enb = haversine(lat_points(j), lon_points(j), lat_eNBs(k), lon_eNBs(k));
                %3D Distance
                distance_to_enb = haversine3d(lat_points(j), lon_points(j),alt_points(j), lat_eNBs(k), lon_eNBs(k),alt_eNBs(k));
                
                path_loss = 20 * log10(distance_to_enb) + 20 * log10(3410*10^6) - 147.55 ; %+ randn * 5;
                RSRP = transmit_power + antenna_gain_enb - path_loss + antenna_gain_uav + adjust;
                
                %awgn_noise = randn * sqrt(0.5 * 10^(noise_power / 10));
                
                
                %constant_path_loss = 75.5 + randn;
                %RSRP2 = transmit_power + antenna_gain_enb - constant_path_loss + antenna_gain_uav + adjust ;                          
                %received_power = RSRP2 + awgn_noise; % 
                received_power = RSRP ; %+ awgn_noise; % 
                
                
                %received_power = RSRP + awgn_noise;


               % received_power = received_power - 30
                %snrgp = snr_gap;
                SNR = received_power - noise_power ;%+ snrgp;
                %C = 1.4 * log2(1 + 10^(SNR/10)); % Capacity in bits per second (bps) %data rate
                %C = C * cod_eff();
                C2 = calculateDataRate(received_power, noise_power);
                C = C2 ; %(C2 + C)/2
                %data_rate = [data_rate; C];
                data_rate(k) = C;
                RSRP_values(i, k) = RSRP;
                SNR_values(i, k) = SNR;
                
                signal_strength_for_handover(k, indx_step) = RSRP;

                rsrp_values_for_biasing (k) = RSRP;

                
% %%%%%Data Collection time
% 
%                 %%%%UAV energy model parameters
%                 delta_uav= 0.012;% a dimensionless efficiency factor related to the UAV's rotor efficiency
%                 p_uav=1.225;%kg/m^3 air density at sea level
%                 s_uav=0.05;% Solidity ratio, which is the ratio of the total rotor blade area to the disk area swept by the rotor
%                 A_uav=0.503;% pi*R^2 0.503 m^2
%                 omiga_uav=300; %rad/s Angular velocity of the rotor blades
%                 R_uav=0.4;%m Radius of the UAV rotor
%                 k_uav=0.1;% Correction factor for induced power, which accounts for additional power required due to various factors such as non-ideal flow conditions
%                 W_uav=20;% Weight of the UAV, which is the force due to gravity (mass times the acceleration due to gravity).
%                 Utip_uav=120;%m/s Tip speed of the rotor blade, calculated as Utip = omiga_uav * R
%                 v0_uav = 4.03;% m/s Induced velocity in hover condition %UAV s mean rotor induced velocity in hover
%                 d0_uav = 0.6; %an initial distance or displacement
%                 P0 = delta_uav/8*p_uav*s_uav*A_uav*omiga_uav^3*R_uav^3;
%                 Pi = (1+k_uav)*W_uav^(3/2)/(2*p_uav*A_uav)^(1/2);
%                 Ph = P0 + Pi;
%                 t_inter =1; p_th = 0.99; z_s = 0.01;
% 
%                 t_c_sensing = abs(t_inter * (log(1-p_th))/ log(1-exp(1)^-z_s*distance_to_enb))
% 
%                 if k==1
%                     t_c(k) = t_c(k)  + t_c_sensing;
%                 elseif k==2
%                     t_c(k) 
%                     t_c(k) = t_c(k) + t_c_sensing;
%                 elseif k==3
%                     t_c(k) = t_c(k) + t_c_sensing;
%                 elseif k==4
%                     t_c(k) = t_c(k) + t_c_sensing;
%                 end
                

                timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
                %fprintf(logFileIDs{k}, '%s, %.6f, %.6f, %.2f, %.2f, %.2f\n', timestamp, lon_points(j), lat_points(j), alt_points(j), RSRP_smooth, SNR_smooth);
                
                % Log GPS data
                gps_data = [gps_data ; sprintf('%.6f,%.6f,%.2f', lat_points(j), lon_points(j), alt_points(j))];
                %fprintf(logFileIDs{k}, '%s, %.6f, %.6f, %.2f, %.2f, %.2f, %.2f, %s\n', timestamp, lon_points(j), lat_points(j), alt_points(j), RSRP, SNR, C, gps_data);
                fprintf(logFileIDs{k}, '%s, %.12f, %.12f, %.12f, %.12f, %.12f, %.12f\n', timestamp, lon_points(j), lat_points(j), alt_points(j), RSRP, SNR, C);
                %fprintf('Horizontal Waypoint %d, Step %d, eNodeB %d: RSRP = %.2f dBm, SNR = %.2f dB, C = %.2f Mbps\n',i-1, j, k, RSRP, SNR, C);
                % Determine the handover points based on the maximum signal strength
                
            end
            indx_step = indx_step + 1;

           %  [~, best_bs] = max(rsrp_values_for_biasing);
           % 
           %  intended_bs = sorted_indices(bs_indx);
           % 
           %  if best_bs == intended_bs
           %      %x=0000000000000
           %      if best_bs == 1
           %          thres = 15;
           %      else
           %          thres = 0;
           %      end
           %      rsrp_values_for_biasing(intended_bs) = rsrp_values_for_biasing(intended_bs) + bias_val + 0.25 + thres
           %  else %add bias
           %      %biassss = rsrp_values_for_handover(best_bs)
           %      %biasss2 = rsrp_values_for_handover(intended_bs)
           %      if bias_val == 0
           %          new_bias = abs(rsrp_values_for_biasing(best_bs) - rsrp_values_for_biasing(intended_bs))
           %      else
           %          new_bias = 0.25;
           %      end
           % 
           %      bias_val = bias_val + new_bias ;
           %    %  rsrp_values_for_handover
           %      rsrp_values_for_biasing(intended_bs) = rsrp_values_for_biasing(intended_bs) + bias_val 
           %  end
           % 
           %  [~, biasing_bs] = max(rsrp_values_for_biasing)            ;
           %  %gps_data
           % 
           %  for k_bias = 1:4
           %      fprintf(logFileIDs_bias{k_bias}, '%s, %.6f, %.6f, %.2f, %.2f\n', timestamp, lon_points(j), lat_points(j), alt_points(j),rsrp_values_for_biasing(k_bias));
           %  end
           % 
           pause(pause_time);
           % % file:///home/mhossen/PhD_research/UAV_Trajectory_Optimization/study/June/summarized%20code/v3/main_data_mule_opt_trajectory.m

        end

        % k1 = t_c(1)  
        % k2 = t_c(2)
        % k3 = t_c(3)
        % k4 = t_c(4)

        % %collect_files_from_bs(sorted_indices(i-1), logTrafficIDs{i-1}) %collect data from LWs
        % data_collection_time_start = tic;
        % if is_data_collection(i) == 1            
        %     dataRate = data_rate(sorted_indices(bs_indx))            
        %     ascending_time_vertically_30m
        %     extra_waypoints_times = extra_waypoints_time
        %     total_time_between_waypoints = extra_waypoints_time + total_time_between_waypoints
        %     %collect_files_from_bs(lat_points(num_steps), lon_points(num_steps), alt_points(num_steps), sorted_lat_eNBs(bs_indx), sorted_lon_eNBs(bs_indx), sorted_indices, sorted_indices(bs_indx), data_sizes(sorted_indices(bs_indx)), total_time_TCP_starts,dataRate) %collect data from LWs
        %     collect_files_from_bs(lat_points(num_steps), lon_points(num_steps), alt_points(num_steps), sorted_lat_eNBs(bs_indx), sorted_lon_eNBs(bs_indx), sorted_indices, sorted_indices(bs_indx), data_sizes(sorted_indices(bs_indx)), total_time_between_waypoints,dataRate, ascending_time_vertically_30m) %collect data from LWs
        % 
        %     bs_indx = bs_indx + 1;
        % elseif is_data_collection(i) == 0
        %     extra_waypoints_time = total_time_between_waypoints
        % end
        % data_collection_time = toc(data_collection_time_start)
        % total_time_TCP_starts = total_time_TCP_starts + data_collection_time
        % %lat,alt,long when UAV collect data -> num_steps
        % 

        for i=1:length(lat_points)
           % i
            if ~is_in_geofence(lat_points(i), lon_points(i))        
                error('Waypoint (%.6f, %.6f) is outside the geofence!', lat_points(i), lon_points(i));
            end
        end

    end
end


set(gcf, 'PaperUnits', 'inches'); % Set paper units to inches
set(gcf, 'PaperSize', [6.5 4.5]); % Set to standard letter size (width, height)
set(gcf, 'PaperPosition', [0 0 8.5 11]); % Ensure the figure fits within the page

% Save the figure as a PDF
print(gcf, 'fixed.pdf', '-dpdf', '-r300', '-bestfit');
print(gcf, 'fixed.png', '-dpng', '-r300');
saveas(gcf, 'fixed.fig');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Handover %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off
% % Plot signal strengths and handovers
% [~, best_base_station] = max(signal_strength_for_handover);
% signal_strength = signal_strength_for_handover;
% best_bs =  best_base_station;
% figure;
% hold on;
% colors = ['r', 'g', 'b', 'm'];
% t = 1: length(signal_strength);
% 
% for i = 1:4
%     plot(t, signal_strength(i, :), 'Color', colors(i), 'DisplayName', ['LW ', num2str(i)],'LineWidth', 1.5);
% end
% plot(t, best_bs, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Handover');
% legend;
% xlabel('Time (s)');
% ylabel('RSRP (dBm)');
% title('Signal Strength from Base Stations and Handover Points');
% grid on;
% 
% % Display handover points
% disp('Handover Points (Time, Best BS):');
% disp([t; best_bs]);
% 
% UAV_trajectory = traj_lat_lon_alt
% % Overlay handover points on the UAV trajectory plot
% figure;
% %plot(UAV_trajectory(1, :), UAV_trajectory(2, :), 'b-', 'LineWidth', 2);
% x=UAV_trajectory(:,1)
% y=UAV_trajectory(:,2)
% plot(x, y, 'b-', 'LineWidth', 1.5)
% signal_strength = [x'; y']
% hold on;
% %plot(BS_positions(:, 1), BS_positions(:, 2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% scatter(signal_strength(1, best_bs == 1), signal_strength(2, best_bs == 1), 'r*');
% scatter(signal_strength(1, best_bs == 2), signal_strength(2, best_bs == 2), 'g*');
% scatter(signal_strength(1, best_bs == 3), signal_strength(2, best_bs == 3), 'c*');
% scatter(signal_strength(1, best_bs == 4), signal_strength(2, best_bs == 4), 'k*');
% legend('UAV Trajectory', 'LW1', 'Handover to LW2', 'Handover to LW3', 'Handover to LW4', 'Handover to BS 4');
% xlabel('Longitude');
% ylabel('Latitude');
% title('UAV Trajectory with Handovers to Base Stations');
% grid on;
% 

% Calculate total flight time
total_flight_time = toc(flight_start_time);
minutes = floor(total_flight_time / 60);
seconds = mod(total_flight_time, 60);
fprintf('Total Flight Time: %d minutes and %.2f seconds\n', minutes, seconds);



% Close log files
for k = 1:4
    fclose(logFileIDs{k});
%    fclose(logFileIDs_bias{k}); %traffic
%    fclose(logTrafficIDs{k}); %traffic
end


hold(ax, 'off');


% % Output waypoints for reference
% disp('Waypoints (Lat, Lon, Alt):');
% disp(waypoints);
% 
% snr_values = [];
% rsrp_values = [];
% for k=2:1:length(SNR_values)-2
% %for k = sorted_indices+1
%     rsrp_values = [rsrp_values; RSRP_values(k,:)];
%     snr_values  = [snr_values; SNR_values(k,:)];
% end
% 
% snr_values_waypoint_wise = snr_values(:,sorted_indices); %Waypoint_wise: 1,2,3,4...
% rsrp_values_waypoint_wise = rsrp_values(:,sorted_indices); %Waypoint_wise: 1,2,3,4...
% 
% %create legend
% csvwrite("legend_order.csv",generate_legend_order(sorted_indices));
% sorted_indices
% disp('RSRP values (dBm) at each waypoint:');
% csvwrite("rsrp_values.csv",rsrp_values_waypoint_wise);
% disp(rsrp_values_waypoint_wise);
% sorted_indices
% disp('SNR values (dB) at each waypoint:');
% csvwrite("snr_values.csv",snr_values_waypoint_wise);
% disp(snr_values_waypoint_wise);

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
    c = 2 * R * asin(sqrt(a));

    % 3D distance in meters
    distance = sqrt((c)^2 + dAlt^2);
end

function createMGNTrafficFile(sorted_indices_array, indexLW,src_ip, src_port, dst_ip, dst_port, total_time_between_waypoints, data_size_LW, theoretical_data_rate, ascending_time_vertically_30m)
    global start_time

    if isempty(start_time)
        current_time = 0;
        ascending_time_vertically_30m = floor((ascending_time_vertically_30m * 20 /5) * 3); % 3 times higher time needed
    else
        current_time = start_time;
        ascending_time_vertically_30m = 0
    end


    flow_id = indexLW
    dst_port = 5000 + flow_id

    %data_size_LW --> MBytes
    data_size = data_size_LW * 1024*1024*8  % in bits. 1MBytes = 1024 KBytes = 1024*1024 Bytes = 1024*1024*8 bits
    theoretical_data_ratesss = theoretical_data_rate
    data_collection_time = (data_size / theoretical_data_rate) / (1024*1024) 
    % % %%%%UAV energy model parameters
    % delta_uav= 0.012;% a dimensionless efficiency factor related to the UAV's rotor efficiency
    % p_uav=1.225;%kg/m^3 air density at sea level
    % s_uav=0.05;% Solidity ratio, which is the ratio of the total rotor blade area to the disk area swept by the rotor
    % A_uav=0.503;% pi*R^2 0.503 m^2
    % omiga_uav=300; %rad/s Angular velocity of the rotor blades
    % R_uav=0.4;%m Radius of the UAV rotor
    % k_uav=0.1;% Correction factor for induced power, which accounts for additional power required due to various factors such as non-ideal flow conditions
    % W_uav=20;% Weight of the UAV, which is the force due to gravity (mass times the acceleration due to gravity).
    % Utip_uav=120;%m/s Tip speed of the rotor blade, calculated as Utip = omiga_uav * R
    % v0_uav = 4.03;% m/s Induced velocity in hover condition %UAV s mean rotor induced velocity in hover
    % d0_uav = 0.6; %an initial distance or displacement
    % P0 = delta_uav/8*p_uav*s_uav*A_uav*omiga_uav^3*R_uav^3;
    % Pi = (1+k_uav)*W_uav^(3/2)/(2*p_uav*A_uav)^(1/2);
    % Ph = P0 + Pi;
    % t_inter =1; p_th = 0.99; z_s = 0.01;
    % 
    % 
    %distance = total_time_between_waypoints* 5 % 5m/s
    %t_c = abs(t_inter * (log(1-p_th))/ log(1-exp(1)^-z_s*distance))

    travelling_to_waypoint_time = ceil(total_time_between_waypoints * 20 / 5)  %5 m/s in AERPAW. I cannot change %extra 2s for adjustment
    %travelling_to_waypoint_time =  floor(total_time_between_waypoints * 20 / 5)  %5 m/s in AERPAW. I cannot change %extra 2s for adjustment
    %we cannot make it waiting. So, keep it fixed 3 times.

    tcp_flow_sender_start_time = ceil(ascending_time_vertically_30m  + current_time) % + travelling_to_waypoint_time
    tcp_flow_receiver_start_time = tcp_flow_sender_start_time - 10 % before 1sec, receiver will start listening to the port
    %adjusting the simulation and real environment 10 times the data_collection_time

    %tcp_flow_sender_off_time = ceil(tcp_flow_sender_start_time + data_collection_time * 11)
    tcp_flow_sender_off_time = ceil(tcp_flow_sender_start_time + data_collection_time * 9 + 3 * travelling_to_waypoint_time)
    tcp_flow_receiver_off_time = tcp_flow_sender_off_time

    hold_time = ceil(tcp_flow_receiver_off_time - tcp_flow_sender_start_time - travelling_to_waypoint_time)
    %start_time = hold_time + tcp_flow_sender_start_time + travelling_to_waypoint_time
    start_time = tcp_flow_sender_off_time + 1
    %PERIODIC
    mgen_command = sprintf(['' ...
    '#TCP Transmission\n' ...
    '%.1f ON %d TCP DST %s/%d PERIODIC [1 %d] COUNT 1\n' ...    
    '%.1f OFF %d\n' ...        
    ...
    '#TCP Reception\n' ...
    '%.1f LISTEN TCP %d\n' ...
    '%.1f IGNORE TCP %d\n' ...
    '%.1f OFF %d\n' ...        ' ...
    '#Hold_time_in_QGroundControld = %.1f\n' ...
    ], ...
    tcp_flow_sender_start_time,flow_id, dst_ip, dst_port, data_size, tcp_flow_sender_off_time, flow_id, ...
    tcp_flow_receiver_start_time, dst_port,tcp_flow_receiver_off_time, dst_port, tcp_flow_receiver_off_time,flow_id, hold_time);
    mgen_filename = sprintf('mgen_traffic_LW%d_tcp.mgn', indexLW);
    fid = fopen(mgen_filename, 'w');
    fprintf(fid, '%s\n', mgen_command);
    fclose(fid);
end
    
        


% Function to simulate UAV collecting files from base stations
%function collect_files_from_bs(base_station_index, data_length)
%collect_files_from_bs(lat_points(num_steps), lon_points(num_steps), alt_points(num_steps), lat_eNBs(i-1), lon_eNBs(i-1), sorted_indices(i-1), data_sizes(sorted_indices(i-1))) %collect data from LWs
function collect_files_from_bs(UAV_lat, UAV_lon, UAV_alt, eNB_lat, eNB_lon, sorted_indices_array, base_station_index, data_length, total_time_TCP_starts,theoretical_data_rate, ascending_time_vertically_30m) %collect data from LWs
    %UAV_lat, UAV_lon, UAV_alt, eNB_lat, eNB_lon
    sorted_indices_array
    src_ip = '127.0.0.1';
    %dst_ip = '192.168.1.10';
    dst_ip = '127.0.0.1';
    dst_port = 5000;
    src_port = 5001;
    createMGNTrafficFile(sorted_indices_array, base_station_index,src_ip, src_port, dst_ip,dst_port,total_time_TCP_starts, data_length,theoretical_data_rate, ascending_time_vertically_30m);
    fprintf('\nUAV is collecting data from LW%d (Data size: %d Mbytes)\n.........................................................\n',base_station_index,data_length);
%    base_station_index
%    dataRate
    %lw_config ="mgen_traffic_LW"+base_station_index+"_udp.mgn";
    lw_config ="mgen_traffic_LW"+base_station_index+"_tcp.mgn";
    %lw_config ="tcp_sample.mgn";
    mgenCommand ="mgen input "+lw_config+" output mgen_traffic_" + "LW_"+base_station_index+".drc"; %input saveLogLW1.mgn save saveLogLW1.mgn log logLW1.drc"
    %mgenCommand ="mgen input "+lw_config;
    % Execute the command using system
    status = system(mgenCommand);    
end



function datarate = calculateDataRate(received_power, noise_power)
    % Calculate SNR
    snr_gap = 0;
    snr = received_power - noise_power - snr_gap;
    
    % Get spectral efficiency based on SNR
    spectral_efficiency = getSpectralEfficiency(snr);
    
    % Calculate data rate
    datarate = spectral_efficiency * 1.4;  % Bandwidth in MHz
end

function dr_chnge = cod_eff()
    % Define the range
    a = 0.18; % Lower bound
    b = 0.25; % Upper bound
    
    % Generate an array of 10 random numbers in the range [0.1, 0.25]
    randomValues = a + (b - a) * rand;
    dr_chnge = randomValues;
    
    % Display the results
    %disp(randomValues);
end

function snr_gp = snr_gap()
      snr_gp = 0 ; %randn + rand ;%+ randn;
end

function se = getSpectralEfficiency(snr)
    %snr
    % Define spectral efficiency based on SNR ranges
    if snr >= 22.7
        se = 5.55;
    elseif snr >= 21.0
        se = 5.12;
    elseif snr >= 18.7
        se = 4.52;
    elseif snr >= 16.3
       se = 3.90;
    elseif snr >= 14.1
        se = 3.32;
    elseif snr >= 11.7
        se = 2.73;
    elseif snr >= 10.3
        se = 2.41;
    elseif snr >= 8.1
        se = 1.91;
    elseif snr >= 5.9
        se = 1.48;
    elseif snr >= 4.3
        se = 1.18;
    elseif snr >= 2.4
        se = 0.88;
    elseif snr >= 0.2
        se = 0.60;
    elseif snr >= -2.3
        se = 0.38;
    elseif snr >= -4.7
        se = 0.23;
    elseif snr >= -6.7
        se = 0.15;
    else
        se = 0.0;
    end
end
