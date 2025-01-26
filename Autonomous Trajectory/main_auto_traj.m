clc;
clear all;
close all;

% Initialize plot
ax = axes;
hold(ax, 'on');
grid on;
xlabel(ax, 'Longitude');
ylabel(ax, 'Latitude');
zlabel(ax, 'Altitude (m)');

% Load the Geofence from KML File
geofence = kml2struct('AERPAW_UAV_Geofence_Phase_1.kml');

% Extract the geofence boundaries
geofence_lat = geofence.Lat;
geofence_lon = geofence.Lon;

% Define eNodeB locations
lat_eNBs = [35.7275, 35.728056, 35.725, 35.733056];
lon_eNBs = [-78.695833, -78.700833, -78.691667, -78.698333];
alt_eNBs = [10, 10, 10, 10]; % altitudes for eNodeBs

% Plot the eNodeBs
colors = ["m", "g", "b", "r"];
%colors = ["k", "k", "k", "k"];
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
    plot3(ax, lon_eNBs(i), lat_eNBs(i), alt_eNBs(i), colors(i) + sym, 'MarkerSize', 10, 'MarkerFaceColor', colors(i));
end

% Drone Initialization
global lat_drone lon_drone alt_drone
lat_drone = 35.7274823;
lon_drone = -78.6962747;
alt_drone = 0; % Start from ground (altitude = 0)
% Store the launch position for later return
launch_lat = lat_drone;
launch_lon = lon_drone;

% Initial drone position plot
%dronePlot = plot3(lon_drone, lat_drone, alt_drone, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
dronePlot = plot3(lon_drone, lat_drone, alt_drone, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

% Plot the geofence
plot(ax, geofence_lon, geofence_lat, 'r-', 'LineWidth', 2);

% Initialize trajectory data
trajectory_x = lon_drone;
trajectory_y = lat_drone;
trajectory_z = alt_drone; % Track altitude for ascent phase

% Create trajectory line (but no data yet)
%trajectory_line = plot3(trajectory_x, trajectory_y, trajectory_z, 'k--', 'LineWidth', 1.3);
%trajectory_line = plot3(trajectory_x, trajectory_y, trajectory_z, 'k-.-', 'LineWidth', 1.3);
trajectory_line = plot3(trajectory_x, trajectory_y, trajectory_z, 'k-.', 'LineWidth', 1.3);


%% Plot the eNodeBs
%plot(ax, lon_eNBs, lat_eNBs, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'blue');
%hold on; % Ensure all points are plotted on the same figure
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


%% geofence 
% Function to check if a point is inside the geofence
is_in_geofence = @(lat, lon) inpolygon(lon, lat, geofence_lon, geofence_lat);

%% Define minimum safe distance (20 meters)
min_safe_distance = 0.00040; % approx 50m
threshold_distance = 0.0020; %.00250m
initial_step_size = 0.00001; % Initial speed of movement

%% Additional features for horizontal fly
%initial_speed = 10;
max_speed = 10;
min_speed = 5;
% Initialize speed and position
%current_speed = initial_speed;
current_speed = max_speed;
uav_speed = max_speed;
current_position = 0;
time_step = 0.1; % Time step in seconds

% Initialize variables
max_attempts = 20; % Maximum attempts to find a valid position
attempt_count = 0; % Track the number of failed attempts
rotation_angles = [0, pi/4, -pi/4, pi/2, -pi/2]; % Angles to rotate by (90 degrees clockwise, counterclockwise)
%rotation_angles = [pi/4] %, pi/4, pi/4, pi/4]; % Angles to rotate by (90 degrees clockwise, counterclockwise)
max_rotation_attempts = length(rotation_angles); % Max rotations to try for each attempt

% Start timer for flight time
flight_time = 180; % Total flight time in seconds (1 minute)
start_time = tic; % Start the timer


%% Report
logFileIDs = cell(4, 1); %FOR rsrp, snr -- theoretical
logFileIDs_bias = cell(4, 1); %FOR rsrp, snr -- theoretical
for k = 1:4
    logFileIDs{k} = fopen(['LW' num2str(k) '_log.txt'], 'w');
    fprintf(logFileIDs{k}, 'time, longitude, latitude, altitude, speed, rsrp, snr, datarate, Prx\n');

    %logFileIDs_bias{k} = fopen(['LW' num2str(k) '_log_bias.txt'], 'w');
    %fprintf(logFileIDs_bias{k}, 'time, Longitude, Latitude, Altitude, rsrp\n');
end


%% Haversine formula
haversine = @(lat1, lon1, lat2, lon2) 2 * 6371000 * ...
    asin(sqrt(sin(deg2rad(lat2 - lat1) / 2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(deg2rad(lon2 - lon1) / 2)^2));

global period_time 
period_time = 1
%% Define a timer object for each function
timer1 = timer('ExecutionMode', 'fixedRate', 'Period',period_time, ...
                'TimerFcn', @(~,~) logUAVCommMetrics(logFileIDs)); % Executes every second
timer2 = timer('ExecutionMode', 'fixedRate', 'Period', 1, ...
                'TimerFcn', @(~,~) my_function2()); % Executes every second

cleanupObj = onCleanup(@() delete_timers(timer1, timer2));

% Start the timers
start(timer1);
   


%% Ascend loop
% Initialize ascent parameters
target_altitude = 30; % Target altitude in meters
current_altitude = 0; % Starting altitude
ascent_completed = false;
time_step = 1; % One-second intervals
base_increment = [0.07, 0.66, 2.25, 4.64, 7.17, 9.67, 12.15, 14.64, 17.12, 19.61, 22.11, 24.58, 26.78, 28.39, 29.44]; % Base pattern of altitude increments
% Initialize position change parameters for lat/lon
% lat_increment = 0.0000000157;  % Small change in latitude
% lon_increment = 0.0000000157;  % Small change in longitude
%                %35.0000003157

%len = length(base_increment)
i=1;
while ~ascent_completed    
    %alt_drone = current_altitude;
    % Randomize each increment slightly around the base pattern values
    %if i<len
        fprintf('%d %.4f',i,base_increment(i));
        %alt_increment = base_increment(i) ;%+ rand;
        if i==1
            current_altitude = rand;
            i = i+1;
        else
            current_altitude = current_altitude + 2 + rand
        end
        %alt_drone = alt_increment;
      %  i=i+1;
    %else
     %   alt_increment = 30;
      %  alt_drone = alt_increment;
    %end

    % Update the current altitude
    %current_altitude =  alt_increment;
    %alt_drone = current_altitude;
    
    % Ensure the UAV doesn't exceed the target altitude
    if current_altitude >= target_altitude
        current_altitude = target_altitude;
       % alt_drone = current_altitude;
        ascent_completed = true;
    end

    alt_drone =  current_altitude

 %   alt_drone = current_altitude
    %pause(.5)
    % Update the UAV's plot position
    % Update latitude and longitude for horizontal movement
    % lat_drone = lat_drone + lat_increment;  % Increment latitude
    % lon_drone = lon_drone + lon_increment;  % Increment longitude
    set(dronePlot, 'XData', lon_drone, 'YData', lat_drone, 'ZData', alt_drone);
    alt_drone =  current_altitude

    % Update trajectory
    trajectory_x(end + 1) = lon_drone; 
    trajectory_y(end + 1) = lat_drone;
    trajectory_z(end + 1) = alt_drone; % Add the altitude

    % Update trajectory line plot
    set(trajectory_line, 'XData', trajectory_x, 'YData', trajectory_y, 'ZData', trajectory_z);    
    alt_drone =  current_altitude

    % Apply axis limits
    xlim([-78.701 -78.6914]);
    ylim([35.722 35.7339]);
    xtickangle(45);

    % Update the plot
    drawnow;
    
    %pause(1); % Pause for 1 second to simulate real-time ascent
    alt_drone =  current_altitude
    % Print the altitude for debugging purposes
    fprintf('Current Altitude: %.2f meters\n', current_altitude);
    pause(.75)
end

disp('UAV reached 30 meters altitude. Starting horizontal movement to base station...');


%% UAV fly horizontally

%start(timer2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move UAV towards base stations
while ~isempty(lat_eNBs) % While there are base stations left
    global uav_speed
    current_position = 0;
    % Check if the flight time has exceeded 1 minute
    elapsed_time = toc(start_time);
    if elapsed_time >= flight_time
        disp('Flight time exceeded. Stopping UAV...');
        break; % Break out of the loop to stop the UAV
    end
    
    % Calculate distance from UAV to all base stations
    distances = pdist2([lat_drone, lon_drone], [lat_eNBs(:), lon_eNBs(:)]); %euclidean distance by default
    
    % Find the nearest base station
    [min_distance, nearest_idx] = min(distances);
    %[min_distance, nearest_idx] = max(distances);
    disp(['Visited a lakewheel................... ', num2str(min_distance)]);
    %fprintf('Visited a lakewheel................... %.4f\n', min_distance)

    % Get the coordinates of the nearest base station
    next_lat = lat_eNBs(nearest_idx);
    next_lon = lon_eNBs(nearest_idx);

    if next_lat == 35.733056 && next_lon == -78.698333
        threshold = 0.0030; %280
    else
        threshold = threshold_distance;
    end

    
    % Calculate the direction to move in (lat, lon)
    lat_diff = next_lat - lat_drone;
    lon_diff = next_lon - lon_drone;

    % Normalize the direction vector
    distance = sqrt(lat_diff^2 + lon_diff^2);

    % Gradually reduce the step size as time progresses (slowing down the UAV)
    remaining_time = flight_time - elapsed_time; % Remaining flight time

    % Calculate distance from UAV to the nearest base station
    distance_to_nearest_eNB = min_distance;
    % If the UAV is far from the base station (>100 meters), increase speed

    %remaining_distance = distance_to_target - current_position;
    remaining_distance = distance_to_nearest_eNB; % - current_position;

    if distance_to_nearest_eNB > 0.0009 % Approx 100 meters
        %step_size = initial_step_size * (distance_to_nearest_eNB / 0.001) % Increase speed when far away
        initial_speed = max_speed;
        current_speed = initial_speed %(remaining_distance / distance_to_nearest_eNB) * initial_speed % Scale speed to remaining distance
        uav_speed = current_speed;
        step_size = current_speed * initial_step_size; % Step size in meters
    else
        % If the UAV is near the base station (within 100 meters), slow down
        %step_size = initial_step_size * (distance_to_nearest_eNB / 0.0009) % Slow down as it gets closer
        initial_speed = min_speed; %max_speed /2
        current_speed = initial_speed %(remaining_distance / distance_to_nearest_eNB) * initial_speed; % Scale speed to remaining distance
        uav_speed = current_speed;
        step_size = current_speed * initial_step_size; % Step size in meters
        current_position = current_position + step_size;
    end
    %pause(10)

    %step_size = initial_step_size * (remaining_time / flight_time); % Reduce speed as time progresses

    % Check if the distance is greater than the minimum safe distance
    if distance > min_safe_distance

        previous_distance = Inf;  % Initialize the previous distance with a large value
        % Move UAV in steps towards the base station
        %while distance > min_safe_distance
        while distance < previous_distance
            % Calculate movement proportions
            
            % remaining_distance = distance;
            % current_speed = (remaining_distance / distance_to_nearest_eNB) * initial_speed % Scale speed to remaining distance
            % step_size = current_speed * initial_step_size % Step size in meters
            % current_position = current_position + step_size;

            lat_step = (lat_diff / distance) * step_size;
            lon_step = (lon_diff / distance) * step_size;

            % Check if the new position would be within the geofence
            new_lat_drone = lat_drone + lat_step;
            new_lon_drone = lon_drone + lon_step;

            if is_in_geofence(new_lat_drone, new_lon_drone)
                %disp('uav in geofence')
                % Update UAV position
                lat_drone = new_lat_drone;
                lon_drone = new_lon_drone;

                % Update the UAV's plot position
                set(dronePlot, 'XData', lon_drone, 'YData', lat_drone, 'ZData', alt_drone);

                % Update trajectory
                trajectory_x(end + 1) = lon_drone; 
                trajectory_y(end + 1) = lat_drone;
                trajectory_z(end + 1) = alt_drone; 

                % Update trajectory line plot
                set(trajectory_line, 'XData', trajectory_x, 'YData', trajectory_y,'ZData', trajectory_z);
                %daspect([1 1 1]); 


                % Update the plot
                drawnow; % Ensure the plot updates immediately
                pause(0.1); % Adjust the pause for faster or slower movement

                % Update previous distance with the current distance
                previous_distance = distance;

                % Recalculate distance to the base station
                lat_diff = next_lat - lat_drone;
                lon_diff = next_lon - lon_drone;
                distance = sqrt(lat_diff^2 + lon_diff^2);
                
                %dynamic_hover_threshold = max(0.0009, distance * 0.8) % Example: set it to 80% of the current distance

                %distance_to_nearest_eNB = min_distance;
                % If the UAV is far from the base station (>100 meters), increase speed
                if distance > 0.0009 % Approx 100 meters
                    remaining_distance = distance;
                    initial_speed = max_speed;
                    current_speed = initial_speed %(remaining_distance / distance_to_nearest_eNB) * initial_speed; % Scale speed to remaining distance
                    uav_speed = current_speed;
                    step_size = current_speed * initial_step_size; % Step size in meters
                    current_position = current_position + step_size;

                 %   step_size = initial_step_size * (distance / 0.001) % Increase speed when far away
                else
                    initial_speed = min_speed;
                    remaining_distance = distance;
                    current_speed = initial_speed %(remaining_distance / distance_to_nearest_eNB) * initial_speed; % Scale speed to remaining distance
                    uav_speed = current_speed;
                    step_size = current_speed * initial_step_size; % Step size in meters
                    current_position = current_position + step_size;
                end

                % Hover if the distance is no longer decreasing (i.e., UAV has arrived)
                if distance >= previous_distance %&& distance < threshold 
                    % Set speed to zero to indicate hovering

                    if distance < threshold 
                        step_size = 0; % Reduce speed to zero when hovering
                        speed = 0; % Update speed to zero
    
                        disp('UAV is near the base station. Hovering for 5 seconds...');                    
                        pause(.5); % Hover for .5 seconds
                        lat_eNBs(nearest_idx) = []; % Remove the base station from the list
                        lon_eNBs(nearest_idx) = [];     
                        current_position = 0;
                        break;
                    % else
                    %     %pause(2)
                    %     % Rotate UAV by the next angle in the list
                    %     attempt_count = attempt_count + 1
                    %     rotation_angle = rotation_angles(attempt_count);
                    %     if attempt_count == 4
                    %         attempt_count = 0;
                    %     end
                    %     temp_lat = lat_diff;
                    %     lat_diff = temp_lat * cos(rotation_angle) - lon_diff * sin(rotation_angle); % Rotate
                    %     lon_diff = temp_lat * sin(rotation_angle) + lon_diff * cos(rotation_angle);
                    % 
                    %     % Reset distance for the new direction
                    %     distance = sqrt(lat_diff^2 + lon_diff^2); 
                    %     disp('its mmmmmm problme................')
                    end                    
                end

                % Check if the flight time has exceeded 1 minute during movement
                if toc(start_time) >= flight_time
                    disp('Flight time exceeded during movement. Stopping UAV...');
                    break; % Break the inner loop and stop movement
                end
            else
                % If out of geofence, increase the attempt count %% Handle out-of-geofence cases
                attempt_count = attempt_count + 1;
                disp('UAV is out of geofence. change orientation');

                % Rotate UAV by the next angle in the list
                rotation_angle = rotation_angles(attempt_count);
                if attempt_count == 4
                    attempt_count = 0;
                end
                temp_lat = lat_diff;
                lat_diff = temp_lat * cos(rotation_angle) - lon_diff * sin(rotation_angle); % Rotate
                lon_diff = temp_lat * sin(rotation_angle) + lon_diff * cos(rotation_angle);
                
                % Reset distance for the new direction
                distance = sqrt(lat_diff^2 + lon_diff^2); 


                
            end
        end
    else
        % If within the minimum safe distance, pause for 5 seconds
        disp('UAV is within the minimum safe distance of the base station. Pausing for 5 seconds...');
        %pause(5); % Pause for 5 seconds

        % Remove the current base station from the list
        lat_eNBs(nearest_idx) = [];
        lon_eNBs(nearest_idx) = [];
        
        % Rotate slightly to move towards the next base station
        disp('Changing heading towards the next base station...');
        attempt_count = 0; % Reset the attempt count for the next base station

        % Rotate UAV by the next angle in the list
        rotation_angle = rotation_angles(1); % Reset to the first rotation angle
        temp_lat = lat_diff;
        lat_diff = temp_lat * cos(rotation_angle) - lon_diff * sin(rotation_angle);
        lon_diff = temp_lat * sin(rotation_angle) + lon_diff * cos(rotation_angle);
        
        % Reset distance for the new direction
        distance = sqrt(lat_diff^2 + lon_diff^2); 

        % Move UAV in the new direction
        lat_step = (lat_diff / distance) * step_size; % Move in the new direction
        lon_step = (lon_diff / distance) * step_size;

        new_lat_drone = lat_drone + lat_step;
        new_lon_drone = lon_drone + lon_step;

        % Check if the new position is within the geofence
        if is_in_geofence(new_lat_drone, new_lon_drone)
            % Update UAV's plot position
            lat_drone = new_lat_drone;
            lon_drone = new_lon_drone;

            % Update trajectory
            trajectory_x(end + 1) = lon_drone; 
            trajectory_y(end + 1) = lat_drone;
            trajectory_z(end + 1) = alt_drone; % Add the altitude

            % Update trajectory line plot
            set(trajectory_line, 'XData', trajectory_x, 'YData', trajectory_y, 'ZData', trajectory_z);    
               
            % Update the plot
            drawnow; % Ensure the plot updates immediately
            pause(0.1); % Adjust the pause for faster or slower movement
        else
            disp('Failed to find valid position after heading adjustment.');
        end
    end
end

disp('UAV has completed its flight and is returning to the launch position.');



% Return to the launch position
lat_diff_return = launch_lat - lat_drone;
lon_diff_return = launch_lon - lon_drone;
distance_return = sqrt(lat_diff_return^2 + lon_diff_return^2);
attempt_count = 0;
rnd = 0;
initial_step_size = 0.00005
% Move UAV back to the launch position
while distance_return > 0.00001 %0.0001  % Continue until very close to launch
    global uav_speed
    % Calculate movement step towards the launch position
    dynamic_hover_threshold = max(0.0009, distance_return * 0.8) % Example: set it to 80% of the current distance

                %distance_to_nearest_eNB = min_distance;
                % If the UAV is far from the base station (>100 meters), increase speed
                if distance_return > 0.0009 % Approx 100 meters
                    current_speed = distance_return / 0.001;
                    uav_speed = current_speed;
                    step_size = initial_step_size * current_speed ;%(distance_return / 0.001); % Increase speed when far away
                else
                    % If the UAV is near the base station (within 100 meters), slow down
                    %step_size = initial_step_size * (distance / 0.0009); % Slow down as it gets closer
                    current_speed = distance_return / dynamic_hover_threshold;
                    uav_speed = current_speed;
                    step_size = initial_step_size * current_speed; %(distance_return / dynamic_hover_threshold); % Slow down as it gets closer
                end

    lat_step_return = (lat_diff_return / distance_return) * step_size;
    lon_step_return = (lon_diff_return / distance_return) * step_size;

    % Predict new position
    new_lat_drone = lat_drone + lat_step_return;
    new_lon_drone = lon_drone + lon_step_return;

    % Check if the new position is within the geofence
    if is_in_geofence(new_lat_drone, new_lon_drone)
        % Update UAV position
        lat_drone = new_lat_drone;
        lon_drone = new_lon_drone;

        % Update UAV position in the plot
        set(dronePlot, 'XData', lon_drone, 'YData', lat_drone, 'ZData', alt_drone);

        % Update trajectory
        trajectory_x(end + 1) = lon_drone; 
        trajectory_y(end + 1) = lat_drone;
        trajectory_z(end + 1) = alt_drone; % Add the altitude

        % Update trajectory line plot
        set(trajectory_line, 'XData', trajectory_x, 'YData', trajectory_y, 'ZData', trajectory_z);    
            
        
        % Update the plot
        drawnow;
        pause(0.1);  % Adjust the pause as necessary

        % Recalculate distance to launch position
        lat_diff_return = launch_lat - lat_drone;
        lon_diff_return = launch_lon - lon_drone;
        distance_return = sqrt(lat_diff_return^2 + lon_diff_return^2);

        % Debugging output
        fprintf('UAV Position: (%.6f, %.6f), Distance to Launch: %.6f\n', lat_drone, lon_drone, distance_return);
    else
        % If outside geofence, handle orientation change
        disp('UAV is outside the geofence during return. Adjusting course...');

        % Increment attempt count and reset if it exceeds available rotations
        attempt_count = mod(attempt_count, numel(rotation_angles)) + 1;
        rotation_angle = rotation_angles(attempt_count);

        % Adjust heading towards launch point more aggressively
        heading_angle = atan2(lat_diff_return, lon_diff_return); % Calculate angle towards launch
        rotation_angle = heading_angle + rotation_angles(attempt_count); % Add rotation to new heading

        % Predict new position based on adjusted heading
        lat_step_return = cos(rotation_angle) * step_size;
        lon_step_return = sin(rotation_angle) * step_size;

        % Predict new position
        if rnd == 0
            new_lat_drone = lat_drone - lat_step_return;
            new_lon_drone = lon_drone - lon_step_return;
            rnd = 1;
        else
            new_lat_drone = lat_drone + lat_step_return;
            new_lon_drone = lon_drone + lon_step_return;
            rnd = 0;
        end

        % Check if the new position is within the geofence
        if is_in_geofence(new_lat_drone, new_lon_drone)
            % Update UAV position
            lat_drone = new_lat_drone;
            lon_drone = new_lon_drone;

            % Update UAV position in the plot
            set(dronePlot, 'XData', lon_drone, 'YData', lat_drone, 'ZData', alt_drone);

            % Update trajectory
            trajectory_x(end + 1) = lon_drone; 
            trajectory_y(end + 1) = lat_drone;
            trajectory_z(end + 1) = alt_drone; % Add the altitude

            % Update trajectory line plot
            set(trajectory_line, 'XData', trajectory_x, 'YData', trajectory_y, 'ZData', trajectory_z);    
                        
            % Update the plot
            drawnow;
            pause(0.1);  % Adjust the pause as necessary

            % Recalculate distance to launch position
            lat_diff_return = launch_lat - lat_drone;
            lon_diff_return = launch_lon - lon_drone;
            distance_return = sqrt(lat_diff_return^2 + lon_diff_return^2);    
            % Optional: Set step size based on distance to nearest base station if needed
            if distance_return > 0.001 % Approx 100 meters
                current_speed = distance_return / 0.001; 
                uav_speed = current_speed;
                step_size = initial_step_size * current_speed ; %(distance_return / 0.001); % Increase speed when far away
            else
                % If the UAV is near the base station (within 100 meters), slow down
                current_speed = distance_return / 0.0009; 
                uav_speed = current_speed;
                step_size = initial_step_size * current_speed; %(distance_return / 0.0009); % Slow down as it gets closer
            end
        end
    end
end




%% UAV Descending vertically

% Define the base increment pattern and scaling factor
min_altitude = 0; % Minimum altitude to reach during descent
%time_step = 1; % One-second intervals

current_altitude = alt_drone; % Starting altitude (after ascent)
ascent_completed = true; % Start descent after ascent is completed

% Initialize trajectory data
%trajectory_x = lon_drone;
%trajectory_y = lat_drone;
%trajectory_z = alt_drone; % Start with initial altitude

% Create trajectory line (but no data yet)
%trajectory_line = plot3(trajectory_x, trajectory_y, trajectory_z, 'r-', 'LineWidth', 1);

% Descent loop
while current_altitude > min_altitude    
    current_altitude = current_altitude -2 - rand; % - rand % Decrement altitude
    %current_altitude = alt_increment;     
    
    % Update UAV position and altitude
    alt_drone = current_altitude
    
    % % Update the UAV's plot position
    % lat_drone = lat_drone - lat_increment;  % Increment latitude
    % lon_drone = lon_drone - lon_increment;  % Increment longitude   
    % if lat_drone <= 35.7274803000
    %     lat_drone = 35.7274803000
    % elseif lon_drone <= -78.6962752000
    %     lon_drone = -78.6962752000
    % end
    set(dronePlot, 'XData', lon_drone, 'YData', lat_drone, 'ZData', alt_drone);

    % Update trajectory data
    trajectory_x(end + 1) = lon_drone;
    trajectory_y(end + 1) = lat_drone;
    trajectory_z(end + 1) = alt_drone; % Add updated altitude to trajectory

    % Update trajectory line plot
    set(trajectory_line, 'XData', trajectory_x, 'YData', trajectory_y, 'ZData', trajectory_z);

    % Update the plot
    drawnow;
    pause(1); % Pause for 1 second to simulate real-time descent

    % Print the altitude for debugging purposes
    fprintf('Current Altitude: %.2f meters\n', current_altitude);
end

disp('UAV reached minimum altitude. Descent completed.');

% Set axes to have equal scaling
axis equal;

set(gcf, 'PaperUnits', 'inches'); % Set paper units to inches
set(gcf, 'PaperSize', [6.5 4.5]); % Set to standard letter size (width, height)
set(gcf, 'PaperPosition', [0 0 8.5 11]); % Ensure the figure fits within the page

% Save the figure as a PDF
print(gcf, 'autonomous.pdf', '-dpdf', '-r300', '-bestfit');
print(gcf, 'autonomous.png', '-dpng', '-r300');
saveas(gcf, 'autonomous.fig');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Asynchronous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% programming
% Stop the timers after completing the movement
% Close log files
for k = 1:4
    fclose(logFileIDs{k});
    %fclose(logFileIDs_bias{k}); %traffic
%    fclose(logTrafficIDs{k}); %traffic
end

stop(timer1);
stop(timer2);
delete(timer1);
delete(timer2);

disp('All tasks completed.');

% Function: logUAVCommMetrics
%log: The function logs data to a file.
%UAV: The context involves a UAV (drone).
%Comm: Refers to the communication (signal strength, noise, etc.).
%Metrics: The function calculates and logs key metrics such as RSRP, SNR, and data rate.
function logUAVCommMetrics(logFileIDs)
            
    haversine = @(lat1, lon1, lat2, lon2) 2 * 6371000 * ...
        asin(sqrt(sin(deg2rad(lat2 - lat1) / 2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(deg2rad(lon2 - lon1) / 2)^2));

    % AWGN Channel Simulation and SNR, RSRP Calculation
    % Simulation parameters
    transmit_power = 10; % dBm
    noise_power = -90; % dBm
    
    % Antenna gains 
    antenna_gain_uav = 2; % dBi  omnidirectional antenna
    antenna_gain_enb = 10; % dBi directional antenna

    lat_eNBs = [35.7275, 35.728056, 35.725, 35.733056]
    lon_eNBs = [-78.695833, -78.700833, -78.691667, -78.698333];

    global lat_drone lon_drone uav_speed alt_drone;


    lat_points = lat_drone;
    lon_points = lon_drone;
    speed_uav = uav_speed;
    alt_points = alt_drone;
    alt_eNBs = 10;
    gps_data = [];
    data_rate = [];
    %pause(.0001)
    for k = 1:4
                %2D Distance
                %distance_to_enb = haversine(lat_points, lon_points, lat_eNBs(k), lon_eNBs(k))
                %3D Distance
                distance_to_enb = haversine3d(lat_points, lon_points, alt_points, lat_eNBs(k), lon_eNBs(k),alt_eNBs);                

                path_loss = 20 * log10(distance_to_enb) + 20 * log10(3410*10^6) - 147.55 ; %+ randn * 5; %This formula computes the path loss using the free-space path loss model. The constant factors relate to frequency and environmental effects
                %offset =0;
                if k==1
                    offset = 8;
                elseif k==2
                    offset = 5;
                elseif k==3
                    offset = 6;
                elseif k==4
                    offset = 7;
                end
                
                offset = 0;
                
                RSRP = ceil(transmit_power + antenna_gain_enb - path_loss + antenna_gain_uav + offset);
                %awgn_noise = (randn) * sqrt(0.5 * 10^(noise_power / 10))%% Convert noise power to watts
                received_power = RSRP ; %+ awgn_noise;
               % received_power = received_power - 30
                snr_gap = -5;
                SNR = received_power - noise_power + snr_gap;

                C = calculateDataRate(received_power, noise_power)
                %C = 1.4 * log2(1 + 10^(SNR/10)); % Capacity in bits per second (bps) %data rate
                %C = C * 0.25


                %fprintf('Theoretical Maximum Data Rate (Shannon Capacity): %.2f Mbps\n', C_mbps);
                timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
                gps_data = [gps_data ; sprintf('%.6f,%.6f', lat_points, lon_points)];
                fprintf(logFileIDs{k}, '%s, %.12f, %.12f, %.12f, %.12f, %.12f, %.12f, %.12f, %.12f\n', timestamp, lon_points, lat_points, alt_points, speed_uav, RSRP, SNR, C, received_power);
                %fprintf('%s, %.6f, %.6f, %.2f, %.2f, %.2f, %.2f\n', timestamp, lon_points, lat_points, alt_points, RSRP, SNR, C);               
    end

    disp('Executing my_function...');
    
end

function my_function2()
    % Your code for my_function2
    disp('Executing my_function2...');
end

function delete_timers(timer1, timer2)
    % Stop and delete timers
    if strcmp(timer1.Running, 'on')
        stop(timer1);
    end
    delete(timer1);
    
    if strcmp(timer2.Running, 'on')
        stop(timer2);
    end
    delete(timer2);
    
    disp('Timers stopped and deleted.');
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
    c = 2 * R * asin(sqrt(a));

    % 3D distance in meters
    distance = sqrt((c)^2 + dAlt^2);
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

function se = getSpectralEfficiency(snr)
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
