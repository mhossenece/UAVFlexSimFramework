function [waypoints, is_data_collection] = create_waypoints(waypoints,sorted_indices,sorted_lat_eNBs,sorted_lon_eNBs,num_waypoints,geofence_lon,geofence_lat,min_lon,max_lon,min_lat,max_lat,lon_drone,lat_drone,altitude,waypoint_index)

is_in_geofence = @(lat, lon) inpolygon(lon, lat, geofence_lon, geofence_lat);


    for i = 1:length(sorted_lat_eNBs)
            if waypoint_index > num_waypoints + 1
                break;
            end
    
            lat = sorted_lat_eNBs(i);
            lon = sorted_lon_eNBs(i);

            lon2 = lon; lat2 = lat;       
            
            while true
                %ddd = i
                 if is_in_geofence(lat, lon)
                    % waypoint_indexFF = waypoint_index
                     is_data_collection(waypoint_index-1) = true;
                     waypoints(waypoint_index, :) = [lat, lon, altitude];                    
                     waypoint_index = waypoint_index + 1;
                   %
                     %way = waypoint_index
                     break;
                 else
                    % Adjust latitude and longitude based on the LW index
                    switch sorted_indices(i)
                        %i
                        case 4 % LW4
                            %lat = lat - rand() / 1000 - rand/1000;% -rand/9000;
                            lat = 35.73 + rand/700
                            lon = lon + randn/4000;
                            %waypoint_indexss = waypoint_index
                           % is_data_collection(waypoint_index-1) = true;
                           % lon = min_lon + (max_lon - min_lon) * rand();
                            %if lon > lon_drone 
                            %    while true
                             %       lon = min_lon + (max_lon - min_lon) * rand();
                                    %lon = lon_drone + (max_lon - min_lon) * rand() %+ 4 * rand/1000 + rand /1000;
                              %      if lon < lon_drone
                               %         break
                                %    end
                            %    end
                           % end
                        case 3 % LW3
                            %x_threshold = lon_drone;
                            % if sorted_indices(1) == 3 %LW3 as a first waypoint
                            %     x_threshold = 
                            % end
                            lat_LW3 = lat;
                            %lat1 = lat_drone - randn / 1000 -rand/1000;
                          %  sorted_indices(i-1)
                            %if lat1 < lat_drone 
                            if sorted_indices(i-1) ==1   %&& sorted_indices(i-1) ~=4)
                                while true
                                    %lon = min_lon + (max_lon - min_lon) * rand()
                                   %lon =  lon_drone +  rand /1000 + rand/1000 + rand/1000;
                                   lon1 = lon_drone - 1* rand/1000 ;
                                   lat1 = lat_drone - 1* rand/ 1000;
                                   if lon1 < lon_drone && lat1 < lat_drone
                                       %lat1
                                       %lon1
                                        if is_in_geofence(lat1, lon1)
                                        %    w1 = waypoint_index;
                                            waypoints(waypoint_index, :) = [lat1, lon1, altitude];  
                                            is_data_collection(waypoint_index-1) = false;
                                            waypoint_index = waypoint_index + 1;
                                            break
                                        end
                                   end                                   
                                end
                            end
                            %end % end of the waypoint so that it is 
                            % 
                            
                            while true
                                %lat = lat_LW3;
                                %lat = lat - rand / 800;
                                lat = 35.724 + rand/1000;
                                %lon = lon_drone + rand /1000;
                                lon = -78.694 - rand /500;
                                if is_in_geofence(lat, lon)
                                    %w2 = waypoint_index;
                                    is_data_collection(waypoint_index-1) = true;
                                    waypoints(waypoint_index, :) = [lat, lon, altitude];  
                                    waypoint_index = waypoint_index + 1;
                                    break
                                end
                            end
                            
                            if i+1 < 5
                                if sorted_indices(i+1)==1 
                                    while true
                                        %lat = lat_LW3;
                                        lat = lat_drone - rand / 1000;
                                        lon = lon_drone - rand /1000;
                                        if is_in_geofence(lat, lon)
                                            w3 = waypoint_index;
                                            is_data_collection(waypoint_index-1) = false;
                                            waypoints(waypoint_index, :) = [lat, lon, altitude];  
                                            waypoint_index = waypoint_index + 1;                                    
                                            break
                                        end
                                    end
                                end
                            end

                            break

                           % lat = 0


                            % while true
                            %     lat = min_lat + (max_lat - min_lat) * rand();
                            %     lon = min_lon + (max_lon - min_lon) * rand();
                            %     if is_in_geofence(lat, lon)
                            %         waypoints(i, :) = [lat, lon, altitude];
                            %     break;
                            %     end
                            % 
                            %     if lat < lat_drone 
                            %     while true
                            %         %lon = min_lon + (max_lon - min_lon) * rand()
                            %        lon =  lon_drone +  rand /1000 + rand/1000 + rand/1000;
                            %        %lon = lon + randn/9000
                            % 
                            %  
                            % if lon > lon_drone 
                            %          break
                            %        end
                            %     end
                        %    end

                %end
                            
                            %lon = min_lon + (max_lon - min_lon) * rand();
    
                        case 2 % LW2                            
                            lon = lon + rand() / 1000 + rand/ 1000 + rand/ 1000;
                            lat = lat + randn /2000;
%                            is_data_collection(waypoint_index) = true;


                        case 1 % LW1
                            %case1 = waypoint_index 
                            while true
                                lon = lon2 - rand() / 1000;
                                if lon > lon_drone && lon < lon2 && lon < -78.6962
                                    break
                                end
                            end
                            lat = lat_drone + randn/5000;
 %                            while true
 %                                lat = min_lat + (max_lat - min_lat) * rand();
 %                                if lat < 35.7285 && lat > 35.7275
 %                                    break
 %                                end
 %                            end
 % %                           is_data_collection(waypoint_index) = true;
                    end
                end
            end
        end
    
        if waypoint_index <= num_waypoints + 1
            %fprintf('Extra waypoints.');
            for i = waypoint_index:num_waypoints + 1
                while true
                    lat = min_lat + (max_lat - min_lat) * rand();
                    lon = min_lon + (max_lon - min_lon) * rand();
                    if is_in_geofence(lat, lon)
                        waypoints(i, :) = [lat, lon, altitude];
                        is_data_collection(waypoint_index) = false;
                        break;
                    end
                end
            end
        end
end