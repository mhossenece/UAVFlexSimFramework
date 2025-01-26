function all_waypoints_in_csv(savTrajWayps, traj, num_waypoints)

    filename = 'waypoints_all_trajectory.csv';

    header = 'latitude,longitude,altitude';
    fid = fopen(filename, 'w');

    if fid == -1
        error('Cannot open file: %s', filename);
    end

    fprintf(fid, '%s\n', header);

    data_idx = 1;
    j = 1;
    x = 0;

    for i = 1:(size(savTrajWayps, 1) + traj)
        i = i + x;
        if i == 1 
            header = sprintf('Path %d, ......, .......', j);
            fprintf(fid,'%s\n',header);
            j = j + 1;   

        elseif mod(i, num_waypoints+4) == 0            
            header = sprintf('Path %d, ......, .......', j);
            fprintf(fid,'%s\n',header);
            j = j + 1;
            x = x + 1;
        else
            fprintf(fid, '%.6f,%.6f,%.6f\n', savTrajWayps(data_idx, :));
            data_idx = data_idx + 1;
        end
    end
    fclose(fid);
end