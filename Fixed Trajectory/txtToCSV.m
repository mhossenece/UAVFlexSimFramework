% Specify the input and output file names
clc
clear all
close all
%input_file = 'LW4_log_bias.txt';
%output_file = 'LW4_log_bias.csv';

input_file = 'LW1_log.txt';
output_file = 'LW1_log.csv';

% Open the input file for reading
fid_in = fopen(input_file, 'r');

% Read the first line (header)
header = fgetl(fid_in);

% Open the output file for writing
fid_out = fopen(output_file, 'w');

% Write the header to the output file
fprintf(fid_out, '%s\n', header);

% Read and process each line from the input file
while ~feof(fid_in)
    % Read a line from the input file
    line = fgetl(fid_in);
    
    % Remove square brackets
    cleaned_line = strrep(line, '[', '');
    cleaned_line = strrep(cleaned_line, ']', '');
    
    % Write the cleaned line to the output file
    fprintf(fid_out, '%s\n', cleaned_line);
end

% Close the input and output files
fclose(fid_in);
fclose(fid_out);

disp('File conversion complete.');
