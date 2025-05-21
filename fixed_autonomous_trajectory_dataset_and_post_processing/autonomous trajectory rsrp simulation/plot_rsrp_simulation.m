clc;
clear;
close all;

% Load data from CSV files for LW1 to LW4
% data_LW1 = readtable('lw1_merged_icc.csv');
% data_LW2 = readtable('lw2_merged_icc.csv');
% data_LW3 = readtable('lw3_merged_icc.csv');
% data_LW4 = readtable('lw4_merged_icc.csv');

data_LW1 = readtable('LW1_log.csv');
data_LW2 = readtable('LW2_log.csv');
data_LW3 = readtable('LW3_log.csv');
data_LW4 = readtable('LW4_log.csv');

% Extract latitude, longitude, and rsrp for each LW
latitude_LW1 = data_LW1.latitude;
longitude_LW1 = data_LW1.longitude;
rsrp_LW1 = data_LW1.rsrp;

latitude_LW2 = data_LW2.latitude;
longitude_LW2 = data_LW2.longitude;
rsrp_LW2 = data_LW2.rsrp;

latitude_LW3 = data_LW3.latitude;
longitude_LW3 = data_LW3.longitude;
rsrp_LW3 = data_LW3.rsrp ;

latitude_LW4 = data_LW4.latitude;
longitude_LW4 = data_LW4.longitude;
rsrp_LW4 = data_LW4.rsrp;


% LW positions
lat_eNBs = [35.7275, 35.728056, 35.725, 35.733056];
lon_eNBs = [-78.695833, -78.700833, -78.691667, -78.698333];

% Create a figure with a 1x4 layout
figure('Position', [150, 150, 1240, 300]); % Keep the figure height as is
t = tiledlayout(1, 4, 'TileSpacing', 'tight', 'Padding', 'tight');

% Define the colormap and color range
caxis_range = [min([rsrp_LW1; rsrp_LW2; rsrp_LW3; rsrp_LW4]), max([rsrp_LW1; rsrp_LW2; rsrp_LW3; rsrp_LW4])];
caxis_range1 = [min([rsrp_LW1]), max([rsrp_LW1])];
caxis_range2 = [min([rsrp_LW2]), max([rsrp_LW2])];
caxis_range3 = [min([rsrp_LW3]), max([rsrp_LW3])];
caxis_range4 = [min([rsrp_LW4]), max([rsrp_LW4])];

% Function to plot a colored line
plot_colored_line = @(x, y, c) surface([x'; x'], [y'; y'], zeros(2, length(x)), [c'; c'], ...
                                       'EdgeColor', 'interp', 'LineWidth', 2, 'FaceColor', 'none');

add_caption = @(index) text(mean(xlim), min(ylim) - 0.003,     sprintf('(%c) LW%d', char('a' + index - 1),index),     'FontSize', 10, 'FontWeight', 'normal',     'HorizontalAlignment', 'center',    'Units', 'data', 'VerticalAlignment', 'top');
yangle = 0
xangle = 0 ;%15
rotation = 38

% LW1 subplot (a)
ax1 = nexttile;
hold on;
plot_colored_line(longitude_LW1, latitude_LW1, rsrp_LW1);
colors = ["k", "k", "k", "k"];
legend_entries = cell(1, 4); % Allocate 5 entries: 4 for eNodeBs + 1 for geofence
%colors = ["m", "g", "b", "k"];
%colors = ["r", "g", "b", "k"]; % High-contrast colors
% colors = [
%     0.5, 0, 0;   % Dark Red
%     0.5, 0, 0.5; % Dark Magenta
%     0, 0, 0.5;   % Dark Blue
%     0, 0, 0      % Black
% ];
h = gobjects(1, 4);  % Initialize array to store plot handles
for i = 1:4
        if i==1
        sym = "s"; %square
    elseif i==2
        sym = "^"; %triangle
    elseif i==3
        sym = "d"; %diamond
    elseif i==4
        sym = "x"; %cross
        end    
    %plot(lon_eNBs(i), lat_eNBs(i), colors(i)+sym, 'MarkerSize', 10) %, 'MarkerEdgeColor', 'k');
    h(i) = plot(lon_eNBs(i), lat_eNBs(i), colors(i) + sym, 'MarkerSize', 8); 
    %h(i) = plot(lon_eNBs(i), lat_eNBs(i), 'Marker', sym, 'MarkerSize', 8, 'Color', colors(i, :));

    legend_entries{i} = ['LW' num2str(i)];  % Store legend entry
end
%plot(lon_eNBs, lat_eNBs, '^', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
%text(lon_eNBs, lat_eNBs, {'LW1', 'LW2', 'LW3', 'LW4'}, 'VerticalAlignment', 'bottom');
%text(lon_eNBs, lat_eNBs, {'LW1', 'LW2', 'LW3', 'LW4'}, 'VerticalAlignment', 'bottom',     'FontSize', 14, 'FontWeight', 'normal'); % Adjust FontSize as needed
%set(gca, 'FontSize', 14, 'FontWeight', 'normal');

%xlim([-78.701 -78.6914]);
%ylim([35.724 35.7333]);
xtickangle(xangle);
ytickangle(yangle);
% Plot eNodeBs and their coverage ranges in 3D
%for i = 1:4
%    plot3(ax, lon_eNBs(i), lat_eNBs(i), alt_eNBs(i), colors(i) + sym, 'MarkerSize', 10, 'MarkerFaceColor', colors(i));
 %   legend_entries{i} = ['LW' num2str(i)]; % Store the legend entry for this plot
    %text(ax, lon_eNBs(i), lat_eNBs(i), alt_eNBs(i), ['LW' num2str(i)], 'Color', colors(i), 'FontSize', 10, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
%end
%legend(legend_entries(1:end)); % Exclude the first entry (LW1)
ylabel('Latitude');
%caxis(caxis_range1);
cb = colorbar; %('Direction','reverse');
%cb.Position = cb.Position + [0, 0, 0, -.1]  % [left, bottom, width, height]
cb.Location = 'eastoutside'
cb.Label.String = 'RSRP (dBm)';
%cb.Label.String = 'dBm';
cb.Label.Rotation = rotation; % Make the text horizontal
cb.Label.Position = cb.Label.Position + [-2 -14.5 0]; % Adjust the position

grid on;
xlabel('Longitude');
add_caption(1); % Add caption (a)
% Add legend and associate it with the plot handles
%legend(h, legend_entries, 'Location', 'best','FontSize', 6.5); % Use 'h' to link handles with legend entries

xtickangle(xangle);

axis(ax1, 'equal');



% Calculate combined limits for equal axis length
xLimits1 = xlim(ax1); % Get x-axis limits of subplot 1
yLimits1 = ylim(ax1); % Get y-axis limits of subplot 1


% Find the max range for both axes
xCombinedLimits = [min([xLimits1(1)]), max([xLimits1(2)])];
yCombinedLimits = [min([yLimits1(1)]), max([yLimits1(2)])];

% Ensure equal scaling by extending the shorter axis range
xRange = diff(xCombinedLimits);
yRange = diff(yCombinedLimits);
maxRange = max(xRange, yRange);


% Adjust the range to reduce equally by a percentage (e.g., 10%)
reductionFactor = 0; % 10% reduction
newRange = maxRange * (1 - reductionFactor);


% % Adjust both axes to have equal limits
% xCenter = mean(xCombinedLimits);
% yCenter = mean(yCombinedLimits);
% xNewLimits = [xCenter - maxRange / 2, xCenter + maxRange / 2];
% yNewLimits = [yCenter - maxRange / 2, yCenter + maxRange / 2];

% Calculate new limits for both axes
xCenter = mean(xCombinedLimits);
yCenter = mean(yCombinedLimits);
xNewLimits = [xCenter - newRange / 2, xCenter + newRange / 2];
yNewLimits = [yCenter - newRange / 2, yCenter + newRange / 2];

xNewLimits = xNewLimits
yNewLimits = yNewLimits

% Apply new limits to both subplots
xlim(ax1, xNewLimits);
ylim(ax1, yNewLimits);

%legend(h, legend_entries, 'Location', 'best', 'FontSize', 6.5); % Use 'h' to link handles with legend entries%plot(lon_eNBs, lat_eNBs, '^', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
legend(h, legend_entries, 'Location', 'best', 'FontSize', 6.5,'Box','off','Position', [0.095 0.6 0.2 0.22]); 

% LW2 subplot (b)
ax2 = nexttile;
hold on;
plot_colored_line(longitude_LW2, latitude_LW2, rsrp_LW2);
%colors = ["k", "k", "k", "k"];
legend_entries = cell(1, 4); % Allocate 5 entries: 4 for eNodeBs + 1 for geofence
%colors = ["m", "g", "b", "k"];
h = gobjects(1, 4);  % Initialize array to store plot handles
for i = 1:4
        if i==1
        sym = "s"; %square
    elseif i==2
        sym = "^"; %triangle
    elseif i==3
        sym = "d"; %diamond
    elseif i==4
        sym = "x"; %cross
        end    
    %plot(lon_eNBs(i), lat_eNBs(i), colors(i)+sym, 'MarkerSize', 10) %, 'MarkerEdgeColor', 'k');
    h(i) = plot(lon_eNBs(i), lat_eNBs(i), colors(i) + sym, 'MarkerSize', 8); 
    legend_entries{i} = ['LW' num2str(i)];  % Store legend entry
end
%legend(h, legend_entries, 'Location', 'best'); % Use 'h' to link handles with legend entries
%plot(lon_eNBs, lat_eNBs, '^', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
%text(lon_eNBs, lat_eNBs, {'LW1', 'LW2', 'LW3', 'LW4'}, 'VerticalAlignment', 'bottom', 'FontSize', 14, 'FontWeight', 'normal'); % Adjust FontSize as needed
% Set font size and bold for the axes
%set(gca, 'FontSize', 14, 'FontWeight', 'normal');

%xlim([-78.701 -78.6914]);
%ylim([35.724 35.7333]);
xtickangle(xangle);
%xtickangle(xangle);

axis(ax2, 'equal');



% Calculate combined limits for equal axis length
xLimits1 = xlim(ax2); % Get x-axis limits of subplot 1
yLimits1 = ylim(ax2); % Get y-axis limits of subplot 1


% Find the max range for both axes
xCombinedLimits = [min([xLimits1(1)]), max([xLimits1(2)])];
yCombinedLimits = [min([yLimits1(1)]), max([yLimits1(2)])];

% Ensure equal scaling by extending the shorter axis range
xRange = diff(xCombinedLimits);
yRange = diff(yCombinedLimits);
maxRange = max(xRange, yRange);


% Adjust the range to reduce equally by a percentage (e.g., 10%)
reductionFactor = 0; % 10% reduction
newRange = maxRange * (1 - reductionFactor);


% % Adjust both axes to have equal limits
% xCenter = mean(xCombinedLimits);
% yCenter = mean(yCombinedLimits);
% xNewLimits = [xCenter - maxRange / 2, xCenter + maxRange / 2];
% yNewLimits = [yCenter - maxRange / 2, yCenter + maxRange / 2];

% Calculate new limits for both axes
xCenter = mean(xCombinedLimits);
yCenter = mean(yCombinedLimits);
xNewLimits = [xCenter - newRange / 2, xCenter + newRange / 2];
yNewLimits = [yCenter - newRange / 2, yCenter + newRange / 2];

xNewLimits = xNewLimits
yNewLimits = yNewLimits

% Apply new limits to both subplots
xlim(ax2, xNewLimits);
ylim(ax2, yNewLimits);

%legend(h, legend_entries, 'Location', 'best', 'FontSize', 6.5); % Use 'h' to link handles with legend entries%plot(lon_eNBs, lat_eNBs, '^', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
legend(h, legend_entries, 'Location', 'best', 'FontSize', 6.5,'Box','off','Position', [0.32 0.6 0.2 0.22]); 
%legend(h, legend_entries, 'Location', 'best', 'FontSize', 6.5) %, 'NumColumns', 2); % Use 'h' to link handles with legend entries%plot(lon_eNBs, lat_eNBs, '^', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');

set(gca, 'YTickLabel', []); % Remove y-tick labels
%caxis(caxis_range2);
cb = colorbar; %('Direction','reverse');
%cb.Position = cb.Position + [0, 0, 0, -.1]  % [left, bottom, width, height]
cb.Location = 'eastoutside'
cb.Label.String = 'RSRP (dBm)';
%cb.Label.String = 'dBm';
cb.Label.Rotation = rotation; % Make the text horizontal
%cb.Label.FontSize = 6;
cb.Label.Position = cb.Label.Position + [-2 -18.5 0]; % Adjust the position

grid on;
xlabel('Longitude');
add_caption(2); % Add caption (b)

% LW3 subplot (c)
ax3 = nexttile;
hold on;
plot_colored_line(longitude_LW3, latitude_LW3, rsrp_LW3);
%colors = ["k", "k", "k", "k"];
legend_entries = cell(1, 4); % Allocate 5 entries: 4 for eNodeBs + 1 for geofence
%colors = ["m", "g", "b", "k"];
h = gobjects(1, 4);  % Initialize array to store plot handles
for i = 1:4
        if i==1
        sym = "s"; %square
    elseif i==2
        sym = "^"; %triangle
    elseif i==3
        sym = "d"; %diamond
    elseif i==4
        sym = "x"; %cross
        end    
    %plot(lon_eNBs(i), lat_eNBs(i), colors(i)+sym, 'MarkerSize', 10) %, 'MarkerEdgeColor', 'k');
    h(i) = plot(lon_eNBs(i), lat_eNBs(i), colors(i) + sym, 'MarkerSize', 8); 
    legend_entries{i} = ['LW' num2str(i)];  % Store legend entry
end
%legend(h, legend_entries, 'Location', 'best'); % Use 'h' to link handles with legend entries%plot(lon_eNBs, lat_eNBs, '^', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
%xlim([-78.701 -78.6914]);
%ylim([35.724 35.7333]);
xtickangle(xangle);

axis(ax3, 'equal');



% Calculate combined limits for equal axis length
xLimits1 = xlim(ax3); % Get x-axis limits of subplot 1
yLimits1 = ylim(ax3); % Get y-axis limits of subplot 1


% Find the max range for both axes
xCombinedLimits = [min([xLimits1(1)]), max([xLimits1(2)])];
yCombinedLimits = [min([yLimits1(1)]), max([yLimits1(2)])];

% Ensure equal scaling by extending the shorter axis range
xRange = diff(xCombinedLimits);
yRange = diff(yCombinedLimits);
maxRange = max(xRange, yRange);


% Adjust the range to reduce equally by a percentage (e.g., 10%)
reductionFactor = 0; % 10% reduction
newRange = maxRange * (1 - reductionFactor);


% % Adjust both axes to have equal limits
% xCenter = mean(xCombinedLimits);
% yCenter = mean(yCombinedLimits);
% xNewLimits = [xCenter - maxRange / 2, xCenter + maxRange / 2];
% yNewLimits = [yCenter - maxRange / 2, yCenter + maxRange / 2];

% Calculate new limits for both axes
xCenter = mean(xCombinedLimits);
yCenter = mean(yCombinedLimits);
xNewLimits = [xCenter - newRange / 2, xCenter + newRange / 2];
yNewLimits = [yCenter - newRange / 2, yCenter + newRange / 2];

xNewLimits = xNewLimits
yNewLimits = yNewLimits

% Apply new limits to both subplots
xlim(ax3, xNewLimits);
ylim(ax3, yNewLimits);

%legend(h, legend_entries, 'Location', 'best', 'FontSize', 6.5); % Use 'h' to link handles with legend entries%plot(lon_eNBs, lat_eNBs, '^', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
legend(h, legend_entries, 'Location', 'best', 'FontSize', 6.5,'Box','off','Position', [0.555 0.6 0.2 0.22]); 
%text(lon_eNBs, lat_eNBs, {'LW1', 'LW2', 'LW3', 'LW4'}, 'VerticalAlignment', 'bottom',     'FontSize', 14, 'FontWeight', 'normal'); % Adjust FontSize as needed
set(gca, 'YTickLabel', []); % Remove y-tick labels
caxis(caxis_range3);
cb = colorbar; %('Direction','reverse');
%cb.Position = cb.Position + [0, 0, 0, -.1]  % [left, bottom, width, height]
cb.Location = 'eastoutside'
cb.Label.String = 'RSRP (dBm)';
%cb.Label.String = 'dBm';
cb.Label.Rotation = rotation; % Make the text horizontal
cb.Label.Position = cb.Label.Position + [-2 -21.3 0]; % Adjust the position

grid on;
xlabel('Longitude');
% Set font size and bold for the axes
%set(gca, 'FontSize', 14, 'FontWeight', 'normal');
%cb = colorbar; %('Direction','reverse');

%cb.Position = cb.Position + [0, 0, 0, -.1]  % [left, bottom, width, height]
%cb.Location = 'northoutside'

add_caption(3); % Add caption (c)

% LW4 subplot (d)
ax4 = nexttile;
hold on;
plot_colored_line(longitude_LW4, latitude_LW4, rsrp_LW4);
%colors = ["k", "k", "k", "k"];
legend_entries = cell(1, 4); % Allocate 5 entries: 4 for eNodeBs + 1 for geofence
%colors = ["m", "g", "b", "k"];
h = gobjects(1, 4);  % Initialize array to store plot handles
for i = 1:4
        if i==1
        sym = "s"; %square
    elseif i==2
        sym = "^"; %triangle
    elseif i==3
        sym = "d"; %diamond
    elseif i==4
        sym = "x"; %cross
        end    
    %plot(lon_eNBs(i), lat_eNBs(i), colors(i)+sym, 'MarkerSize', 10) %, 'MarkerEdgeColor', 'k');
    h(i) = plot(lon_eNBs(i), lat_eNBs(i), colors(i) + sym, 'MarkerSize', 8); 
    legend_entries{i} = ['LW' num2str(i)];  % Store legend entry
end

%text(lon_eNBs, lat_eNBs, {'LW1', 'LW2', 'LW3', 'LW4'}, 'VerticalAlignment', 'bottom',     'FontSize', 14, 'FontWeight', 'normal'); % Adjust FontSize as needed
set(gca, 'YTickLabel', []); % Remove y-tick labels
%set(gca, 'FontSize', 14, 'FontWeight', 'normal');
%caxis(caxis_range4);

%xlim([-78.701 -78.6914]);
%ylim([35.724 35.7333]);
%xtickangle(xangle);


% Ensure equal aspect ratio and equal axis lengths
%axis('equal');
axis(ax4, 'equal');



% Calculate combined limits for equal axis length
xLimits1 = xlim(ax4); % Get x-axis limits of subplot 1
yLimits1 = ylim(ax4); % Get y-axis limits of subplot 1


% Find the max range for both axes
xCombinedLimits = [min([xLimits1(1)]), max([xLimits1(2)])];
yCombinedLimits = [min([yLimits1(1)]), max([yLimits1(2)])];

% Ensure equal scaling by extending the shorter axis range
xRange = diff(xCombinedLimits);
yRange = diff(yCombinedLimits);
maxRange = max(xRange, yRange);


% Adjust the range to reduce equally by a percentage (e.g., 10%)
reductionFactor = 0; % 10% reduction
newRange = maxRange * (1 - reductionFactor);


% % Adjust both axes to have equal limits
% xCenter = mean(xCombinedLimits);
% yCenter = mean(yCombinedLimits);
% xNewLimits = [xCenter - maxRange / 2, xCenter + maxRange / 2];
% yNewLimits = [yCenter - maxRange / 2, yCenter + maxRange / 2];

% Calculate new limits for both axes
xCenter = mean(xCombinedLimits);
yCenter = mean(yCombinedLimits);
xNewLimits = [xCenter - newRange / 2, xCenter + newRange / 2];
yNewLimits = [yCenter - newRange / 2, yCenter + newRange / 2];

xNewLimits = xNewLimits
yNewLimits = yNewLimits

% Apply new limits to both subplots
xlim(ax4, xNewLimits);
ylim(ax4, yNewLimits);

%legend(h, legend_entries, 'Location', 'best', 'FontSize', 6.5,'Box','on'); % Use 'h' to link handles with legend entries%plot(lon_eNBs, lat_eNBs, '^', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
legend(h, legend_entries, 'Location', 'best', 'FontSize', 6.5,'Box','off','Position', [0.79 0.6 0.2 0.22]); 

grid on;
xlabel('Longitude');
add_caption(4); % Add caption (d)

% Add colorbar
cb = colorbar; %('Direction','reverse');
%cb.Position = cb.Position + [0, 0, 0, -.1]  % [left, bottom, width, height]
cb.Location = 'eastoutside'

%cb.Position = [0.9616 -.04 0.0129 0.7824]

%westoutside northoutside southoutside
cb.Label.String = 'RSRP (dBm)';
%cb.Label.String = 'dBm';
cb.Label.Rotation = rotation; % Make the text horizontal
cb.Label.Position = cb.Label.Position + [-2 -7 0]; % Adjust the position
% Adjust the colorbar position to move it slightly down (south)

% Set label font size and weight
%cb.Label.FontSize = 14; % Adjust the font size (change to your preference)
cb.Label.FontWeight = 'normal'; % Make the label bold

%colormap(jet);
%colormap('parula');

%colormap('turbo')
colormap(flipud(copper))

%colormap('hsv')

%cp = flipud(hot);
% % Gradually decrease the color intensity to light gray in the last 20 rows
% start_color = cp(1:20, :);  % Get the color at the 20th last position
% end_color = [0.75, 0.75, 0.75];  % Define the target light gray color
% 
% % Use linspace on each RGB component separately
% cp(1:20, 1) = linspace(start_color(1), end_color(1), 20);  % Red channel
% cp(1:20, 2) = linspace(start_color(2), end_color(2), 20);  % Green channel
% cp(end-1:end, 3) = linspace(start_color(3), end_color(3), 20);  % Blue channel
%colormap(flipud(hot))

% Define custom colormap with RGB values
custom_cmap1 = [
    0, 0, 1;    % Blue
    0, 1, 0;    % Green
    1, 1, 0;    % Yellow
    1, 0, 0     % Red
]; 

custom_cmap = [
    0, 0, 0;    % Black
    1, 1, 1;    % White
    1, 0, 0;    % Red
    0, 1, 0;    % Green
    0, 0, 1;    % Blue
    0.5, 0.5, 0; % Yellow (dark)
    0.5, 0, 0.5; % Purple (dark)
];


%colormap(mymap);

%caxis([-75 -44]);
%caxis(caxis_range); % Set the same color limits for the colorbar

%print(gcf, 'your_figure.png', '-dpng', '-r300');
%print(gcf, '1.pdf', '-dpdf', '-r300','-bestfit');
% Adjust paper size to fit the figure properly
% Set the aspect ratio equal
%axis equal;

% Set paper units to inches
set(gcf, 'PaperUnits', 'inches');

% Specify the desired paper size (width, height)
set(gcf, 'PaperSize', [8.5 2.3]); 

% Specify the figure's position on the paper ([left, bottom, width, height])
set(gcf, 'PaperPosition', [0 0 8.5 2.3]); 
%set(gcf, 'PaperPosition', [0 0 8.5 11]); % Ensure the figure fits within the page

% Save the figure as a PDF with high resolution
print(gcf, 'autonomous_rsrp_simulation.pdf', '-dpdf', '-r300', '-bestfit');
