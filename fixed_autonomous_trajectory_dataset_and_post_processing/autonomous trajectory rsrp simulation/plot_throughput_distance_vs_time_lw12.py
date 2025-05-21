import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.lines as mlines

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great-circle distance between two points 
    on the Earth specified in decimal degrees.
    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371  # Earth radius in km
    return c * r

def convert_time_to_seconds(time_column):
    """
    Convert timestamps to seconds relative to the first timestamp.
    """
    return (time_column - time_column.iloc[0]).dt.total_seconds()

def calculate_distance_to_reference(data, lon_column, lat_column, ref_lon, ref_lat):
    """
    Calculate distance from a fixed reference point (e.g., initial position or BS).
    """
    distances = []
    for i in range(len(data)):
        curr_lon, curr_lat = data[lon_column].iloc[i], data[lat_column].iloc[i]
        distance = haversine(ref_lon, ref_lat, curr_lon, curr_lat) * 1000  # Convert km to meters
        distances.append(distance)
    data['distance_to_reference'] = distances

def plot_scatter(csv_files, x_column, y_column, y1_color, colors, markers, use_bs_reference=False, fig_num=1):
    """
    Plot scatter graphs with throughput vs. time and distance to reference vs. time.
    """
    fig, ax1 = plt.subplots(figsize=(8, 6))  # Increased figure size for better spacing
    #fig, ax1 = plt.subplots(figsize=(8, 7))  # Increased figure size for better spacing
    fontTxt = 22

    # Secondary y-axis for distance
    ax2 = ax1.twinx()
    ax2.set_ylabel('Distance (meters)', color='black', fontsize=fontTxt, fontweight='normal')
    ax2.tick_params(axis='y', labelcolor='black', labelsize=fontTxt)

    handles = []
    distance_handles = []

    # Define launch locations (LW1, LW2, LW3, LW4)
    #lat_LWs = [35.7275, 35.728056, 35.725, 35.733056]  # Latitude for LW1, LW2, LW3, LW4
    #lon_LWs = [-78.695833, -78.700833, -78.691667, -78.698333]  # Longitude for LW1, LW2, LW3, LW4

    lat_LWs = [35.7275, 35.728056] #, 35.725, 35.733056]  # Latitude for LW1, LW2, LW3, LW4
    lon_LWs = [-78.695833, -78.700833] #, -78.691667, -78.698333]  # Longitude for LW1, LW2, LW3, LW4
    
    for idx, csv_file in enumerate(csv_files):
        data = pd.read_csv(csv_file)
        data.columns = data.columns.str.strip()

        if 'time' in data.columns:
            data['time'] = pd.to_datetime(data['time'])
            data['time_in_seconds'] = convert_time_to_seconds(data['time'])
            x_column_internal = 'time_in_seconds' if x_column == 'time' else x_column
        else:
            x_column_internal = x_column

        lon_column, lat_column = "longitude", "latitude"
        if x_column_internal not in data.columns or y_column not in data.columns:
            print(f"Required columns '{x_column_internal}' or '{y_column}' missing in {csv_file}. Skipping.")
            continue

        # Use appropriate reference points for each LW
        if lon_column in data.columns and lat_column in data.columns:
            ref_lon, ref_lat = lon_LWs[idx], lat_LWs[idx]  # For LW1, LW2, LW3, LW4
            calculate_distance_to_reference(data, lon_column, lat_column, ref_lon, ref_lat)
        else:
            print(f"Longitude/latitude columns missing in {csv_file}. Skipping distance calculation.")

        color = colors[idx % len(colors)]
        marker = markers[idx % len(markers)]
        marker_indices = data.index[::1]

        # Plot throughput data
        ax1.scatter(data[x_column_internal].iloc[marker_indices], data[y_column].iloc[marker_indices],
                    color=color, marker=marker, label=f'LW{idx + 1} Throughput', facecolors='none', s=50)

        # Plot distance data with dashed and solid lines
        if idx == 0:  # LW1
            line1, = ax2.plot(data[x_column_internal], data['distance_to_reference'], color='darkblue', linestyle='-', linewidth=1.75)
        elif idx == 1:  # LW2
            line2, = ax2.plot(data[x_column_internal], data['distance_to_reference'], color='darkgreen', linestyle='--', linewidth=1.75)

        legend_marker = mlines.Line2D([], [], color=color, marker=marker, markersize=12, markerfacecolor='none',
                                      linewidth=0)
        handles.append(legend_marker)

        if idx == 0:
            distance_marker = mlines.Line2D([], [], color='darkblue', linestyle='-', linewidth=1.7,
                                            label='Distance to LW1')
            distance_handles.append(distance_marker)
        elif idx == 1:
            distance_marker = mlines.Line2D([], [], color='darkgreen', linestyle='--', linewidth=1.7,
                                            label='Distance to LW2')
            distance_handles.append(distance_marker)

    ax1.set_ylabel("Throughput (Mbits/sec)", color='black', fontsize=fontTxt, fontweight='normal')
    ax1.set_xlabel('Time (seconds)', fontsize=fontTxt, fontweight='normal')
    ax1.tick_params(axis='y', labelcolor='black', labelsize=fontTxt)
    ax1.tick_params(axis='x', labelcolor='black', labelsize=fontTxt)

    # Add arrows and labels to point to LW1, LW2, LW3, LW4 results
    if fig_num == 1:
        #ax1.annotate('Distance to LW1', xy=(1, 4.9), xytext=(1, 5.2), arrowprops=dict(facecolor='darkblue', arrowstyle='->', lw=2), fontsize=fontTxt, color='darkblue')        
        ax1.annotate('Distance to LW1', xy=(62.5, 5), xytext=(68, 4.75), arrowprops=dict(facecolor='darkblue', edgecolor='darkblue', arrowstyle='->', lw=2), fontsize=fontTxt-2, color='darkblue')        
                    #arrowprops=dict(facecolor='black', shrink=0.05), fontsize=fontTxt)
        ax1.annotate('Distance to LW2', xy=(86.7, 7.2), xytext=(85, 7.4), arrowprops=dict(facecolor='darkgreen', edgecolor='darkgreen',arrowstyle='->', lw=2), fontsize=fontTxt-2, color='darkgreen')
        
        #ax1.annotate('TP at LW1', xy=(109, 7.2), xytext=(100, 6.7), arrowprops=dict(facecolor='darkblue', edgecolor='darkblue', arrowstyle='->', lw=2), fontsize=fontTxt, color='darkblue')        
        
        #ax1.annotate('TP at LW2', xy=(107.4, 5.7), xytext=(104, 5.4), arrowprops=dict(facecolor='darkblue', edgecolor='darkblue', arrowstyle='->', lw=2), fontsize=fontTxt, color='darkblue')
    #else:
    #    ax1.annotate('LW3', xy=(60, 20), xytext=(80, 40),
    #                arrowprops=dict(facecolor='black', shrink=0.05), fontsize=fontTxt)
    #    ax1.annotate('LW4', xy=(100, 35), xytext=(120, 50),
    #                arrowprops=dict(facecolor='black', shrink=0.05), fontsize=fontTxt)

    ax1.set_xlim(-3, 165)

    # Combine all handles, including distance lines
    all_handles = handles + distance_handles
    #all_labels = [f'LW{idx + 1} Throughput' for idx in range(len(csv_files))]# + \
                 #['Distance to LW1)', 'Distance to LW2']

    all_labels = [f'TP at LW{idx + 1}' for idx in range(len(csv_files))]# + \
                 #['Distance to LW1)', 'Distance to LW2']
    ax1.legend(handles=all_handles + [line1, line2], labels=all_labels, loc='center right', fontsize=fontTxt-3, frameon=True, ncol=1,
               bbox_to_anchor=(1, .6))

    plt.tight_layout()
    plt.savefig('throughput_distance_vs_time_lw1_lw2.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.show()

# Example usage
if __name__ == "__main__":
    # For LW1 and LW2
    csv_files_LW1_LW2 = ["LW1_log.csv", "LW2_log.csv"]
    colors = ['darkblue', 'darkgreen']
    markers = ['o', 's']
    plot_scatter(csv_files_LW1_LW2, x_column="time", y_column="datarate", y1_color='blue', colors=colors, markers=markers, use_bs_reference=False, fig_num=1)

    # For LW3 and LW4
  #  csv_files_LW3_LW4 = ["LW3_log.csv", "LW4_log.csv"]
  #  colors = ['darkred', 'indigo']
  #  markers = ['^', 'D']
  #  plot_scatter(csv_files_LW3_LW4, x_column="time", y_column="datarate", y1_color='blue', colors=colors, markers=markers, use_bs_reference=False, fig_num=2)
