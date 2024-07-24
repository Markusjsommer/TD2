import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Function to parse the results from a text file
def parse_results(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    data = {
        'Partition': [],
        'User Time (seconds)': [],
        'System Time (seconds)': [],
        'Elapsed Time (seconds)': [],
        'Maximum Memory (kbytes)': []
    }
    
    for i in range(0, len(lines), 6):
        partition_line = lines[i]
        user_time_line = lines[i + 1]
        system_time_line = lines[i + 2]
        elapsed_time_line = lines[i + 3]
        memory_line = lines[i + 4]
        
        partition = partition_line.split(' ')[2].strip(':')
        user_time = float(user_time_line.split(': ')[1])
        system_time = float(system_time_line.split(': ')[1])
        elapsed_time = float(elapsed_time_line.split(': ')[1])
        max_memory = float(memory_line.split(': ')[1])
        
        data['Partition'].append(partition)
        data['User Time (seconds)'].append(user_time)
        data['System Time (seconds)'].append(system_time)
        data['Elapsed Time (seconds)'].append(elapsed_time)
        data['Maximum Memory (kbytes)'].append(max_memory)
    
    return pd.DataFrame(data)

# Parse the results
td_path = 'td_results.txt'
td2_path = 'td2_results.txt'

td_df = parse_results(td_path)
td2_df = parse_results(td2_path)

# Convert partition to numeric and sort
td_df['Partition'] = td_df['Partition'].str.replace('partition_', '').str.replace(':', '').astype(int)
td2_df['Partition'] = td2_df['Partition'].str.replace('partition_', '').str.replace(':', '').astype(int)

td_df = td_df.sort_values('Partition')
td2_df = td2_df.sort_values('Partition')

# Create the Seaborn plots
sns.set_theme(style="whitegrid")

# Function to create and save line plots
def create_line_plot(data1, data2, y, title, ylabel, filename):
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=data1, x='Partition', y=y, marker='o', label='TD')
    sns.lineplot(data=data2, x='Partition', y=y, marker='o', label='TD2')
    plt.title(title)
    plt.xlabel('Partition (%) of MANE')
    plt.ylabel(ylabel)
    plt.xticks(ticks=data1['Partition'], labels=[f'{x}%' for x in data1['Partition']])
    plt.legend()
    plt.savefig(filename)
    plt.show()

# Plot User Time
create_line_plot(td_df, td2_df, 'User Time (seconds)', 'Average User Time (seconds)', 'User Time (seconds)', 'user_time.png')

# Plot System Time
create_line_plot(td_df, td2_df, 'System Time (seconds)', 'Average System Time (seconds)', 'System Time (seconds)', 'system_time.png')

# Plot Elapsed Time
create_line_plot(td_df, td2_df, 'Elapsed Time (seconds)', 'Average Elapsed Time (seconds)', 'Elapsed Time (seconds)', 'elapsed_time.png')

# Plot Maximum Memory
create_line_plot(td_df, td2_df, 'Maximum Memory (kbytes)', 'Average Maximum Memory (kbytes)', 'Maximum Memory (kbytes)', 'max_memory.png')