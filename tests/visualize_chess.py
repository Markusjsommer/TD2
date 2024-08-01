import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def read_benchmark_data(directories):
    data = []
    for directory in directories:
        file_path = os.path.join('results/chess_benchmark', directory, 'time_stats.txt')
        with open(file_path, 'r') as file:
            stats = file.readlines()[1].split(', ') # Skip the first line
            user_time = float(stats[0])
            system_time = float(stats[1])
            elapsed_time = float(stats[2])
            max_memory = float(stats[3])
            benchmark_name = directory.split('/')[-1]  # Assuming directory name is the benchmark name
            data.append([benchmark_name, user_time, system_time, elapsed_time, max_memory])
    return data

def save_to_tsv(data, output_file):
    df = pd.DataFrame(data, columns=["Benchmark", "User Time (seconds)", "System Time (seconds)", "Elapsed Time (seconds)", "Maximum Memory (kbytes)"])
    # Convert time columns to minutes
    time_columns = ["User Time (seconds)", "System Time (seconds)", "Elapsed Time (seconds)"]
    for column in time_columns:
        df[column] = (df[column] / 60).round(3)
        # Rename the columns to reflect the new units
        df.rename(columns={column: column.replace("(seconds)", "(minutes)")}, inplace=True)
    # Convert memory columns to MB
    memory_columns = ["Maximum Memory (kbytes)"]
    for column in memory_columns:
        df[column] = (df[column] / 1024).round(3)
        # Rename the columns to reflect the new units
        df.rename(columns={column: column.replace("(kbytes)", "(MB)")}, inplace=True)
    df.to_csv(output_file, sep='\t', index=False)
    return df

def visualize_data(df):
    df_melted = df.melt(id_vars=["Benchmark"], var_name="Metric", value_name="Value")
    plt.figure(figsize=(12, 8))
    sns.barplot(x="Benchmark", y="Value", hue="Metric", data=df_melted)
    plt.title('Benchmark Time and Memory Statistics')
    plt.xticks(rotation=45)
    plt.ylabel('Value')
    plt.legend(loc='upper right')
    plt.savefig('tests/results/consolidated_time_stats.png')
    plt.show()

# List of directories containing the time_stats.txt files
directories = [
    'td',
    'td2_multi1',
    'td2_multi16',
    'td2_multi32',
    'td2_multi64',
    'td2_multi16_mem2',
    'td2_multi16_mem5',
]

# Read data, save to TSV, and visualize
data = read_benchmark_data(directories)
df = save_to_tsv(data, 'tests/results/consolidated_time_stats.tsv')
visualize_data(df)
