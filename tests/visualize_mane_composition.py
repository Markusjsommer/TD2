import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 

def read_compositions(name, file):
    with open(file, 'r') as f:
        lines = f.readlines()
    counts = {'complete': 0, 'internal': 0, '5prime_partial': 0, '3prime_partial': 0}
    for line in lines:
        if line.startswith('>'):
            orf_type = line.split('type:')[1].split()[0]
            counts[orf_type] += 1
    print(counts)
    total = sum(counts.values())
    print(total)
    with open('results/orf_composition.csv', 'a') as f:
        f.write(f'{name},{counts["complete"]},{counts["3prime_partial"]},{counts["5prime_partial"]},{counts["internal"]},{total}\n')

def generate_graphs():
    df = pd.read_csv('results/orf_composition.csv')
    
    print(df.head())
    
    # # Generate pie charts for each tool
    # for index, row in df.iterrows():
    #     tool_name = row['Tool']
    #     counts = row[1:-1]  # Exclude the Tool and Total columns
    #     labels = counts.index
        
    #     plt.figure(figsize=(6, 6))
    #     plt.pie(counts, labels=labels, autopct='%1.1f%%', startangle=140)
    #     plt.title(f'Distribution of {tool_name}')
    #     plt.savefig(f'results/{tool_name}_composition.png', dpi=300)

    # Generate a bar plot for total counts comparison
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Tool', y='total', data=df, color='grey')
    plt.title('Comparison of Total Counts Across Tools')
    plt.ylabel('Total Count')
    plt.xticks(rotation=45, ha="right")

    # Label the actual count on the bar
    for index, row in df.iterrows():
        tool_name = row['Tool']
        count = row['total']
        plt.text(index, count, str(count), ha='center', va='bottom')
        
    plt.tight_layout()
    plt.savefig('results/total_counts_comparison.png', dpi=300)
    
if __name__ == "__main__":
    generate_graphs()