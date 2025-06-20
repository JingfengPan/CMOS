import numpy as np
import pickle
import matplotlib.pyplot as plt

def get_valid_input(prompt, valid_values=None, value_type=str):
    while True:
        try:
            value = value_type(input(prompt).strip())
            if valid_values is None or value in valid_values:
                return value
            print(f"Invalid input. Please choose from: {sorted(valid_values)}")
        except ValueError:
            print(f"Please enter a valid {value_type.__name__}")

tensor = np.load('tensor.npy', allow_pickle=True)
with open('mappings.pkl', 'rb') as f:
    mappings = pickle.load(f)

valid_features = set(mappings['feature_to_idx'].keys())
valid_source_nums = sorted(set(k[0] for k in valid_features))
valid_gene_nums = sorted(set(k[1] for k in valid_features))

source_to_genes = {}
for source, gene in valid_features:
    if source not in source_to_genes:
        source_to_genes[source] = []
    source_to_genes[source].append(gene)

print("\nAvailable options:")
print("Valid source-gene combinations:")
for source in sorted(source_to_genes.keys()):
    print(f"Source {source}: genes {sorted(source_to_genes[source])}")
print("\nModalities: Promoter, Protein")
print(f"\nCells: {sorted(mappings['cell_to_idx'].keys())}")
print(f"\nSamples: {sorted(mappings['sample_to_idx'].keys())}\n")

while True:
    source_num = get_valid_input("Enter source number: ", 
                               valid_values=valid_source_nums,
                               value_type=int)
    gene_num = get_valid_input("Enter gene number: ", 
                             valid_values=source_to_genes[source_num],
                             value_type=int)
    modality = get_valid_input("Enter modality (Promoter/Protein): ",
                             valid_values=['Promoter', 'Protein'])
    cell = get_valid_input("Enter cell name: ",
                         valid_values=mappings['cell_to_idx'].keys())
    sample = get_valid_input("Enter sample number: ", 
                           valid_values=mappings['sample_to_idx'].keys(),
                           value_type=int)

    feature_idx = mappings['feature_to_idx'][(source_num, gene_num)]
    modality_idx = mappings['modality_to_idx'][modality]
    cell_idx = mappings['cell_to_idx'][cell]
    sample_idx = mappings['sample_to_idx'][sample]

    rates = [tensor[feature_idx, modality_idx, cell_idx, time_idx, sample_idx]
             for time_idx in mappings['time_to_idx'].values()]
    if not any(not np.isnan(rate) for rate in rates):
        print("\nWarning: This combination does not exist in the tensor.")
        retry = input("Would you like to try a different combination? (y/n): ").lower()
        if retry != 'y':
            break
        continue

    time_points = []
    expression_rates = []
    for time, time_idx in mappings['time_to_idx'].items():
        rate = tensor[feature_idx, modality_idx, cell_idx, time_idx, sample_idx]
        if not np.isnan(rate):
            time_points.append(time)
            expression_rates.append(rate)

    sorted_indices = sorted(range(len(time_points)), key=lambda i: time_points[i])
    sorted_times = [time_points[i] for i in sorted_indices]
    sorted_rates = [expression_rates[i] for i in sorted_indices]

    plt.figure(figsize=(10, 6))
    plt.plot(sorted_times, sorted_rates, marker='.')
    plt.xlabel('Time')
    plt.ylabel('Expression Rate')
    plt.title(f'Gene ({source_num}, {gene_num}), modality {modality} in cell {cell}, sample {sample}')
    
    y_min = min(sorted_rates)
    y_max = max(sorted_rates)
    
    plot_filename = f"plots/gene_{source_num}_{gene_num}_{modality}_cell_{cell}_sample_{sample}.png"
    plt.savefig(plot_filename, bbox_inches='tight')
    plt.close()
    
    print(f"\nPlot has been saved to: {plot_filename}")
    print(f"Number of time points: {len(sorted_times)}")
    print(f"Expression rate range: {y_min} to {y_max}")
    
    another = input("\nWould you like to create another plot? (y/n): ").lower()
    if another != 'y':
        break