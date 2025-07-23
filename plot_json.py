import matplotlib.pyplot as plt
import json
import os
import numpy as np

def get_valid_input(prompt, valid_values=None, value_type=str):
    while True:
        try:
            value = value_type(input(prompt).strip())
            if valid_values is None or value in valid_values:
                return value
            print(f"Invalid input. Please choose from: {sorted(valid_values)}")
        except ValueError:
            print(f"Please enter a valid {value_type.__name__}")

def gather_all_cells(json_dir, n_samples=8):
    all_cells = set()
    for sample_num in range(1, n_samples+1):
        json_path = os.path.join(json_dir, f"sample_{sample_num}_alive.json")
        if not os.path.exists(json_path):
            continue
        with open(json_path, 'r') as f:
            data = json.load(f)
        all_cells.update(data.keys())
    return sorted(all_cells)

def get_cell_specific_genes(cell, json_dir, n_samples=8):
    protein_genes = set()
    promoter_genes = set()
    for sample_num in range(1, n_samples+1):
        json_path = os.path.join(json_dir, f"sample_{sample_num}_alive.json")
        if not os.path.exists(json_path):
            continue
        with open(json_path, 'r') as f:
            data = json.load(f)
        if cell not in data:
            continue
        for t in data[cell]:
            entry = data[cell][t]
            if 'proteins' in entry:
                protein_genes.update(entry['proteins'].keys())
            if 'promoters' in entry:
                promoter_genes.update(entry['promoters'].keys())
    return sorted(protein_genes), sorted(promoter_genes)

def plot_combined_across_samples(cell, modality, gene_name=None, json_dir='json', n_samples=8):
    """Plot all samples on the same graph for surface area, volume, proteins, and promoters"""
    modality_folder_map = {
        'Proteins Gene Expression Rate': 'proteins',
        'Promoters Gene Expression Rate': 'promoters',
        'Surface Area': 'surface_area',
        'Volume': 'volume',
        'Contacting Area with Neighbours': 'contacting_area',
    }
    
    folder = os.path.join('plots', modality_folder_map[modality])
    os.makedirs(folder, exist_ok=True)
    
    plt.figure(figsize=(12, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, n_samples))
    
    for sample_num in range(1, n_samples+1):
        json_path = os.path.join(json_dir, f"sample_{sample_num}_alive.json")
        if not os.path.exists(json_path):
            print(f"Sample {sample_num}: JSON file not found, skipping.")
            continue
        with open(json_path, 'r') as f:
            data = json.load(f)
        if cell not in data:
            print(f"Sample {sample_num}: Cell '{cell}' not found, skipping.")
            continue
            
        time_points = sorted(data[cell].keys(), key=lambda x: int(x))
        plot_times = []
        plot_values = []
        
        if modality == 'Proteins Gene Expression Rate':
            for t in time_points:
                entry = data[cell][t]
                rate = entry['proteins'].get(gene_name, None)
                if rate is not None:
                    plot_times.append(int(t))
                    plot_values.append(rate)
        elif modality == 'Promoters Gene Expression Rate':
            for t in time_points:
                entry = data[cell][t]
                rate = entry['promoters'].get(gene_name, None)
                if rate is not None:
                    plot_times.append(int(t))
                    plot_values.append(rate)
        elif modality == 'Surface Area':
            for t in time_points:
                entry = data[cell][t]
                val = entry.get('surface_area', None)
                if val is not None:
                    plot_times.append(int(t))
                    plot_values.append(val)
        elif modality == 'Volume':
            for t in time_points:
                entry = data[cell][t]
                val = entry.get('volume', None)
                if val is not None:
                    plot_times.append(int(t))
                    plot_values.append(val)
        
        if plot_times:
            plt.plot(plot_times, plot_values, marker='.', label=f'Sample {sample_num}', 
                    color=colors[sample_num-1], linewidth=2, markersize=6)
    
    # Set labels and title
    if modality == 'Proteins Gene Expression Rate':
        y_label = f"Protein {gene_name} Expression Rate"
        title = f"{y_label} in cell {cell} across all samples"
    elif modality == 'Promoters Gene Expression Rate':
        y_label = f"Promoter {gene_name} Expression Rate"
        title = f"{y_label} in cell {cell} across all samples"
    elif modality == 'Surface Area':
        y_label = "Surface Area"
        title = f"{y_label} in cell {cell} across all samples"
    elif modality == 'Volume':
        y_label = "Volume"
        title = f"{y_label} in cell {cell} across all samples"
    
    plt.xlabel('Time')
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    plot_filename = f"{modality_folder_map[modality]}_cell_{cell}_combined"
    if gene_name:
        plot_filename += f"_{gene_name}"
    plot_filename += ".png"
    plot_path = os.path.join(folder, plot_filename)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Combined plot saved to {plot_path}")

def plot_contacting_area_group(cell, json_dir='json', n_samples=8):
    """Create a group of plots for contacting area across all samples"""
    folder = os.path.join('plots', 'contacting_area')
    os.makedirs(folder, exist_ok=True)
    
    fig, axes = plt.subplots(2, 4, figsize=(36, 16))
    axes = axes.flatten()
    
    for sample_num in range(1, n_samples+1):
        json_path = os.path.join(json_dir, f"sample_{sample_num}_alive.json")
        if not os.path.exists(json_path):
            print(f"Sample {sample_num}: JSON file not found, skipping.")
            continue
        with open(json_path, 'r') as f:
            data = json.load(f)
        if cell not in data:
            print(f"Sample {sample_num}: Cell '{cell}' not found, skipping.")
            continue
            
        time_points = sorted(data[cell].keys(), key=lambda x: int(x))
        all_neighbours = set()
        
        # Collect all neighbours
        for t in time_points:
            entry = data[cell][t]
            contacting_area = entry.get('contacting_area', {})
            all_neighbours.update(contacting_area.keys())
        
        if not all_neighbours:
            print(f"Sample {sample_num}: No contacting area data for cell '{cell}'.")
            continue
        
        # Plot for this sample
        ax = axes[sample_num-1]
        colors = plt.cm.tab10(np.linspace(0, 1, len(all_neighbours)))
        
        for i, neighbour in enumerate(sorted(all_neighbours)):
            neighbour_times = []
            neighbour_values = []
            for t in time_points:
                entry = data[cell][t]
                contacting_area = entry.get('contacting_area', {})
                val = contacting_area.get(neighbour, None)
                if val is not None:
                    neighbour_times.append(int(t))
                    neighbour_values.append(val)
            
            if neighbour_times:
                ax.plot(neighbour_times, neighbour_values, marker='.', 
                       label=neighbour, color=colors[i], linewidth=1, markersize=3)
        
        ax.set_xlabel('Time')
        ax.set_ylabel('Contacting Area')
        ax.set_title(f'Sample {sample_num}')
        ax.grid(True, alpha=0.3)
        
        # Add legend if there are neighbours - place outside the plot
        if all_neighbours:
            ax.legend(fontsize='x-small', loc='center left', bbox_to_anchor=(1.02, 0.5))
    
    # Hide empty subplots
    for i in range(n_samples, 8):
        axes[i].set_visible(False)
    
    plt.suptitle(f'Contacting Area with Neighbours in cell {cell} across all samples', fontsize=16)
    plt.tight_layout()
    
    # Save plot
    plot_filename = f"contacting_area_cell_{cell}_group.png"
    plot_path = os.path.join(folder, plot_filename)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Contacting area group plot saved to {plot_path}")

def plot_modality_group(cell, modality, gene_name=None, json_dir='json', n_samples=8):
    """Create a group of plots for other modalities across all samples"""
    modality_folder_map = {
        'Proteins Gene Expression Rate': 'proteins',
        'Promoters Gene Expression Rate': 'promoters',
        'Surface Area': 'surface_area',
        'Volume': 'volume',
    }
    
    folder = os.path.join('plots', modality_folder_map[modality])
    os.makedirs(folder, exist_ok=True)
    
    fig, axes = plt.subplots(2, 4, figsize=(32, 12))
    axes = axes.flatten()
    
    for sample_num in range(1, n_samples+1):
        json_path = os.path.join(json_dir, f"sample_{sample_num}_alive.json")
        if not os.path.exists(json_path):
            print(f"Sample {sample_num}: JSON file not found, skipping.")
            continue
        with open(json_path, 'r') as f:
            data = json.load(f)
        if cell not in data:
            print(f"Sample {sample_num}: Cell '{cell}' not found, skipping.")
            continue
            
        time_points = sorted(data[cell].keys(), key=lambda x: int(x))
        plot_times = []
        plot_values = []
        
        # Collect data for this sample
        if modality == 'Proteins Gene Expression Rate':
            for t in time_points:
                entry = data[cell][t]
                rate = entry['proteins'].get(gene_name, None)
                if rate is not None:
                    plot_times.append(int(t))
                    plot_values.append(rate)
        elif modality == 'Promoters Gene Expression Rate':
            for t in time_points:
                entry = data[cell][t]
                rate = entry['promoters'].get(gene_name, None)
                if rate is not None:
                    plot_times.append(int(t))
                    plot_values.append(rate)
        elif modality == 'Surface Area':
            for t in time_points:
                entry = data[cell][t]
                val = entry.get('surface_area', None)
                if val is not None:
                    plot_times.append(int(t))
                    plot_values.append(val)
        elif modality == 'Volume':
            for t in time_points:
                entry = data[cell][t]
                val = entry.get('volume', None)
                if val is not None:
                    plot_times.append(int(t))
                    plot_values.append(val)
        
        # Plot for this sample
        ax = axes[sample_num-1]
        if plot_times:
            ax.plot(plot_times, plot_values, marker='.', linewidth=2, markersize=6)
        
        # Set labels and title
        if modality == 'Proteins Gene Expression Rate':
            y_label = f"Protein {gene_name} Expression Rate"
        elif modality == 'Promoters Gene Expression Rate':
            y_label = f"Promoter {gene_name} Expression Rate"
        elif modality == 'Surface Area':
            y_label = "Surface Area"
        elif modality == 'Volume':
            y_label = "Volume"
        
        ax.set_xlabel('Time')
        ax.set_ylabel(y_label)
        ax.set_title(f'Sample {sample_num}')
        ax.grid(True, alpha=0.3)
    
    # Hide empty subplots
    for i in range(n_samples, 8):
        axes[i].set_visible(False)
    
    # Set overall title
    if modality == 'Proteins Gene Expression Rate':
        title = f"Protein {gene_name} Expression Rate in cell {cell} across all samples"
    elif modality == 'Promoters Gene Expression Rate':
        title = f"Promoter {gene_name} Expression Rate in cell {cell} across all samples"
    elif modality == 'Surface Area':
        title = f"Surface Area in cell {cell} across all samples"
    elif modality == 'Volume':
        title = f"Volume in cell {cell} across all samples"
    
    plt.suptitle(title, fontsize=16)
    plt.tight_layout()
    
    # Save plot
    plot_filename = f"{modality_folder_map[modality]}_cell_{cell}_group"
    if gene_name:
        plot_filename += f"_{gene_name}"
    plot_filename += ".png"
    plot_path = os.path.join(folder, plot_filename)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Group plot saved to {plot_path}")

def main():
    json_dir = 'json'
    all_cells = gather_all_cells(json_dir)
    while True:
        print("\nAvailable cells:")
        print(", ".join(all_cells))
        cell = get_valid_input("Enter cell name: ", valid_values=all_cells)
        
        modality_map = {
            'prot': 'Proteins Gene Expression Rate',
            'prom': 'Promoters Gene Expression Rate',
            'sa': 'Surface Area',
            'v': 'Volume',
            'ca': 'Contacting Area with Neighbours',
        }
        shortcuts = ['prot', 'prom', 'sa', 'v', 'ca']
        modalities = [
            'Proteins Gene Expression Rate',
            'Promoters Gene Expression Rate',
            'Surface Area',
            'Volume',
            'Contacting Area with Neighbours'
        ]
        print("\nAvailable modalities:")
        for m, s in zip(modalities, shortcuts):
            print(f"  {m} ({s})")
        modality_input = get_valid_input("Enter modality (shortcut): ", valid_values=shortcuts)
        modality = modality_map[modality_input]
        
        print("\nPlot style:")
        if modality == 'Contacting Area with Neighbours':
            print("  1. Group plot (subplots for each sample)")
            plot_style = get_valid_input("Enter plot style (1): ", valid_values=['1'])
        else:
            print("  1. Combined plot (all samples on one graph)")
            print("  2. Group plot (subplots for each sample)")
            plot_style = get_valid_input("Enter plot style (1/2): ", valid_values=['1', '2'])
        gene_name = None
        
        # Only show genes present for the selected cell
        cell_proteins, cell_promoters = get_cell_specific_genes(cell, json_dir)
        if modality == 'Proteins Gene Expression Rate':
            print(f"\nAvailable protein genes for {cell}:")
            print(", ".join(cell_proteins))
            gene_name = get_valid_input("Enter protein gene name: ", valid_values=cell_proteins)
        elif modality == 'Promoters Gene Expression Rate':
            print(f"\nAvailable promoter genes for {cell}:")
            print(", ".join(cell_promoters))
            gene_name = get_valid_input("Enter promoter gene name: ", valid_values=cell_promoters)
        
        # Execute the appropriate plotting function
        if modality == 'Contacting Area with Neighbours':
            # Only group plot available for contacting area
            plot_contacting_area_group(cell, json_dir=json_dir)
        else:
            # For other modalities, handle combined and group options
            if plot_style == '1':
                plot_combined_across_samples(cell, modalities[shortcuts.index(modality_input)], gene_name, json_dir=json_dir)
            elif plot_style == '2':
                plot_modality_group(cell, modalities[shortcuts.index(modality_input)], gene_name, json_dir=json_dir)
        
        another = input("\nWould you like to create another plot? (y/n): ").strip().lower()
        if another != 'y':
            break

if __name__ == '__main__':
    main()