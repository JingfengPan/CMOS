import matplotlib.pyplot as plt
import json
import os

def get_valid_input(prompt, valid_values=None, value_type=str):
    while True:
        try:
            value = value_type(input(prompt).strip())
            if valid_values is None or value in valid_values:
                return value
            print(f"Invalid input. Please choose from: {sorted(valid_values)}")
        except ValueError:
            print(f"Please enter a valid {value_type.__name__}")

def gather_all_cells_and_genes(json_dir, n_samples=8):
    all_cells = set()
    all_proteins = set()
    all_promoters = set()
    for sample_num in range(1, n_samples+1):
        json_path = os.path.join(json_dir, f"sample_{sample_num}_alive.json")
        if not os.path.exists(json_path):
            continue
        with open(json_path, 'r') as f:
            data = json.load(f)
        all_cells.update(data.keys())
        # Get genes from the first available time point for each cell
        for cell in data:
            for t in data[cell]:
                entry = data[cell][t]
                if 'proteins' in entry:
                    all_proteins.update(entry['proteins'].keys())
                if 'promoters' in entry:
                    all_promoters.update(entry['promoters'].keys())
                break
    return sorted(all_cells), sorted(all_proteins), sorted(all_promoters)

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

def plot_across_samples(cell, modality, gene_name=None, json_dir='json', n_samples=8):
    modality_folder_map = {
        'Proteins Gene Expression Rate': 'proteins',
        'Promoters Gene Expression Rate': 'promoters',
        'Surface Area': 'surface_area',
        'Volume': 'volume',
        'Contacting Area with Neighbours': 'contacting_area',
    }
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
        y_label = ''
        plot_times = []
        plot_values = []
        folder = os.path.join('plots', modality_folder_map[modality])
        os.makedirs(folder, exist_ok=True)
        if modality == 'Proteins Gene Expression Rate':
            y_label = f"Protein {gene_name} Expression Rate"
            for t in time_points:
                entry = data[cell][t]
                rate = entry['proteins'].get(gene_name, None)
                if rate is not None:
                    plot_times.append(int(t))
                    plot_values.append(rate)
            if not plot_times:
                print(f"Sample {sample_num}: No data to plot for cell '{cell}' and modality '{modality}'.")
                continue
            plt.figure(figsize=(10, 6))
            plt.plot(plot_times, plot_values, marker='.')
            plt.xlabel('Time')
            plt.ylabel(y_label)
            title = f"{y_label} in cell {cell}, sample {sample_num}"
            if gene_name:
                title = f"{y_label} ({gene_name}) in cell {cell}, sample {sample_num}"
            plt.title(title)
            plot_filename = f"{modality_folder_map[modality]}_cell_{cell}_sample_{sample_num}"
            if gene_name:
                plot_filename += f"_{gene_name}"
            plot_filename += ".png"
            plot_path = os.path.join(folder, plot_filename)
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()
            print(f"Sample {sample_num}: Plot saved to {plot_path}")
        elif modality == 'Promoters Gene Expression Rate':
            y_label = f"Promoter {gene_name} Expression Rate"
            for t in time_points:
                entry = data[cell][t]
                rate = entry['promoters'].get(gene_name, None)
                if rate is not None:
                    plot_times.append(int(t))
                    plot_values.append(rate)
            if not plot_times:
                print(f"Sample {sample_num}: No data to plot for cell '{cell}' and modality '{modality}'.")
                continue
            plt.figure(figsize=(10, 6))
            plt.plot(plot_times, plot_values, marker='.')
            plt.xlabel('Time')
            plt.ylabel(y_label)
            title = f"{y_label} in cell {cell}, sample {sample_num}"
            if gene_name:
                title = f"{y_label} ({gene_name}) in cell {cell}, sample {sample_num}"
            plt.title(title)
            plot_filename = f"{modality_folder_map[modality]}_cell_{cell}_sample_{sample_num}"
            if gene_name:
                plot_filename += f"_{gene_name}"
            plot_filename += ".png"
            plot_path = os.path.join(folder, plot_filename)
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()
            print(f"Sample {sample_num}: Plot saved to {plot_path}")
        elif modality == 'Surface Area':
            y_label = "Surface Area"
            for t in time_points:
                entry = data[cell][t]
                val = entry.get('surface_area', None)
                if val is not None:
                    plot_times.append(int(t))
                    plot_values.append(val)
            if not plot_times:
                print(f"Sample {sample_num}: No data to plot for cell '{cell}' and modality '{modality}'.")
                continue
            plt.figure(figsize=(10, 6))
            plt.plot(plot_times, plot_values, marker='.')
            plt.xlabel('Time')
            plt.ylabel(y_label)
            title = f"{y_label} in cell {cell}, sample {sample_num}"
            plt.title(title)
            plot_filename = f"{modality_folder_map[modality]}_cell_{cell}_sample_{sample_num}.png"
            plot_path = os.path.join(folder, plot_filename)
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()
            print(f"Sample {sample_num}: Plot saved to {plot_path}")
        elif modality == 'Volume':
            y_label = "Volume"
            for t in time_points:
                entry = data[cell][t]
                val = entry.get('volume', None)
                if val is not None:
                    plot_times.append(int(t))
                    plot_values.append(val)
            if not plot_times:
                print(f"Sample {sample_num}: No data to plot for cell '{cell}' and modality '{modality}'.")
                continue
            plt.figure(figsize=(10, 6))
            plt.plot(plot_times, plot_values, marker='.')
            plt.xlabel('Time')
            plt.ylabel(y_label)
            title = f"{y_label} in cell {cell}, sample {sample_num}"
            plt.title(title)
            plot_filename = f"{modality_folder_map[modality]}_cell_{cell}_sample_{sample_num}.png"
            plot_path = os.path.join(folder, plot_filename)
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()
            print(f"Sample {sample_num}: Plot saved to {plot_path}")
        elif modality == 'Contacting Area with Neighbours':
            # For each time point, plot contacting area with each neighbour as a separate line
            y_label = "Contacting Area with Neighbours"
            all_neighbours = set()
            for t in time_points:
                entry = data[cell][t]
                contacting_area = entry.get('contacting_area', {})
                all_neighbours.update(contacting_area.keys())
            if not all_neighbours:
                print(f"Sample {sample_num}: No contacting area data for cell '{cell}'.")
                continue
            # For each neighbour, collect their contacting area over time
            neighbour_to_times = {n: [] for n in all_neighbours}
            neighbour_to_values = {n: [] for n in all_neighbours}
            for t in time_points:
                entry = data[cell][t]
                contacting_area = entry.get('contacting_area', {})
                for n in all_neighbours:
                    val = contacting_area.get(n, None)
                    if val is not None:
                        neighbour_to_times[n].append(int(t))
                        neighbour_to_values[n].append(val)
            plt.figure(figsize=(14, 6))
            for n in sorted(all_neighbours):
                if neighbour_to_times[n]:
                    plt.plot(neighbour_to_times[n], neighbour_to_values[n], marker='.', label=n)
            plt.xlabel('Time')
            plt.ylabel(y_label)
            title = f"{y_label} in cell {cell}, sample {sample_num}"
            plt.title(title)

            if any(neighbour_to_times[n] for n in all_neighbours):
                ncol = min(4, max(1, len(all_neighbours)//10 + 1))
                plt.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize='x-small', ncol=ncol, borderaxespad=0.)
            plt.tight_layout(rect=[0, 0, 0.8, 1])
            plot_filename = f"{modality_folder_map[modality]}_cell_{cell}_sample_{sample_num}.png"
            plot_path = os.path.join(folder, plot_filename)
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()
            print(f"Sample {sample_num}: Plot saved to {plot_path}")
        else:
            print(f"Unknown modality: {modality}")
            continue

def main():
    json_dir = 'json'
    all_cells, all_proteins, all_promoters = gather_all_cells_and_genes(json_dir)
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
        plot_across_samples(cell, modalities[shortcuts.index(modality_input)], gene_name, json_dir=json_dir)
        another = input("\nWould you like to create another plot? (y/n): ").strip().lower()
        if another != 'y':
            break

if __name__ == '__main__':
    main()