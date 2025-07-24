import pandas as pd
import numpy as np
import pickle
import json
from collections import defaultdict

# Load tensor and mappings
print("Loading tensor and mappings...")
tensor = np.load('tensor/tensor.npy', allow_pickle=True)
with open('tensor/mappings.pkl', 'rb') as f:
    mappings = pickle.load(f)

# Load cell ID to name mapping
def load_name_dict():
    df = pd.read_csv('data/additional/name_dictionary.csv', header=None)[1:]
    return {str(row[0]): str(row[1]) for _, row in df.iterrows()}

# Load lineage tree data
def load_lineage_trees():
    # Load parent relationships
    parent_df = pd.read_csv('data/additional/lineage_tree_parent.csv')
    parent_dict = {str(row['child']): str(row['parent']) for _, row in parent_df.iterrows()}
    
    return parent_dict

# Load cell fate data
def load_cell_fate_data():
    df = pd.read_csv('data/Cell Fate.csv')
    cell_fate_dict = {}
    for _, row in df.iterrows():
        cell_name = str(row['Cell identity']).strip("'")  # Remove single quotes
        cell_lineage = str(row['Cell lineage']).strip("'")  # Remove single quotes
        cell_fate = str(row['Cell fate']).strip("'")  # Remove single quotes
        cell_fate_dict[cell_name] = {
            'cell_lineage': cell_lineage,
            'cell_fate': cell_fate
        }
    return cell_fate_dict

name_dict = load_name_dict()
parent_dict = load_lineage_trees()
cell_fate_dict = load_cell_fate_data()

for sample_num in range(1, 9):
    print(f"\nProcessing sample {sample_num}...")
    # 1. Load alive time points for each cell
    lifecycle_path = f'data/additional/WT_Sample{sample_num}/WT_Sample{sample_num}_lifescycle.csv'
    cell_lifecycles = {}
    with open(lifecycle_path, 'r') as f:
        for line in f:
            # Remove whitespace and split by comma
            parts = line.split(',')
            if not parts or not parts[0]:
                continue
            cell_id = str(parts[0])
            cell_name = name_dict.get(cell_id, cell_id)
            alive_times = [int(x) for x in parts[1:] if x]
            if alive_times:
                cell_lifecycles[cell_name] = list(range(alive_times[0], alive_times[-1] + 1))

    # 2. Load surface and volume
    surface_path = f'data/additional/WT_Sample{sample_num}/WT_Sample{sample_num}_surface.csv'
    volume_path = f'data/additional/WT_Sample{sample_num}/WT_Sample{sample_num}_volume.csv'
    surface_df = pd.read_csv(surface_path, index_col=0)
    volume_df = pd.read_csv(volume_path, index_col=0)

    # 3. Load contact area (Stat)
    stat_path = f'data/additional/WT_Sample{sample_num}/WT_Sample{sample_num}_Stat.csv'
    stat_df = pd.read_csv(stat_path)
    stat_time_points = [str(tp) for tp in stat_df.columns[2:]]

    # Preprocess Stat for fast lookup
    cell_time_to_neighbors = defaultdict(dict)
    for _, row in stat_df.iterrows():
        c1, c2 = str(row['cell1']), str(row['cell2'])
        for t_str in stat_time_points:
            val = row.get(t_str, None)
            if pd.notna(val) and val > 0:
                cell_time_to_neighbors[(c1, t_str)][c2] = float(val)
                cell_time_to_neighbors[(c2, t_str)][c1] = float(val)

    # 4. Build the output structure
    output = {}
    for cell, times in cell_lifecycles.items():
        output[cell] = {}
        
        # Get parent for this cell
        parent = parent_dict.get(cell, None)
        
        for i, t in enumerate(times):
            t_str = str(t)
            
            # Calculate age (first time point is age 0)
            age = i
            
            # Surface and volume
            surface = None
            volume = None
            if cell in surface_df.columns and t in surface_df.index:
                val = surface_df.at[t, cell]
                surface = float(val) if pd.notna(val) else None
            if cell in volume_df.columns and t in volume_df.index:
                val = volume_df.at[t, cell]
                volume = float(val) if pd.notna(val) else None

            # Neighbours and contacting area
            contacting_area = cell_time_to_neighbors.get((cell, t_str), {})
            neighbours = list(contacting_area.keys())



            # Gene expression
            proteins = {}
            promoters = {}
            for gene_name, feature_idx in mappings['feature_to_idx'].items():
                for modality, modality_idx in mappings['modality_to_idx'].items():
                    try:
                        sample_idx = mappings['sample_to_idx'][sample_num]
                        time_idx = mappings['time_to_idx'][t]
                        cell_idx = mappings['cell_to_idx'][cell]
                    except KeyError:
                        continue
                    rate = tensor[sample_idx, time_idx, cell_idx, modality_idx, feature_idx]
                    if not np.isnan(rate):
                        if modality == 'Protein':
                            proteins[gene_name] = float(rate)
                        elif modality == 'Promoter':
                            promoters[gene_name] = float(rate)

            # Get cell lineage and fate information
            cell_fate_info = cell_fate_dict.get(cell, {})
            cell_lineage = cell_fate_info.get('cell_lineage', None)
            cell_fate = cell_fate_info.get('cell_fate', None)
            
            output[cell][t_str] = {
                "lifecycle": "alive",
                "age": age,
                "parent": parent,
                "cell_lineage": cell_lineage,
                "cell_fate": cell_fate,
                "proteins": proteins,
                "promoters": promoters,
                "surface_area": surface,
                "volume": volume,
                "neighbours": neighbours,
                "contacting_area": contacting_area
            }

    # 5. Save to file
    out_path = f'json/sample_{sample_num}_alive.json'
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=4)
    print(f"Saved {out_path}")