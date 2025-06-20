import numpy as np
import pickle
import json

tensor = np.load('tensor.npy', allow_pickle=True)
with open('mappings.pkl', 'rb') as f:
    mappings = pickle.load(f)

idx_to_feature = {v: k for k, v in mappings['feature_to_idx'].items()}
idx_to_modality = {v: k for k, v in mappings['modality_to_idx'].items()}
idx_to_cell = {v: k for k, v in mappings['cell_to_idx'].items()}
idx_to_sample = {v: k for k, v in mappings['sample_to_idx'].items()}
idx_to_time = {v: k for k, v in mappings['time_to_idx'].items()}

output = []

for feature_idx in range(tensor.shape[0]):
    source_num, gene_num = idx_to_feature[feature_idx]
    for modality_idx in range(tensor.shape[1]):
        modality = idx_to_modality[modality_idx]
        for cell_idx in range(tensor.shape[2]):
            cell = idx_to_cell[cell_idx]
            for sample_idx in range(tensor.shape[4]):
                sample = idx_to_sample[sample_idx]
                time_points = []
                expression_rates = []
                for time_idx in range(tensor.shape[3]):
                    rate = tensor[feature_idx, modality_idx, cell_idx, time_idx, sample_idx]
                    if not np.isnan(rate):
                        time_points.append(idx_to_time[time_idx])
                        expression_rates.append(float(rate))
                if time_points:
                    output.append({
                        "features": {
                            "source_num": int(source_num),
                            "gene_num": int(gene_num)
                        },
                        "modality": str(modality),
                        "cell": str(cell),
                        "sample": int(sample),
                        "time_points": [int(tp) for tp in time_points],
                        "expression_rates": [float(rate) for rate in expression_rates]
                    })

with open("gene_expressions.json", "w") as f:
    json.dump(output, f, indent=4)

print("All tensor combinations have been saved to gene_expressions.json")