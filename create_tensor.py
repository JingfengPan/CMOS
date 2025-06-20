import os
import pandas as pd
import numpy as np
import re
import pickle
import sys

def parse_filename(filename):
    match = re.match(r'WorkSpace_(\d+)_(\d+)_(\d+)_(\d+)\.csv', filename)
    if match:
        source_num, gene_num, construct_num, sample_num = map(int, match.groups())
        # Modality: 'Promoter' or 'Protein'
        if construct_num == 1:
            modality = 'Promoter'
        elif construct_num == 2:
            modality = 'Protein'
        return source_num, gene_num, modality, sample_num
    return None

def create_tensor(data_dir, csv_files):
    features = set()
    cells = set()
    times = set()
    samples = set()

    for filename in csv_files:
        parsed = parse_filename(filename)
        if parsed is None:
            continue
        source_num, gene_num, modality, sample_num = parsed
        df = pd.read_csv(os.path.join(data_dir, filename))
        df = df[df['Table4'].notna()]  # Filter out rows where Table4 is NaN
        features.add((source_num, gene_num))
        cells.update(df['Table1'].unique())
        times.update(df['Table2'].unique())
        samples.add(sample_num)

    feature_to_idx = {f: i for i, f in enumerate(sorted(features))}
    cell_to_idx = {c: i for i, c in enumerate(sorted(cells))}
    time_to_idx = {t: i for i, t in enumerate(sorted(times))}
    sample_to_idx = {s: i for i, s in enumerate(sorted(samples))}
    modality_to_idx = {'Promoter': 0, 'Protein': 1}

    tensor = np.full((
        len(features),
        2,
        len(cells),
        len(times),
        len(samples)
    ), np.nan)

    for filename in csv_files:
        parsed = parse_filename(filename)
        if parsed is None:
            continue
        source_num, gene_num, modality, sample_num = parsed
        modality_idx = modality_to_idx[modality]
        sample_idx = sample_to_idx[sample_num]
        df = pd.read_csv(os.path.join(data_dir, filename))
        df = df[df['Table4'].notna()]  # Filter out rows where Table4 is NaN
        feature_idx = feature_to_idx[(source_num, gene_num)]
        for _, row in df.iterrows():
            cell = row['Table1']
            time = row['Table2']
            expression_rate = row['Table4']
            cell_idx = cell_to_idx[cell]
            time_idx = time_to_idx[time]
            tensor[feature_idx, modality_idx, cell_idx, time_idx, sample_idx] = expression_rate
    return tensor, {
        'feature_to_idx': feature_to_idx,
        'modality_to_idx': modality_to_idx,
        'cell_to_idx': cell_to_idx,
        'time_to_idx': time_to_idx,
        'sample_to_idx': sample_to_idx
    }

def main():
    data_dir = 'data'
    if len(sys.argv) > 1:
        csv_files = [sys.argv[1]]
    else:
        csv_files = [f for f in os.listdir(data_dir) if f.endswith('.csv')]
    print("Creating tensor...")
    tensor, mappings = create_tensor(data_dir, csv_files)
    print("Saving tensor and mappings...")
    np.save('tensor.npy', tensor, allow_pickle=True)
    with open('mappings.pkl', 'wb') as f:
        pickle.dump(mappings, f)
    print("Done!")
    print(f"Tensor shape: {tensor.shape}")
    print("\nMappings:")
    for key, value in mappings.items():
        print(f"{key}: {len(value)} unique values")

if __name__ == "__main__":
    main() 