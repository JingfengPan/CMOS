import os
import pandas as pd
import numpy as np
import re
import pickle

def parse_filename(filename):
    match = re.match(r'WorkSpace_(\d+)_(\d+)_(\d+)_(\d+)\.csv', filename)
    if match:
        _, _, construct_num, sample_num = map(int, match.groups())
        # Modality: 'Promoter' or 'Protein'
        if construct_num == 1:
            modality = 'Promoter'
        elif construct_num == 2:
            modality = 'Protein'
        return modality, sample_num
    return None

def build_filename_to_gene_map(fileinfo_path):
    filename_to_gene = {}
    with open(fileinfo_path, 'r', encoding='utf-8') as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            filename = parts[0]
            gene = parts[1]
            filename_to_gene[filename] = gene
    return filename_to_gene

def create_tensor(data_dir, csv_files):
    features = set()
    cells = set()
    times = set()
    samples = set()

    fileinfo_path = 'data/raw/FileInfo.txt'
    filename_to_gene = build_filename_to_gene_map(fileinfo_path)

    for filename in csv_files:
        parsed = parse_filename(filename)
        if parsed is None:
            continue
        modality, sample_num = parsed
        gene_name = filename_to_gene.get(filename)
        if gene_name is None:
            continue
        df = pd.read_csv(os.path.join(data_dir, filename))
        # df = df[df['Table4'].notna()]  # Filter out rows where Table4 is NaN
        features.add(gene_name)
        cells.update(df['Table1'].unique())
        times.update(df['Table2'].unique())
        samples.add(sample_num)

    feature_to_idx = {f: i for i, f in enumerate(sorted(features))}
    cell_to_idx = {c: i for i, c in enumerate(sorted(cells))}
    time_to_idx = {t: i for i, t in enumerate(sorted(times))}
    sample_to_idx = {s: i for i, s in enumerate(sorted(samples))}
    modality_to_idx = {'Promoter': 0, 'Protein': 1}

    # Tensor shape: (sample, time, cell, modality, feature)
    tensor = np.full((
        len(samples),
        len(times),
        len(cells),
        2,
        len(features)
    ), np.nan)

    for filename in csv_files:
        parsed = parse_filename(filename)
        if parsed is None:
            continue
        modality, sample_num = parsed
        gene_name = filename_to_gene.get(filename)
        if gene_name is None:
            continue
        modality_idx = modality_to_idx[modality]
        sample_idx = sample_to_idx[sample_num]
        df = pd.read_csv(os.path.join(data_dir, filename))
        df = df[df['Table4'].notna()]
        feature_idx = feature_to_idx[gene_name]
        for _, row in df.iterrows():
            cell = row['Table1']
            time = row['Table2']
            expression_rate = row['Table4']
            cell_idx = cell_to_idx[cell]
            time_idx = time_to_idx[time]
            tensor[sample_idx, time_idx, cell_idx, modality_idx, feature_idx] = expression_rate
    return tensor, {
        'sample_to_idx': sample_to_idx,
        'time_to_idx': time_to_idx,
        'cell_to_idx': cell_to_idx,
        'modality_to_idx': modality_to_idx,
        'feature_to_idx': feature_to_idx,
        'filename_to_gene': filename_to_gene
    }

def main():
    data_dir = 'data/raw'

    csv_files = [f for f in os.listdir(data_dir) if f.endswith('.csv')]
    print("Creating tensor...")
    tensor, mappings = create_tensor(data_dir, csv_files)
    print("Saving tensor and mappings...")
    np.save('tensor/tensor.npy', tensor, allow_pickle=True)
    with open('tensor/mappings.pkl', 'wb') as f:
        pickle.dump(mappings, f)
    print("Done!")
    print(f"\nTensor shape: {tensor.shape}")
    print("\nMappings:")
    for key, value in mappings.items():
        print(f"{key}: {len(value)} unique values")

if __name__ == "__main__":
    main() 