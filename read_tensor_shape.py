import numpy as np

tensor = np.load('tensor.npy', allow_pickle=True)

print(f"Tensor shape: {tensor.shape}")
print("\nDimension meanings:")
print(f"1. Features (source_num, gene_num): {tensor.shape[0]} unique combinations")
print(f"2. Modalities (Promoter/Protein): {tensor.shape[1]} types")
print(f"3. Cells: {tensor.shape[2]} unique cells")
print(f"4. Time points: {tensor.shape[3]} unique time points")
print(f"5. Samples: {tensor.shape[4]} unique samples")
