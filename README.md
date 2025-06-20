1. Download all data into the "data" folder.
2. "create_tensor.py" to 
- build a 5D tensor from all the data,
- save it as "tensor.npy" and,
- save the label‐to‐index mappings (feature, modality, cell, time, sample) as "mappings.pkl".
3. "read_tensor_shape.py" to print out the tensor's shape, which is (614, 2, 974, 195, 8).
4. "plot_tensor.py" to visualize selected tensor samples and save the resulting figures in the "plots" folder.
5."tensor_to_json.py" to convert each tensor sample into JSON and output it as "gene_expressions.json".
6. Since "gene_expressions.json" is too large (3.21GB) and cannot be opened directly, write "read_json.py" using "ijson" library to stream and inspect its contents.


