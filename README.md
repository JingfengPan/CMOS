1. The `create_tensor.py` file creates a tensor from raw data (cells’ protein and promoter gene expression rates over time).
2. The `create_json_alive.py` file combines raw and additional data (cells’ age, parent, surface area, volume, and contacting area with neighbors) and outputs them to JSON files. This version only records the data when a cell is alive.
3. The `create_json_unborn.py` file outputs a comprehensive version of JSON files. It includes all time points of the sample even if the cell is "unborn" or "dead"/"divided". The JSON files also list the two children into which the cell has divided, if any.
4. The `plot_json.py` file visualizes how the five modalities: `surface_area`, `volume`, `contacting_area`, `proteins` and `promoters` vary over time.
5. The `build_lineage_tree.py` file builds two lineage trees: `lineage_tree_parent.csv` and `lineage_tree_children.csv`.
