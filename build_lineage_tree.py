import pandas as pd

# Load all cell names from name_dictionary.csv
name_dict_path = 'data/additional/name_dictionary.csv'
name_df = pd.read_csv(name_dict_path, header=None, skiprows=1)
all_cells = set(name_df[1].values)

# Load the initial lineage tree from lineage_tree_children_beginning.csv
children_path = 'data/additional/lineage_tree_children_beginning.csv'
children_df = pd.read_csv(children_path)

# Build initial parent->children mapping from the provided file
parent_to_children = {}
for _, row in children_df.iterrows():
    parent = row['parent']
    c1 = row['child1'] if pd.notna(row['child1']) else ''
    c2 = row['child2'] if pd.notna(row['child2']) else ''
    if parent:
        parent_to_children[parent] = [c1, c2]

# Build a set of all parents that have explicit children listed
explicit_parents = set(parent_to_children.keys())
# Build a set of all children that are explicitly listed
explicit_children = set()
for v in parent_to_children.values():
    explicit_children.update([x for x in v if x])

# For all cells in name_dictionary.csv, if not already in parent_to_children and not a root, infer children
for cell in all_cells:
    if cell not in parent_to_children:
        # Only add children if this cell is not a leaf (i.e., if its children exist in all_cells)
        children_pairs = []
        # Check for 'a' and 'p' pair
        child_a = cell + 'a'
        child_p = cell + 'p'
        if child_a in all_cells and child_p in all_cells:
            children_pairs = [child_a, child_p]
        # Check for 'l' and 'r' pair
        child_l = cell + 'l'
        child_r = cell + 'r'
        if child_l in all_cells and child_r in all_cells:
            children_pairs = [child_l, child_r]
        # Check for 'd' and 'v' pair
        child_d = cell + 'd'
        child_v = cell + 'v'
        if child_d in all_cells and child_v in all_cells:
            children_pairs = [child_d, child_v]
        if children_pairs:
            parent_to_children[cell] = children_pairs

# Build the full lineage_tree_children.csv
full_children_rows = []
for parent, children in parent_to_children.items():
    c1 = children[0] if len(children) > 0 else ''
    c2 = children[1] if len(children) > 1 else ''
    full_children_rows.append({'parent': parent, 'child1': c1, 'child2': c2})
full_children_df = pd.DataFrame(full_children_rows)
full_children_df = full_children_df.sort_values('parent')
full_children_df.to_csv('data/additional/lineage_tree_children.csv', index=False)

# Build the lineage_tree_parent.csv (child, parent)
child_to_parent = {}
for parent, children in parent_to_children.items():
    for child in children:
        if child:
            child_to_parent[child] = parent

parent_rows = []
for cell in sorted(all_cells):
    parent = child_to_parent.get(cell, '')
    parent_rows.append({'child': cell, 'parent': parent})
parent_df = pd.DataFrame(parent_rows)
parent_df.to_csv('data/additional/lineage_tree_parent.csv', index=False)

# Print statistics
print('Lineage tree children written to data/additional/lineage_tree_children.csv')
print('Lineage tree parent written to data/additional/lineage_tree_parent.csv')
print(f"Total distinct cells in name_dictionary.csv: {len(all_cells)}")
print(f"Total distinct parents in lineage_tree_children.csv: {full_children_df['parent'].nunique()}")
print(f"Total distinct children in lineage_tree_children.csv: {pd.unique(full_children_df[['child1','child2']].values.ravel('K')).size}")
print(f"Total distinct children in lineage_tree_parent.csv: {parent_df['child'].nunique()}")
print(f"Total distinct parents in lineage_tree_parent.csv: {parent_df['parent'].nunique()}")

# Check for missing parents
missing_parents = [cell for cell in all_cells if cell not in child_to_parent]
if missing_parents:
    print(f"Cells with no parent: {missing_parents}")
