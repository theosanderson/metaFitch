import treeswift
import pandas as pd

import argparse
import tqdm
import xopen

parser = argparse.ArgumentParser(description='Metafitch')
parser.add_argument(
    '-t',
    '--tree',
    #file
    required=True,
    help='Newick tree')
parser.add_argument('-m', '--metadata', required=True, help='Metadata file')
parser.add_argument('-o',
                    '--output',
                    type=str,
                    required=True,
                    help='Output file')
parser.add_argument('-f',
                    '--fields',
                    type=str,
                    required=False,
                    help='Comma-separated to be used')

parser.add_argument('-q',
                    '--id_field',
                    type=str,
                    required=False,
                    help='Field to use as ID')

args = parser.parse_args()

tree_file = xopen.xopen(args.tree, 'rt')
tree_string = tree_file.read()
tree_file.close()
print('Reading tree...')
tree = treeswift.read_tree_newick(tree_string)
del tree_string
metadata = pd.read_csv(args.metadata,
                       sep='\t' if args.metadata.endswith('.tsv')
                       or args.metadata.endswith('.txt')
                       or args.metadata.endswith('.tsv.gz') else ',')
print("Metadata began:\n", metadata.head())

# if no id field is specified, use the first column
if args.id_field is None:
    args.id_field = metadata.columns[0]

print(f"Using {args.id_field} as ID")

# If no fields specified use all columns except the id field
if args.fields is None:
    fields = [c for c in metadata.columns if c != args.id_field]
else:
    fields = args.fields.split(',')

print(f"Using {', '.join(fields)} as fields")

metadata.set_index(args.id_field, inplace=True)
label_to_node = tree.label_to_node()
results = {}
results['strain'] = label_to_node.keys()
for field in tqdm.tqdm(fields):
    for node in tqdm.tqdm(tree.traverse_preorder()):
        if node.label is not None and node.label in metadata.index:
            node.character = set([metadata.loc[node.label, field]])

    # Fitch's algorithm step 1
    for node in tqdm.tqdm(tree.traverse_postorder()):
        if node.is_leaf():
            # check it has an attribute character
            if not hasattr(node, 'character'):
                node.character = set()

        else:
            intersection = node.children[0].character
            union = node.children[0].character
            for child in node.children[1:]:
                intersection = intersection & child.character
                union = union | child.character
            if len(intersection) == 0:
                node.character = union
            else:
                node.character = intersection

    # Fitch's algorithm step 2
    for node in tqdm.tqdm(tree.traverse_preorder()):
        if node.is_leaf():
            continue
        else:
            parent = node.parent
            if parent is None:
                continue
            intersection = parent.character & node.character
            if len(intersection) == 0:
                print(
                    f"No intersection between {parent.label} and {node.label}, with characters {parent.character} and {node.character}"
                )

            else:
                node.character = intersection
    results[field] = []
    for label in label_to_node:
        if len(label_to_node[label].character) == 1:
            results[field].append(list(label_to_node[label].character)[0])
        else:
            results[field].append('NA')

# Write results to file
output = pd.DataFrame(results)
output.to_csv(
    args.output,
    sep='\t' if args.output.endswith('.tsv') or args.output.endswith('.txt')
    or args.output.endswith('.tsv.gz') else ',',
    index=False)
