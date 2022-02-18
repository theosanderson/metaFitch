import treeswift
import pandas as pd

import argparse
import tqdm
import xopen

parser = argparse.ArgumentParser(description='NumDescendants')
parser.add_argument(
    '-t',
    '--tree',
    #file
    required=True,
    help='Newick tree')
parser.add_argument('-o',
                    '--output',
                    type=str,
                    required=True,
                    help='Output file')

args = parser.parse_args()

tree_file = xopen.xopen(args.tree, 'rt')
tree_string = tree_file.read()
tree_file.close()
print('Reading tree...')
tree = treeswift.read_tree_newick(tree_string)
# get number of descendants
for node in tqdm.tqdm(tree.traverse_postorder()):
    if node.is_leaf():
        node.num_descendants = 0
    else:
        node.num_descendants = len(node.children) + sum(
            [c.num_descendants for c in node.children])

label_to_node = tree.label_to_node(selection='all')
keys = list(label_to_node.keys())
num_descendants = [label_to_node[k].num_descendants for k in keys]
edge_lengths = [label_to_node[k].edge_length for k in keys]
result = pd.DataFrame({
    'strain': keys,
    'num_descendants': num_descendants,
    'edge_length': edge_lengths
})
result.to_csv(args.output, sep='\t', index=False)
