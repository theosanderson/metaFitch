import treeswift
import pandas as pd

import argparse
from alive_progress import alive_it
import xopen

def get_parser():
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
    parser.add_argument('-a',
                        '--tolerate_ambiguity',
                        action='store_true',
                        help= "Tolerate ambiguity in the reconstructions")
    return parser

def main():
    parser = get_parser()


    args = parser.parse_args()

    tree_file = xopen.xopen(args.tree, 'rt')
    tree_string = tree_file.read()
    tree_file.close()
    print('Reading tree...')
    tree = treeswift.read_tree_newick(tree_string)
    del tree_string
    
    # first just get the top line
    headers = pd.read_csv(args.metadata,
                        sep='\t' if args.metadata.endswith('.tsv')
                        or args.metadata.endswith('.txt')
                        or args.metadata.endswith('.tsv.gz') else ',', 
                        nrows=1)
    if args.fields is not None:
        fields = args.fields.split(',')
    else:
        fields = list(headers.columns)
    if args.id_field is not None:
        id_field = args.id_field
    else:
        id_field = fields[0]

    print('Reading metadata...')

    # Read all columns as strings
    
    metadata = pd.read_csv(args.metadata,
                        sep='\t' if args.metadata.endswith('.tsv')
                        or args.metadata.endswith('.txt')
                        or args.metadata.endswith('.tsv.gz') else ',',
                        usecols=fields)
    

    # check for duplicates in id_field
    if len(metadata[id_field].unique()) != len(metadata):
        print('Duplicates in ID field! Will randomly drop duplicates.')
        metadata = metadata.drop_duplicates(subset=id_field)
    
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
    label_to_node = tree.label_to_node(selection='all')
    results = {}
    results['strain'] = label_to_node.keys()
    for field in alive_it(fields):
        for node in alive_it(tree.traverse_preorder()):
            if node.label is not None and node.label in metadata.index:
                node.character = set([metadata.loc[node.label, field]])

        # Fitch's algorithm step 1
        for node in alive_it(tree.traverse_postorder()):
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
        for node in alive_it(tree.traverse_preorder()):
            if node.is_leaf():
                continue
            else:
                parent = node.parent
                if parent is None:
                    continue
                intersection = parent.character & node.character
                if len(intersection) == 0:
                    pass
                else:
                    node.character = intersection
                if len(node.character) == 0:
                    node.character = parent.character
        results[field] = []
        for label in label_to_node:
            if len(label_to_node[label].character) == 1:
                results[field].append(list(label_to_node[label].character)[0])
            else:
                if args.tolerate_ambiguity:
                    results[field].append(", ".join([str(x) for x in label_to_node[label].character]))
                else:
                    results[field].append("");
                

    # Write results to file
    output = pd.DataFrame(results)
    output.to_csv(
        args.output,
        sep='\t' if args.output.endswith('.tsv') or args.output.endswith('.txt')
        or args.output.endswith('.tsv.gz') else ',',
        index=False)

if __name__ == '__main__':
    main()
