import pandas as pd
from py2cytoscape import cyrest
from urllib import parse
import os

# Create cyrest client
cytoscape = cyrest.cyclient()

# How to select the %id column in the edge table
id = parse.quote("%id")

def get_table():

    """ Returns the edge table from Cytoscape with list of pairs of nodes
        and percent identity between the given nodes """

    id = parse.quote("%id")
    # Define columns to retrieve from Cytoscape
    column_list = ['shared name', id]
    # Retrieve edge table
    e_table = cytoscape.table.getTable(table='edge', columns=column_list)

    return e_table


def get_root_dir():

    for file in os.listdir(os.path.join(os.getcwd(), 'network/')):
        if not file.endswith('xgmml'):
            continue
        else:
            print(file)
            dir = os.path.join(os.getcwd(), file+'/')
            print(dir)

            return dir


def main():

    print('Preparing the edge table from Cytoscape for ClusterONE.\n')

    root_dir = get_root_dir()

    table = get_table()

    # Split the 'shared name' column into two columns
    cone_edge_table = pd.concat([table['shared name'].str.split(',', expand=True), table['%id']], axis=1)
    # id has to be between 0 and 1 for ClusterONE tool
    cone_edge_table['%id'] = cone_edge_table['%id'].divide(100)

    # Save the prepared table to csv
    ready_edge_table = cone_edge_table.to_csv(os.path.join(os.getcwd(), root_dir, 'data', 'ready_edge_table.csv'),
                                              header=False, index=False, sep=' ')
    print(f'The prepared edge table will be written to {root_dir}/data/ready_edge_table.csv\n')

if __name__ == '__main__':
    main()

