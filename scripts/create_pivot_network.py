from py2cytoscape import cyrest
import pandas as pd
import networkx as nx
import os

# Create cyrest client
cytoscape = cyrest.cyclient()

def get_root_dir():
    for file in os.listdir(os.path.join(os.getcwd(), 'network/')):
        if not file.endswith('xgmml'):
            continue
        else:
            dir = os.path.join(os.getcwd(), file + '/')

            return dir


def main():
    root_dir = get_root_dir()

    # Read edge table from cytoscape for full network
    edge_table = cytoscape.table.getTable(columns=['shared name'], table='edge')
    # Read node table
    node_table = pd.read_csv(os.path.join(os.getcwd(), root_dir, 'data', 'node_table.csv'))

    # Create a dataframe with columns source and target columns and split the shared name column into these two
    source_target_list = pd.DataFrame(edge_table['shared name'].str.split(',', 1).to_list(), columns=['source', 'target'])

    # Get only shared name and clusterONE columns from node_table
    reduced_node_table = node_table[['shared name', 'clusterONE']]

    # Add cluster numbers to corresponding source nodes
    df = pd.merge(source_target_list, reduced_node_table, copy=True, left_on=['source'], right_on=['shared name'])
    df = df.rename(columns={'clusterONE': 'clusterONE_source'})
    # Add cluster numbers to corresponding target nodes
    df = pd.merge(df, reduced_node_table, copy=True, left_on=['target'], right_on=['shared name'])
    df = df.rename(columns={'clusterONE': 'clusterONE_target'})

    # Create a new graph using cluster numbers as source and target node names
    source_target_df = df[['clusterONE_source', 'clusterONE_target']]
    print(source_target_df)
    source_target_df.to_csv(os.path.join(os.getcwd(), root_dir, 'results', 'test.csv'))
    # Get rid of self loops
    source_target_df.drop(source_target_df.index[source_target_df['clusterONE_source'] ==
                                                 source_target_df['clusterONE_target']], inplace=True)

    source_target_df.dropna(inplace=True)

    # Draw a graph
    g = nx.from_pandas_edgelist(source_target_df, source='clusterONE_source', target='clusterONE_target')
    # Save graph
    nx.write_graphml_xml(g, path=os.path.join(os.getcwd(), root_dir, 'results', 'pivot.xml'))


if __name__ == '__main__':
    main()