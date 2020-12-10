import pandas as pd
import os
from py2cytoscape import cyrest
import warnings

# Create cyrest client
cytoscape = cyrest.cyclient()


def get_root_dir():
    for file in os.listdir(os.path.join(os.getcwd(), 'network/')):
        if not file.endswith('xgmml'):
            continue
        else:
            dir = os.path.join(os.getcwd(), file + '/')

            return dir


root_dir = get_root_dir()


def count_pfam_interpro(table):

    """ Returns number of sequences per Pfam family and per InterPro defined domain """

    # Count pfam and interpro families
    pfam = table['PFAM'].str[1:-1].value_counts()
    interpro_domain = table['InterPro (Domain)'].str.get(0).value_counts()

    return pfam, interpro_domain


def count_taxonomy(table):

    """ Returns number of occurences of a given taxonomical categories in a network.
        Input is a node table """

    # Count of taxonomical categories for the entire network
    kingdom_count = table['Kingdom'].value_counts()
    phylum_count = table['Phylum'].value_counts()

    return kingdom_count, phylum_count


def get_general_info(table):

    """ Get general information about the network based on the node table.
        Returns number of sequences (number of accession numbers in the table),
        number of defined sequences (from SwissProt), number of clusters with
        defined sequences and the number of unknown clusters (clusters that
        contain only sequences from TrEMBL).

        Input is a node table extracted from the root network or a subnetwork. """

    # Count number of clusters
    clust_num = len(table.index.get_level_values(0).unique().sort_values()) - 1

    # Get number of defined (swissprot) and unknown (trembl) clusters
    swiss_table = table[table['UniProt Annotation Status'] == 'SwissProt']['UniProt Annotation Status']
    # number of clusters with swissprot entries
    swiss_num = len(swiss_table.index.get_level_values(0).sort_values().unique()) - 1
    trembl_num = clust_num - swiss_num  # number of unknown clusters

    # Get number of swissprot and trembl sequences in the network
    swiss_acc = swiss_table.index.get_level_values(1)
    swiss_seq_num = len(swiss_acc)  # number of sequences from swissprot
    seq_num = len(table.index.get_level_values(1))  # total number of sequences in the network

    # Turn swiss prot cluster values into format that can be writen into a text file
    swiss_cluster_numbers = swiss_table.index.get_level_values(0).sort_values().unique().to_frame(index=False).to_string()

    return clust_num, swiss_num, trembl_num, seq_num, swiss_cluster_numbers, swiss_acc, swiss_seq_num


def sequence_info(table):

    """ Returns inforamtion about sequence length for each of the proteins in the network """

    avg_seq_length = table['Sequence Length'].mean()
    min_seq_length = table['Sequence Length'].min()
    max_seq_length = table['Sequence Length'].max()
    mode_seq_length = table['Sequence Length'].mode()
    median_seq_length = table['Sequence Length'].median()

    return avg_seq_length, min_seq_length, max_seq_length, mode_seq_length, median_seq_length


def number_of_edges(table):

    """ Returns number of edges in cluster """

    # number of edges in the network
    count_number_of_edges = len(table.index.get_level_values(0))

    return count_number_of_edges


def percentage_info(table):

    """ Returns information about percentage identity based on edge etwork """

    # percentage identity
    avg_percentage_id = table['%id'].mean()
    min_percentage_id = table['%id'].min()
    mode_percentage_id = table['%id'].mode()
    median_percentage_id = table['%id'].median()

    return avg_percentage_id, min_percentage_id, mode_percentage_id, median_percentage_id


def alignment_score_info(table):

    """ Returns information about alignment scores between connected nodes of a cluster """

    avg_alignment_score = table['alignment_score'].mean()
    min_alignment_score = table['alignment_score'].min()
    max_alignment_score = table['alignment_score'].max()
    mode_alignment_score = table['alignment_score'].mode()
    median_alignment_score = table['alignment_score'].median()

    return avg_alignment_score, min_alignment_score, max_alignment_score, mode_alignment_score, median_alignment_score


def alignment_length_info(table):

    """ Returns information about alignment lengths between connected nodes of a cluster """

    avg_alignment_length = table['alignment_len'].mean()
    min_alignment_length = table['alignment_len'].min()
    max_alignment_length = table['alignment_len'].max()
    mode_alignment_length = table['alignment_len'].mode()
    median_alignment_length = table['alignment_len'].median()

    return avg_alignment_length, min_alignment_length, max_alignment_length, mode_alignment_length, median_alignment_length


def cluster_conectivity_node(table, network):

    """ Returns information about how well is the rest of the cluster connected
        to its defined sequences (if there are any defiende sequences present).

        Inputs are a node table for the given network and the network id """

    # Number of nodes
    table = table.reset_index()
    total_number_of_nodes = len(table.index.get_level_values(0))

    # Select all nodes from SwissProt in the given cluster
    swiss_nodes_select = cytoscape.network.select(nodeList='UniProt Annotation Status:SwissProt',
                                                  firstNeighbors='any', network=network)
    if swiss_nodes_select is not None:
        # count the total number of selected nodes
        total_number_of_selected_nodes = len(swiss_nodes_select['nodes'])
        # count the number of SwissProt nodes
        swiss_prot = table[table['UniProt Annotation Status'] == 'SwissProt']['UniProt Annotation Status']
        number_of_swiss_nodes = len(swiss_prot.index.get_level_values(0))
        number_of_connected_trembl_nodes = total_number_of_selected_nodes - number_of_swiss_nodes
        deselect = cytoscape.network.deselect(edgeList='all', nodeList='all')
        return total_number_of_nodes, total_number_of_selected_nodes, number_of_swiss_nodes, number_of_connected_trembl_nodes
    else:
        total_number_of_selected_nodes = 0
        number_of_swiss_nodes = 'There are no nodes from SwissProt in this cluster'
        number_of_connected_trembl_nodes = 0
        deselect = cytoscape.network.deselect(edgeList='all', nodeList='all')

        return total_number_of_nodes, total_number_of_selected_nodes, number_of_swiss_nodes, number_of_connected_trembl_nodes


def dataframe_cluster_conectivity_edge(edge_table_column_names, network, selected_edges):

    """ Returns a dataframe containing values from a particular column of the edge table (alignment length,
        alignment score, percentage id ...)

        Inputs are a network suid and a  column name for which valuse should be extracted. """

    # Get information about selected edges (only if the cluster has defined nodes)
    if selected_edges is None:
        return 'There are no defined proteins in this cluster.'
    else:
        edge_information = cytoscape.edge.get_attribute(edgeList='selected', columnList=edge_table_column_names,
                                                        network=network)
        df = pd.DataFrame(edge_information)

        return df


def cluster_connectivity_edge(df, column_name):

    """ Returns information about edges connected to defined nodes in a network

        Inputs are a dataframe containg information about all of the edges connected
        to defined nodes in the network and a name of the column to be analyzed. """

    # Information about edges
    average = df[column_name].mean()
    minimum = df[column_name].min()
    maximum = df[column_name].max()
    mode = df[column_name].mode()
    median = df[column_name].median()

    return average, minimum, maximum, mode, median


def num_selected_edges(df):

    """ Retruns a number of selected edges.

        Input is a dataframe containing all selected edges """

    # Count number of edges going from and to defined proteins
    number_of_selected_edges = len(df.index.get_level_values(0))

    return number_of_selected_edges


def general_info_per_cluster(out_dir, node_table, cluster_name):

    """ Writes general information about a subnetwork to a file (based on a node table)"""

    pfam, interpro_domain = count_pfam_interpro(table=node_table)
    kingdom_count, phylum_count = count_taxonomy(table=node_table)
    clust_num, swiss_num, trembl_num, seq_num, swiss_cluster_numbers, swiss_acc, swiss_seq_num = get_general_info(table=node_table)
    avg_seq_length, min_seq_length, max_seq_length, mode_seq_length, median_seq_length = sequence_info(table=node_table)

    with open(out_dir + cluster_name + '.txt', 'a') as f:
        f.write('\n')
        f.write('{}\n'.format(cluster_name))
        f.write('\n')
        f.write('NODE INFORMATION\n')
        f.write('\n')
        f.write('Number of sequences: {}\n'.format(str(seq_num)))
        f.write('\n')
        if swiss_seq_num > 0:
            f.write('Number of defined sequences: {}\n'.format(str(swiss_seq_num)))
            f.write('Defined sequences are {}\n'.format(' '.join(swiss_acc)))
        else:
            f.write('There are no defined sequences in this cluster, only sequences from TrEMBL.\n')
        f.write('Number of sequences from TrEMBL: {}\n'.format(str(seq_num - swiss_seq_num)))
        f.write('\n')
        f.write('Taxonomy composition of the cluster on Kingdom level:\n')
        f.write(kingdom_count.to_string())
        f.write('\n')
        f.write('\n')
        f.write('Taxonomy composition of the cluster on Phylum level:\n')
        f.write(phylum_count.to_string())
        f.write('\n')
        f.write('\n')
        f.write('Pfam families in this cluster:\n')
        f.write(pfam.to_string())
        f.write('\n')
        f.write('\n')
        f.write(200*'#')
        f.write('\n')
        f.write(200*'#')
        f.write('\n')
        f.write('\n')
        f.write('SEQUENCE LENGTH\n')
        f.write('\n')
        f.write('Average sequence length: {}\n'.format(str(round(avg_seq_length, 2))))
        f.write('Minimum sequence length in the cluster: {}\n'.format(str(min_seq_length)))
        f.write('Maximum sequence length in the cluster: {}\n'.format(str(max_seq_length)))
        f.write('The most common sequence length: {}\n'.format(str(pd.Series(mode_seq_length).array[0])))
        f.write('Median sequence length: {}\n'.format(str(median_seq_length)))
        f.write('\n')
        f.write(200*'#')
        f.write('\n')
        f.write(200*'#')
        f.write('\n')
        f.write('\n')


def write_edge_info(out_dir, node_table, edge_table, network_id, cluster_name):

    """ Writes inforamtion about edges of a cluster to a file

          Inputs are directory where the files are to br written, a node and a corresponding edge table """

    num_edges = number_of_edges(table=edge_table)
    avg_percentage_id, min_percentage_id, mode_percentage_id, median_percentage_id = percentage_info(table=edge_table)
    avg_alignment_score, min_alignment_score, max_alignment_score, mode_alignment_score, median_alignment_score = alignment_score_info(table=edge_table)
    avg_alignment_length, min_alignment_length, max_alignment_length, mode_alignment_length, median_alignment_length = alignment_length_info(table=edge_table)
    total_number_of_nodes, total_number_of_selected_nodes, number_of_swiss_nodes, number_of_connected_trembl_nodes = cluster_conectivity_node(table=node_table, network=network_id)
    edge_table_columns = cytoscape.edge.list_attributes()
    print('EDGE TABLE COLUMNS')
    print(edge_table_columns)

    with open(out_dir + cluster_name + '.txt', 'a') as f:
        f.write('EDGE INFORMATION\n')
        f.write('\n')
        f.write('Number of edges: {}\n'.format(str(num_edges)))
        f.write('\n')
        f.write('Average percentage identity: {}\n'.format(str(round(avg_percentage_id, 2))))
        f.write('Minimum percentage identity in the cluster: {}\n'.format(str(min_percentage_id)))
        f.write('Most common percentage identity:\n'.format(str(pd.Series(mode_percentage_id).array[0])))
        f.write('Median percentage identity: {}\n'.format(str(median_percentage_id)))
        f.write('\n')
        f.write('Average alignment score: {}\n'.format(str(round(avg_alignment_score, 2))))
        f.write('Minimum alignment score in the cluster: {}\n'.format(str(min_alignment_score)))
        f.write('Maximum alignment score in the cluster: {}\n'.format(str(max_alignment_score)))
        f.write('Most common alignment score: {}\n'.format(str(pd.Series(mode_alignment_length).array[0])))
        f.write('Median alignment score: {}\n'.format(str(median_alignment_score)))
        f.write('\n')
        f.write('Average alignment length: {}\n'.format(str(round(avg_alignment_length, 2))))
        f.write('Minimum alignment length in the cluster: {}\n'.format(str(min_alignment_length)))
        f.write('Maximum alignment length in the cluster: {}\n'.format(str(max_alignment_length)))
        f.write('Most common alignment length: {}\n'.format(str(pd.Series(mode_alignment_length).array[0])))
        f.write('Median alignment length: {}\n'.format(str(median_alignment_length)))
        f.write('\n')
        f.write('NODE CONNECTIVITY INFORMATION\n')
        f.write('\n')
        f.write('We transfer annotation from manually curated sequences (SwissProt) inside a cluster to all TrEMBL sequences.\n')
        f.write('It is therefore important to see the number of TrEMBL sequences directly connected to SwissProt sequences\n')
        f.write('and to see what is the percentage identity in these established connections.\n')
        f.write('\n')
        if total_number_of_selected_nodes is 0:
            f.write(number_of_swiss_nodes)
            f.write('\n')
        else:
            f.write('Number of sequences in the cluster: {}\n'.format(str(total_number_of_nodes)))
            f.write('Number of TrEMBL sequence in the first neighborhood to SwissProt sequences: {}\n'.format(
                str(number_of_connected_trembl_nodes)))
            f.write('Number of SwissProt sequenecs: {}\n'.format(str(number_of_swiss_nodes)))
            f.write('\n')
            f.write('EDGE CONNECTIVITY INFORMATION\n')
            f.write('\n')
            # Select defined nodes and edges connected to these nodes
            selected_edges = cytoscape.network.select(nodeList='UniProt Annotation Status:SwissProt',
                                                      adjacentEdges=True, network=network_id)
            for edge_table_column_name in edge_table_columns:
                if edge_table_column_name == '%id':
                    dataframe = dataframe_cluster_conectivity_edge(edge_table_column_name, network_id, selected_edges)
                    average, minimum, maximum, mode, median = cluster_connectivity_edge(dataframe, edge_table_column_name)
                    f.write('Average percentage identity for edges connected to defined nodes: {}\n'.format(str(round(average, 2))))
                    f.write('Minimum percentage identity for edges connected to defined nodes: {}\n'.format(str(minimum)))
                    f.write('Mode percentage identity in the first neighborhood to defined proteins: {}\n'.format(str(pd.Series(mode).array[0])))
                    f.write('Median percentage identity in the first neighborhood to defined proteins: {}\n'.format(str(median)))
                    f.write('\n')
                elif edge_table_column_name == 'alignment_score':
                    dataframe = dataframe_cluster_conectivity_edge(edge_table_column_name, network_id, selected_edges)
                    average, minimum, maximum, mode, median = cluster_connectivity_edge(dataframe, edge_table_column_name)
                    f.write('Average alignment score for edges connected to defined nodes: {}\n'.format(str(round(average, 2))))
                    f.write('Minimum alignment score for edges connected to defined nodes: {}\n'.format(str(minimum)))
                    f.write('Maximum alignment score for edges connected to defined nodes: {}\n'.format(str(maximum)))
                    f.write('Mode alignment score in the first neighborhood to defined proteins: {}\n'.format(str(pd.Series(mode).array[0])))
                    f.write('Median alignment score in the first neighborhood to defined proteins: {}\n'.format(str(median)))
                    f.write('\n')
                elif edge_table_column_name == 'alignment_len':
                    dataframe = dataframe_cluster_conectivity_edge(edge_table_column_name, network_id, selected_edges)
                    average, minimum, maximum, mode, median = cluster_connectivity_edge(dataframe, edge_table_column_name)
                    f.write('Average alignment length for edges connected to defined nodes: {}\n'.format(str(round(average, 2))))
                    f.write('Minimum alignment length for edges connected to defined nodes: {}\n'.format(str(minimum)))
                    f.write('Maximum alignment length for edges connected to defined nodes: {}\n'.format(str(maximum)))
                    f.write('Mode alignment length in the first neighborhood to defined proteins: {}\n'.format(str(pd.Series(mode).array[0])))
                    f.write('Median alignment length in the first neighborhood to defined proteins: {}\n'.format(str(median)))
                    f.write('\n')
                else:
                    continue

            number_of_selected_edges = num_selected_edges(df=edge_table)
            f.write('Number of edges directly connected to defined nodes in the cluster: {}\n'.format(str(number_of_selected_edges)))
            deselect_edges = cytoscape.network.deselect(nodeList='all', edgeList='all')


def main():

    # Ignore future warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)

    data_directory = os.path.join(os.getcwd(), root_dir, 'data', 'Clusters', 'tables/')
    output_directory = os.path.join(os.getcwd(), root_dir, 'results', 'Clusters/')

    print(f'Writinng information for defined clusters to {output_directory}\n')

    for file in os.listdir(data_directory):
        if file.endswith("edge_table.csv"):

                # Extract network id corresponding to the given cluster
                start = file.find('[') + 1
                end = file.find(']')
                subnetwork_id = file[start:end]

                # Get node and edge table dataframes
                edge_table_df = pd.read_csv(data_directory + file)
                node_table_df = pd.read_csv(data_directory + file[:start-1] + 'node_table.csv')
                # Set table multiindex
                node_table_df = node_table_df.set_index(keys=['clusterONE', 'shared name'])

                # Isolate cluster number for current cluster
                cluster_number_start = file.find('r_') + 2
                cluster_number_end = file.find('_[')
                cluster_number = file[cluster_number_start:cluster_number_end]
                cluster_name = "Cluster_" + cluster_number + '_info'

                # Write general information about clustrs to files
                general_info_per_cluster(out_dir=output_directory, node_table=node_table_df, cluster_name=cluster_name)
                # Write edge information
                write_edge_info(out_dir=output_directory, node_table=node_table_df, edge_table=edge_table_df,
                                network_id=subnetwork_id, cluster_name=cluster_name)

                print('Cluster {} done'.format(cluster_name))


if __name__ == '__main__':
    main()