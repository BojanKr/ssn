import pandas as pd
from py2cytoscape import cyrest
import os
import math

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

def visualize_predictions(predictions_file, pivot_network):

    """ Creates and saves visualization of predictions for each cluster """

    print('views')
    print(cytoscape.view.list())

    with open(predictions_file, 'r') as f:
        for line in f.readlines():
            if line.startswith('Cluster list'):
                cluster_list = line.split(': ')[1]
                cluster_list = cluster_list.split(', ')
                cluster_list = [x.strip() for x in cluster_list]
                print(predictions_file, cluster_list)
                # Convert list values to float in order to match cell values in cytoscape table
                if len(cluster_list) > 1:
                    cluster_list = [float(x) for x in cluster_list]
                    cluster_list = [str(x) for x in cluster_list]
                elif len(cluster_list) == 1:
                    if cluster_list[0] == '':
                        continue
                    else:
                        cluster_list = [float(x) for x in cluster_list]
                        cluster_list = [str(x) for x in cluster_list]
    print('ovde')
    print(cluster_list)

    # Update visual styles
    NODE_BORDER_WIDTH = cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_BORDER_WIDTH',
                                                           mappingType='discrete',
                                                           mappingColumn='shared name',
                                                           discrete=[cluster_list, ['10'] * len(cluster_list)])
    NODE_BORDER_PAINT = cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_BORDER_PAINT',
                                                           mappingType='discrete',
                                                           mappingColumn='shared name',
                                                           discrete=[cluster_list, ['#FF1919'] * len(cluster_list)])
    NODE_LABEL = cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL", mappingType="passthrough",
                                                    mappingColumn="shared name")
    NODE_LABEL_FONT_SIZE = cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL_FONT_SIZE", mappingType="discrete",
                                                              mappingColumn="shared name",
                                                              discrete=[cluster_list, ['15'] * len(cluster_list)])
    cytoscape.vizmap.update_style(title='pivot_SSN', mappings=[NODE_BORDER_WIDTH, NODE_BORDER_PAINT, NODE_LABEL, NODE_LABEL_FONT_SIZE])
    cytoscape.vizmap.apply(styles="pivot_SSN")

    # Export view as svg
    output_file = predictions_file[:-3] + 'png'
    cytoscape.view.export(outputFile=output_file, options='PNG')


def predictions(file, e_value, directory, cluster):

    """ Writes predictions to file based on minimum E-value that indicates significance """

    df = pd.read_fwf(file)
    node_table = pd.read_csv(os.path.join(os.getcwd(), root_dir, 'data', 'node_table.csv'))
    predictions_file = directory + '/' + 'predictions.txt'
    predictions_file_summary = os.path.join(os.getcwd(), root_dir, 'results', 'function_prediction.txt')

    # Initiate the file to write predictions
    with open(predictions_file, 'a') as f:
        f.write(f'MAST search hits for {cluster} fingerprint whose E-value is better than the set threshold\n')
        f.write('\n')

    # Initiate the file to write summary of function prediction to
    with open(predictions_file_summary, 'a') as f:
        f.write(f'Function prediction based on searching the {cluster} against the dataset.\n')
        f.write('\n')
        f.write('\n')
        f.write('Clusters that might have the same or similar function compared to this cluster:\n')

    # Write all hits better tha the set E-value threshold to the file
    for index,row in df.iterrows():
        if float(row['E-VALUE']) < float(e_value):
            acc_number = row['SEQUENCE NAME'].split('|')[1] + '\n'
            with open(predictions_file, 'a') as f:
                f.write(acc_number)

    # Write cluster numbers found to match the fingerprint searched against the dataset
    with open(predictions_file, 'r') as f:
        acc_list = f.readlines()
        acc_list = [acc.strip() for acc in acc_list]

    # Keep only rows of node table that contain significant hits
    significant_hits_df = node_table.loc[node_table['shared name'].isin(acc_list)]
    # Create a list of cluster to which the hits correspond
    cluster_list = list(set(significant_hits_df['clusterONE'].tolist()))
    # Remove cluster number belonging to the cluster whose fingerprint was used in search
    if float(cluster[-1:]) in cluster_list:
        cluster_list.remove(float(cluster[-1:]))
    # Remove possible nan values from the list
    clean_cluster_list = [x for x in cluster_list if not math.isnan(x)]

    # Write cluster numbers which are found to contain similar motifs to those searched against the dataset
    with open(predictions_file, 'a') as f:
        f.write('\n')
        f.write(f'Cluster found to contain motifs similar to those of cluster {cluster[-1:]} are:\n')
        f.write('Cluster list: ')
        # Format cluster numbers before writing to file
        clust_list_to_write = list(map(int, clean_cluster_list))
        clust_list_to_write = list(map(str, clust_list_to_write))
        # Write
        f.write(', '.join(clust_list_to_write) + '\n')
        f.write('\n')

        # Write accession numbers that correspond to each hit cluster
        for clust_num in clean_cluster_list:
            hits = significant_hits_df.loc[significant_hits_df['clusterONE'] == clust_num]
            hits_list = hits['shared name'].tolist()
            with open(predictions_file, 'a') as f:
                f.write(f'Cluster {str(int(clust_num))}:\n')
                f.write(', '.join(hits_list))
                f.write('\n')

            with open(predictions_file_summary, 'a') as f:
                f.write(f'Cluster {str(clust_num)}\n')

    with open(predictions_file_summary, 'a') as f:
        f.write('\n')
        f.write('\n')


def sequences_of_interest(cluster, file):

    """ Writes hits from mast search to a file """

    lines = []
    out_file = os.path.join(os.getcwd(), root_dir, 'results', 'mast', cluster) + '/' + 'sequences.txt'
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('SECTION I'):
                break
        for line in f:
            if line.startswith('SECTION II: MOTIF DIAGRAMS'):
                break
            lines.append(line)
            with open(out_file, 'a') as f:
                if line.startswith('SEQUENCE'):
                    f.write(line)
                elif line.startswith('tr|'):
                    f.write(line)
                elif line.startswith('sp|'):
                    f.write(line)
                else:
                    continue

    return out_file


def main():

    # Remove all existing networks in the current cytoscape session
    for network in cytoscape.network.list()['networks']:
        network_id = 'SUID:' + str(network)
        cytoscape.network.destroy(network=network_id)

    mast_directory = os.path.join(os.getcwd(), root_dir, 'results', 'mast/')

    # Load pivot network to cytoscape (it will be used for visualiziing mast rezults)
    pivot = os.path.join(os.getcwd(), root_dir, 'results', 'pivot.xml')
    pivot_cy = cytoscape.network.load_file(pivot)
    # Get SUID for pivot network
    pivot_suid = 'SUID:' + str(pivot_cy['networks'][0])
    # Apply layout to pivot networks
    cytoscape.layout.apply_preferred(network=pivot_suid)

    # Create a visual style
    default_style_dict = {"NODE_SHAPE": "ellipse", "NODE_SIZE": "60", "NODE_FILL_COLOR": "#99CCFF",
                          "EDGE_TRANSPARENCY": "200", 'NODE_LABEL_FONT_SIZE':'20'}
    defaults = cytoscape.vizmap.simple_defaults(default_style_dict)
    NODE_LABEL = cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL", mappingType="passthrough",
                                                    mappingColumn="shared name")
    cytoscape.vizmap.create_style(title="pivot_SSN", defaults=defaults, mappings=[NODE_LABEL])
    cytoscape.vizmap.apply(styles="pivot_SSN")

    # Write and visualize cluster function prediction based on mast search
    for dir in os.listdir(mast_directory):
        if dir.startswith('Cluster'):
            mast_file = mast_directory  + dir + '/' + 'mast.txt'
            e_value = input(f'E-value for {dir}: ')
            output_file = sequences_of_interest(cluster = dir, file=mast_file)

            # Read the sequences of interest into a dataframe
            predictions(file=output_file, e_value=e_value, directory=mast_directory+'/'+dir, cluster=dir)

            # Visualize predictions
            predictions_file = mast_directory + '/' + dir + '/' + 'predictions.txt'
            visualize_predictions(predictions_file=predictions_file, pivot_network=pivot_suid)


if __name__ == '__main__':
    main()