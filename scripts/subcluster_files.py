import pandas as pd
import os


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

    root_dir = get_root_dir()

    print('Saving accession numbers of proteins that belong to defined clusters.\n')

    # Read a list of defined clusters (clusters with sequences from SwissProt)
    with open(os.path.join(os.getcwd(), root_dir, 'data', 'swiss_cluster_numbers.txt'), 'r') as f:
        cluster_list = f.readline().split(', ')

    # Read node table
    node_table = pd.read_csv(os.path.join(os.getcwd(), root_dir, 'data', 'node_table.csv'))

    list_of_clusters = []
    # Save accession numbers corresponding to defined clusters
    for i in cluster_list:
        defined_cluster = node_table[node_table['clusterONE'] == float(i)]
        defined_cluster = defined_cluster['shared name']
        out_file = os.path.join(os.getcwd(), root_dir, 'data', 'Clusters', 'Cluster_' + str(int(float(i))) + '.txt')
        with open(out_file, 'a') as f:
            for acc_num in defined_cluster:
                f.write(acc_num + '\n')
        list_of_clusters.append(i)

    cluster_list_df = pd.Series(list_of_clusters)
    cluster_list_df.to_csv(os.path.join(os.getcwd(), root_dir, 'data', 'cluster_list.csv'))

    print('Files with accession numbers written to {}'.format(os.path.join(root_dir, 'data', 'Clusters/')))

if __name__ == '__main__':
    main()
