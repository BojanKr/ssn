from py2cytoscape import cyrest
from urllib import parse
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


root_dir = get_root_dir()


def get_node_table_columns():

    """ Checks if the table is a node table ot some other table from Cytoscape and if it is node table, returns,
        the list of columns is returned. """

    # List all table ids from Cytoscape
    tables = cytoscape.table.list_tables()

    for table_id in tables['tables']:
        table_name = 'SUID:' + str(table_id)
        columns = cytoscape.table.list_columns(table_name)
        print('KOLONE')
        print(columns)
        # check if column names belong to node table
        if 'PDB' in columns:
            try:
                columns.remove('name')
            except:
                continue
            return columns
        else:
            continue


def get_edge_table_columns():

    """ Checks if the table is a node table ot some other table from Cytoscape and if it is node table, returns,
        the list of columns is returned. """

    # List all table ids from Cytoscape
    tables = cytoscape.table.list_tables()

    for table_id in tables['tables']:
        table_name = 'SUID:' + str(table_id)
        columns = cytoscape.table.list_columns(table_name)
        percent_id = parse.quote("%id")
        # check if column names belong to node table
        if 'alignment_score' in columns:
            try:
                columns.remove('name')
            except:
                continue
            columns.remove('%id')
            columns.append(percent_id)
            return columns
        else:
            continue


def write_subnetworks_edge(network, column_names, table, name):

    """ Writes edge tables of all defined clusters to csv """

    edge_table = cytoscape.table.getTable(columns=column_names, network=network, table=table)
    file_name = name + '_edge_table.csv'
    file_path = os.path.join(os.getcwd(), root_dir, 'data', 'Clusters', 'tables/') + file_name
    edge_table.to_csv(file_path)


def write_subnetworks_node(network, column_names, table, name):

    """ Writes node tables of all defined clusters to csv """

    node_table = cytoscape.table.getTable(columns=column_names, network=network, table=table)
    file_name = name + '_node_table.csv'
    file_path = os.path.join(os.getcwd(), root_dir, 'data', 'Clusters', 'tables/') + file_name
    node_table.to_csv(file_path)


def main():
    # Get network id for root network
    networkSUID_list = cytoscape.network.list()
    root_suid = networkSUID_list['networks'][0]
    network_id = 'SUID:' + str(root_suid)

    # To track which subnetworks have already been processed
    cluster_list = []
    subnetworkSUID_list = [root_suid]
    # To extract node and edge tables from each subnetwork
    subnetworkSUID_dict = {}

    for file in os.listdir(os.path.join(os.getcwd(), root_dir, 'data', 'Clusters/')):
        if file.endswith('txt'):
            # Adding cluster numbers that were processed to the list
            cluster_list.append(file[:-4])
            cluster = os.path.join(os.getcwd(), root_dir, 'data', 'Clusters/') + file

            # Create a list of COLUMN:VALUE pairs for creating subnetworks
            node_list = []
            with open(cluster, 'r') as f:
                for line in f.readlines():
                    pair = 'shared name:' + line.strip()
                    node_list.append(pair)

            # Argument for subnetwork (nodeList)
            nodes = ','.join(node_list)

            # Create subnetworks corresponding to defined clusters
            subnetwork = cytoscape.network.create(networkName=file, nodeList=nodes, source=network_id)
            subnetworkSUID = cytoscape.network.list()
            for suid in subnetworkSUID['networks']:
                if suid not in subnetworkSUID_list:
                    cluster = file[:-4]
                    subnetworkSUID_dict[cluster] = str(suid)
                    subnetworkSUID_list.append(suid)

    with open(os.path.join(os.getcwd(), root_dir, 'data', 'test.txt'), 'w') as f:
        f.write(str(subnetworkSUID_dict))

    # Get list of columns for edge and node tables (used to retrieve corresponding tables)
    edge_table_columns = get_edge_table_columns()
    node_table_columns = get_node_table_columns()

    for key, value in subnetworkSUID_dict.items():
        # Write node and edge tables for each of the defined clusters to csv
        subnetwork_id = 'SUID:' + str(value)
        name_edge = key + '_' + '[' + 'SUID:' + str(value) + ']'
        write_subnetworks_edge(network=subnetwork_id, column_names=edge_table_columns, table='edge', name=name_edge)
        print('Edge for {} done'.format(key))
        write_subnetworks_node(network=subnetwork_id, column_names=node_table_columns, table='node', name=key)
        print('Node for {} done'.format(key))
        print('Done {}'.format(key))

    # Set the root network as the current for future use
    cytoscape.network.set_current(network=network_id)


if __name__ == '__main__':
    main()
