import os
import pandas as pd
from py2cytoscape import cyrest
import matplotlib
from matplotlib import colors
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math

# Create cyrest client
cytoscape = cyrest.cyclient()


def table_multiindex(df, index_columns):

    """ Sets multiindex to the given table and returns the new df"""

    df.set_index(keys=index_columns, inplace=True)


def get_number_of_clusters(df):

    """ Returns total number of clusters in the SSN """

    number = len(df.index.get_level_values(0).unique().sort_values()) - 1

    return number


def swiss_prot_df(df):

    """ Retruns a dataframe containing only rows with sequences from SwissProt """

    swiss_table = df[df['UniProt Annotation Status'] == 'SwissProt']['UniProt Annotation Status']

    return swiss_table


def swiss_prot_number(df):

    """ Returns the number of clusters that contain at least one sequence from SwissProt """

    swiss_num = len(df.index.get_level_values(0).sort_values().unique())

    return swiss_num


def trembl_number(number_of_clusters, swiss_num):

    """ Returns the number of clusters that contain only seqeunces from TrEMBL """

    trembl_num = number_of_clusters - swiss_num

    return trembl_num


def swissprot_seq_num(df, root_dir):

    """ Return the number of sequences from SwissProt
        and the list of corresponding accession numbers and cluster numbers """

    # Get number of defined clusters
    swiss_seq_num = len(df.index.get_level_values(1))
    # Get cluster numbers of defined clusters
    swiss_clust_num = df.index.get_level_values(0).sort_values().unique().tolist()
    # Remove nan values from the list
    for number in swiss_clust_num:
        if math.isnan(number) is True:
            swiss_clust_num.remove(number)
    # Save the list to a file
    with open(os.path.join(os.getcwd(), root_dir, 'data', 'swiss_cluster_numbers.txt'), 'w') as f:
        f.write(', '.join(str(int(num)) for num in swiss_clust_num))
    # Get a list of accession numbers
    swiss_acc_list = df.index.get_level_values(1).sort_values().unique().to_series().to_string()

    return swiss_seq_num, swiss_clust_num, swiss_acc_list


def seqence_number(df):

    """ Returns total number of sequences from the SSN """

    seq_num = len(df.index.get_level_values(1))

    return seq_num


def ec(df):

    """ Returns all rows from the node table that contain EC numbers with all four digits
        and returns a list of clusters that contain seqeunces from SwissProt and four digit EC numbers """

    # Get a table without the None values for the EC column
    ec_df = df[df['EC'] != 'None'][['UniProt Annotation Status', 'EC']].sort_index()

    # To get rid of empty strings in the EC column
    ec_df = ec_df[ec_df['EC'] != '']

    # To keep only those EC numbers that are assigned to the manually curated UniProt entries
    ec_df = ec_df[ec_df['UniProt Annotation Status'] == 'SwissProt']

    # Get rid of anything but the EC number in each row of the EC column
    # (it is possible that some proteins have more than one EC number assigned, which has to be checked later)
    ec_df['EC'] = ec_df['EC'].str.split(' ').str.get(0)

    # Get rid of all x.x.x.- numbers
    # (they appear in multiple clusters since the first three or less numbers
    # tend to be shared between many members of the superfamily)
    for index, row in ec_df.iterrows():
        if row.EC[-1] == '-':
            ec_df.drop(index, inplace=True)

    # List of cluster numbers that contain sequences from SwissProt
    swiss_clusters = ec_df.index.get_level_values(0).unique()

    return ec_df, swiss_clusters


def count_taxonomy(df):

    """ Returns number of occurences of a given taxonomical categories in a network """

    # Count of taxonomical categories for the entire network
    kingdom_count = df['Kingdom'].value_counts()
    phylum_count = df['Phylum'].value_counts()

    return kingdom_count, phylum_count


def count_pfam(df):

    """ Returns number of sequences per Pfam family """

    # Count pfam families
    # Exclude all empty strings in PFAM column
    pfam_table = df[df['PFAM'] != "['']"]
    # Exclude all None strings from PFAM column
    pfam_table = pfam_table[pfam_table['PFAM'] != "['None']"]
    pfam_table = pfam_table.explode('PFAM')

    # Count
    pfam = pfam_table['PFAM'].value_counts()

    return pfam


def count_interpro_family(df):

    """ Returns number of sequences per InterPro family"""

    # Count InterPro families (per InterPro Family)
    # Exclude all empty strings
    interpro_family_table = df[df['InterPro (Family)'] != "['']"]
    # Exclude all None strings
    interpro_family_table = interpro_family_table[interpro_family_table['InterPro (Family)'] != "['None']"]
    interpro_family_table = interpro_family_table.explode('InterPro (Family)')
    # Count
    interpro_family = interpro_family_table['InterPro (Family)'].value_counts()

    return interpro_family


def count_interpro_domain(df):

    """ Returns number of sequences per InterPro defined domain """

    # Count InterPro domains (per InterPro Domain)
    # Exclude all empty strings
    interpro_domain_table = df[df['InterPro (Domain)'] != "['']"]
    # Exclude all None strings
    interpro_domain_table = interpro_domain_table[interpro_domain_table['InterPro (Domain)'] != "['None']"]
    interpro_domain_table = interpro_domain_table.explode('InterPro (Domain)')
    # Count
    interpro_domain = interpro_domain_table['InterPro (Domain)'].value_counts()

    return interpro_domain


def get_pdb_codes_per_cluster(df, outfile):

    """ Writes a file with all PDB codes corresponding to given accession numbers
        and returns a list of all cluster numbers that contain proteins with solved strucutres """

    # Get list of cluster numbers with solved pdb strucutres
    pdb_clust_list = []
    cluster_numbers = df.index.get_level_values(0)
    for number in cluster_numbers.dropna().unique():
        if int(number) not in pdb_clust_list:
            pdb_clust_list.append(int(number))

    # Write pdb codes to file
    with open(outfile, 'a') as f:
        for index, row in df.sort_index().iterrows():
            if not math.isnan(index[0]):
                f.write(index[1] + ' (cluster ' + str(int(index[0])) + ')' + ': ' + ', '.join(row['PDB']) + '\n')

    return pdb_clust_list


def count_pdb(cluster_list, df):

    """ Returns number of clusters with solved structures
        and number of proteins with solved structure """

    number_of_pdb_clusters = len(cluster_list)
    number_of_pdb_seq = len(df.index.get_level_values(0))

    return number_of_pdb_clusters, number_of_pdb_seq


def get_pdb_table(df):

    """ Returns a dataframe with column containing all PDB codes in the SSN """

    p_df = df['PDB'].to_frame()
    p_df.dropna(inplace=True)

    for index, row in p_df.iterrows():
        if row['PDB'][0] == 'None':
            p_df.drop(index=index, inplace=True)
        elif row['PDB'][0] == '':
            p_df.drop(index=index, inplace=True)
        elif row['PDB'][0] == 'nan':
            p_df.drop(index=index, inplace=True)
        elif math.isnan(index[0]):
            p_df.drop(index=index, inplace=True)
        else:
            continue

    return p_df


def viz_color_clusters(color_list, style, cluster_numbers, output_file):

    """ Color distinct clusters in the SSN.
        Export the SSN image. """

    NODE_FILL_COLOR = cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_FILL_COLOR', mappingType='discrete',
                                                         mappingColumn='clusterONE',
                                                         discrete=[cluster_numbers, color_list])
    cytoscape.vizmap.update_style(title=style, mappings=[NODE_FILL_COLOR])
    cytoscape.vizmap.apply(styles=style)
    cytoscape.view.export(outputFile=output_file, options='PNG')


def viz_UniProt(color_list, style, output_file):

    """ Visualize all nodes from SwissProt in the SSN.
        Export the SSN image. """

    NODE_FILL_COLOR = cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_FILL_COLOR', mappingType='discrete',
                                                         mappingColumn='UniProt Annotation Status',
                                                         discrete=[['SwissProt', 'TrEMBL'], color_list])
    NODE_SIZE = cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_SIZE', mappingType='discrete',
                                                   mappingColumn='UniProt Annotation Status',
                                                   discrete=[['SwissProt'], ['200']])
    NODE_LABEL = cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL", mappingType="passthrough",
                                                    mappingColumn='clusterONE')
    NODE_LABEL_FONT_SIZE = cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL_FONT_SIZE", mappingType="discrete",
                                                              mappingColumn='UniProt Annotation Status',
                                                              discrete=[['SwissProt'], ['100']])
    cytoscape.vizmap.update_style(title=style, mappings=[NODE_FILL_COLOR, NODE_SIZE, NODE_LABEL, NODE_LABEL_FONT_SIZE])
    cytoscape.vizmap.apply(styles=style)
    cytoscape.view.export(outputFile=output_file, options='PNG')


def viz_pdb(df, style, output_file):

    """ Visualize all nodes that represent proteins with solved structure.
        Export the SSN image. """

    # Get a list of nodes that have PDB codes associated to them
    p_df = df[['PDB', 'shared name', 'clusterONE']]
    p_df['PDB'] = p_df['PDB'].str.get(0)
    p_df['PDB'].replace('', value=np.nan, inplace=True)
    p_df['PDB'].replace('None', value=np.nan, inplace=True)
    p_df.dropna(inplace=True)
    pdb_list = p_df['shared name'].tolist()

    NODE_BORDER_WIDTH = cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_BORDER_WIDTH', mappingType='discrete',
                                                           mappingColumn='shared name',
                                                           discrete=[pdb_list, ['50'] * len(pdb_list)])
    NODE_SIZE = cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_SIZE', mappingType='discrete',
                                                   mappingColumn='shared name',
                                                   discrete=[pdb_list, ['300'] * len(pdb_list)])
    NODE_FILL_COLOR = cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_FILL_COLOR', mappingType='discrete',
                                                         mappingColumn='shared name',
                                                         discrete=[pdb_list, ['#ff4019'] * len(pdb_list)])
    NODE_LABEL = cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL", mappingType="passthrough",
                                                    mappingColumn='clusterONE')
    NODE_LABEL_FONT_SIZE = cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL_FONT_SIZE", mappingType="discrete",
                                                               mappingColumn='shared name',
                                                               discrete=[pdb_list, ['100']*len(pdb_list)])

    cytoscape.vizmap.update_style(title=style, mappings=[NODE_BORDER_WIDTH, NODE_SIZE, NODE_FILL_COLOR, NODE_LABEL, NODE_LABEL_FONT_SIZE])
    cytoscape.vizmap.apply(styles=style)
    cytoscape.view.export(outputFile=output_file, options='PNG')


def viz_taxonomy(df, color_list, style, output_file):

    """ Color nodes by Taxonomy (Kingdom).
        Export the SSN image. """

    NODE_FILL_COLOR = cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_FILL_COLOR', mappingType='discrete',
                                                         mappingColumn='Kingdom',
                                                         discrete=[list(set(df['Kingdom'].dropna().tolist())), color_list])
    NODE_SIZE = cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_SIZE', mappingType='discrete',
                                                   mappingColumn='Sequence Source',
                                                   discrete=[['USER'], ['60'] * len(list(set(df['Sequence Source'])))])
    NODE_LABEL_FONT_SIZE = cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL_FONT_SIZE",
                                                              mappingType="discrete",
                                                              mappingColumn='Sequence Source',
                                                              discrete=[['USER'], ['20']])
    cytoscape.vizmap.update_style(title=style, mappings=[NODE_SIZE, NODE_FILL_COLOR, NODE_LABEL_FONT_SIZE])
    cytoscape.vizmap.apply(styles=style)
    cytoscape.view.export(outputFile=output_file, options='PNG')


def define_colors_for_discrete_mapping(df, column):

    """ Returns a list of colors for each of unique key in the selected column """

    # Replace all empty strings in the column by nan so they can be removed together with other nan values
    df[column].replace('', value=np.nan, inplace=True)
    df[column].replace('None', value=np.nan, inplace=True)
    df[column].replace('NA', value=np.nan, inplace=True)

    # Color list in which to store color codes
    color_list = []

    # Generate color codes
    cmap = matplotlib.cm.get_cmap(lut=len(list(set(df[column].dropna().tolist()))), name='hsv')
    if len(list(set(df[column].dropna().tolist()))) > 4:
        bounds = np.array([i for i in range(1, len(list(set(df[column].dropna().tolist()))) + 1)])
        norm = matplotlib.colors.BoundaryNorm(boundaries=bounds, ncolors=len(list(set(df[column].dropna().tolist()))))
        for i in range(len(list(set(df[column].dropna().tolist())))):
            c = matplotlib.colors.rgb2hex(cmap(norm(i)))
            color_list.append(c)
    else:
        # Use these two colors if there are only two values for which to generate color codes
        color_list = ['#ff4019', '#4da6ff', '#87ff5c', '#59d3ff']

    return color_list


def merge_clusterONE(df):

    """ Cleans clusterONE results and loads them to Cytoscape
        merging them with the node table of the current session """

    dictionary = {}

    df.set_index('shared name', inplace=True)
    df['clusterONE'] = df['clusterONE'].astype(str)

    for j in df.index.tolist():
        v = df.loc[j, 'clusterONE']
        if type(v) == str:
            dictionary[j] = v
        else:
            dictionary[j] = v[1]

    clean_clusterONE_ser = pd.Series(dictionary)
    clean_clusterONE_ser = clean_clusterONE_ser.to_frame('clusterONE')
    clean_clusterONE_ser['clusterONE'] = clean_clusterONE_ser['clusterONE'].astype(int)

    # Load clusterONE results to Cytoscape node table
    cytoscape.table.loadTableData(clean_clusterONE_ser)


def get_column_names():

    """ Checks if the table is a node table ot some other table from Cytoscape and if it is node table, returns,
        the list of columns is returned. """

    # List all table ids from Cytoscape
    tables = cytoscape.table.list_tables()

    for table_id in tables['tables']:
        table_name = 'SUID:' + str(table_id)
        columns = cytoscape.table.list_columns(table_name)
        # check if column names belong to node table
        if 'PDB' in columns:
            try:
                columns.remove('name')
            except:
                continue

            return columns

        else:
            continue


def get_node_table(column_names, network):

    """ Returns the node table from Cytoscape """

    table = cytoscape.table.getTable(columns=column_names, table='node', network=network)

    return table


def check_isofunctionality(df, cluster_number):

    """ Checks whether a cluster is isofunctional

        Inputs are a dataframe with EC numbers of SwissProt members
        in the network and the number of a cluster whose isofunctionality
        is being checked. """

    # Seqeunces that were not found to be in any of the clusters by clusterONE
    # are all grouped under nan and they shouldn't be considered
    if str(cluster_number) != 'nan':
        if df['EC'].loc[cluster_number].nunique() == 1:
            isofunctionality = 'ISOFUNCTIONAL'
            ec_num = df['EC'].loc[cluster_number].unique()[0]

            return isofunctionality, ec_num

        else:
            isofunctionality = 'NOT isofunctional'
            ec_num = ', '.join(df['EC'].loc[cluster_number].unique())

            return isofunctionality, ec_num


def get_visual_properties(df, column, property_name):

    """ Returns a dictionary (value:color) for a given set of nodes """

    # Get a list of all values from the given column
    df[column].replace('', value=np.nan, inplace=True)
    df[column].replace('None', value=np.nan, inplace=True)
    df[column].replace('NA', value=np.nan, inplace=True)
    list_of_values = list(set(df[column].dropna().tolist()))

    dictionary = {}

    for value in list_of_values:
        node_list = column + ':' + value
        properties = cytoscape.node.get_properties(nodeList=node_list, propertyList=property_name)
        color = properties[0]['visualProperties'][0]['value']
        dictionary[value] = color

    return dictionary


def legend(visual_properties, file):

    """ Creates a legend based on value:color dictionary """

    f = lambda m, c: plt.plot([], [], marker=m, color=c, ls="none")[0]
    handles = [f("s", visual_properties[key]) for key in visual_properties.keys()]
    labels = visual_properties.keys()
    fig_legend = plt.legend(handles, labels, loc=10, framealpha=5, frameon=False)
    fig = fig_legend.figure
    plt.gca().set_axis_off()
    fig.canvas.draw()

    fig.savefig(file, dpi=300, bbox_inches='tight')


def write_general_info(file, clust_num, swiss_num, swiss_cluster_numbers, swiss_seq_num, seq_num, trembl_num,
                       swiss_clusters, ec_number, kingdom_count, phylum_count, pfam, interpro_domain, interpro_family,
                       prot_struct_num, clust_struct_num):

    """ Writes results to a specified file. """

    with open(file + 'general_info.txt', 'a') as f:
        f.write('GENERAL INFORMATION ABOUT IDENTIFIED CLUSTERS\n')
        f.write('\n')
        f.write('Number of clusters: {}\n'.format(clust_num))
        f.write('Number of defined clusters: {}\n'.format(swiss_num))
        f.write('\n')
        f.write('Defined clusters are:\n')
        f.write(', '.join(str(int(num)) for num in swiss_cluster_numbers))
        f.write('\n')
        f.write('\nNumber of sequences from SwissProt: {}\n'.format(swiss_seq_num))
        f.write('Percentage of SwissProt seqeunces in the network: {}%\n'.format(round(swiss_seq_num / seq_num * 100, 3)))
        f.write('\n')
        f.write('Number of nondefined clusters: {}\n'.format(trembl_num))
        f.write('\n')
        f.write(200 * '#')
        f.write('\n')
        f.write(200 * '#')
        f.write('\n')
        f.write('\n')
        f.write('CHECKING CLUSTER ISOFUNCTIONALITY\n')
        f.write('\n')
        for cluster_number in swiss_clusters:
            if str(cluster_number).lower() != 'nan':
                isofunctionality, ec_num = check_isofunctionality(ec_number, cluster_number)
                if isofunctionality == 'ISOFUNCTIONAL':
                    f.write('Cluster {} is {}\n'.format(str(cluster_number), isofunctionality))
                    f.write('EC number in the cluster: {}\n'.format(ec_num))
                    f.write('\n')
                elif isofunctionality == 'NOT isofunctional':
                    f.write('Cluster {} is {}\n'.format(str(cluster_number), isofunctionality))
                    f.write('EC number in the cluster: {}\n'.format(ec_num))
                    f.write('\n')
        f.write('\n')
        f.write(200 * '#')
        f.write('\n')
        f.write(200 * '#')
        f.write('\n')
        f.write('\n')
        f.write('TAXONOMY INFORMATION\n')
        f.write('\n')
        f.write('Number of representatives per kingdom is:\n')
        if kingdom_count.size == 0:
            f.write('There is no information about Kingdom taxonomical category')
        else:
            f.write(kingdom_count.to_string())
        f.write('\n')
        f.write('\n')
        f.write('Number of representatives per phylum is:\n')
        if phylum_count.size == 0:
            f.write('There is no information about Phylum taxonomical category')
        else:
            f.write(phylum_count.to_string())
        f.write('\n')
        f.write('\n')
        f.write(200 * '#')
        f.write('\n')
        f.write(200 * '#')
        f.write('\n')
        f.write('\n')
        f.write('PFAM AND INTERPRO FAMILIES\n')
        f.write('\n')
        f.write('Number of seqeunces per Pfam family:\n')
        f.write(pfam.to_string())
        f.write('\n')
        f.write('\n')
        f.write('Number of seqeunces per InterPro defined domain:\n')
        f.write(interpro_domain.to_string())
        f.write('\n')
        f.write('\n')
        f.write('Number of seqeunces per InterPro defined family:\n')
        f.write(interpro_family.to_string())
        f.write('\n')
        f.write('\n')
        f.write(200 * '#')
        f.write('\n')
        f.write(200 * '#')
        f.write('\n')
        f.write('\n')
        f.write('PDB IN CLUSTERS\n')
        f.write('\n')
        f.write('Number of proteins with solved structure is: {}\n'.format(prot_struct_num))
        f.write('Number of clusters that contain solved structures is {}\n'.format(clust_struct_num))
        f.write('List of PDBs written to {}\n'.format(file + 'pdb.txt'))
        f.write('\n')
        f.write(200 * '#')
        f.write('\n')
        f.write(200 * '#')
        f.write('\n')
        f.write('\n')
        f.write('GENERATED IMAGES\n')
        f.write('\n')
        f.write('Image 1 is an image of the full SSN topology.\n')
        f.write('\n')
        f.write('Image 2 is an image of the full SSN with each separate cluster, identified by ClusterONE algorithm,'
                'is colored differently.\n')
        f.write('\n')
        f.write('Image 3 shows the position of sequences from SwissProt and TrEMBL. Each SwissProt node is labeled'
                ' with a corresponding cluster number. '
                'Color coding for this image is explained in Image_3_legend.png.\n')
        f.write('\n')
        f.write('Image 4 shows the distribution of the Kingdom taxonomical category in the SSN.'
                'Color coding is explained in Image_4_legend.png. \n')
        f.write('\n')
        f.write('Image 5 shows proteins with solved structure in the SSN. Nodes that represent proteins proteins '
                'with solved structure are colored red and labeled with their corresponding cluster number.\n')
        f.write('\n')
        f.write(200 * '#')
        f.write('\n')
        f.write(200 * '#')
        f.write('\n')


def main():

    print('Extracting general information about the full SSN.\n')

    root_dir = os.listdir(os.path.join(os.getcwd(), 'network/'))[0]

    # clusterONE results dataframe
    clusterONE_results = pd.read_csv(os.path.join(os.getcwd(), root_dir, 'results', 'sorted_clusterONE_results.csv'))

    # Merge clusterONE results with Cytoscape node table
    merge_clusterONE(df=clusterONE_results)

    # Get the list ids for networks in the current cytoscape session
    networkSUID_list = cytoscape.network.list()
    # Get the network name for the given netwok (SUID:id)
    network_name = 'SUID:' + str(networkSUID_list['networks'])
    column_names = get_column_names()
    # Retrieve the table from Cytoscape
    node_table = get_node_table(column_names=column_names, network=network_name)
    node_table.to_csv(os.path.join(os.getcwd(), root_dir, 'data', 'node_table.csv'))

    # Create visual style
    cytoscape.vizmap.apply(styles='default')
    NODE_LABEL = cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL", mappingType="passthrough",
                                                    mappingColumn="clusterONE")
    default_style_dict = {"NODE_SHAPE": "ellipse", "NODE_SIZE": "60", "NODE_FILL_COLOR": "#00F00",
                          "EDGE_TRANSPARENCY": "200"}
    defaults = cytoscape.vizmap.simple_defaults(default_style_dict)
    cytoscape.vizmap.create_style(title="SSN", defaults=defaults, mappings=[NODE_LABEL])
    # Apply style
    cytoscape.vizmap.apply(styles="SSN")
    # First image of the full SSN
    cytoscape.view.export(outputFile=os.path.join(os.getcwd(), root_dir, 'results', 'images', 'Image_1.png'), options='PNG')

    # Color clusters
    clust_numbers = list(set(clusterONE_results['clusterONE'].tolist()))
    clust_colors = define_colors_for_discrete_mapping(node_table, column='clusterONE')
    out_file = os.path.join(os.getcwd(), root_dir, 'results', 'images', 'Image_2.png')
    viz_color_clusters(color_list=clust_colors, style='SSN', cluster_numbers=clust_numbers, output_file=out_file)

    # Visualize UniProt
    uniprot_colors = define_colors_for_discrete_mapping(node_table, 'UniProt Annotation Status')
    uniprot_out_file = os.path.join(os.getcwd(), root_dir, 'results', 'images', 'Image_3.png')
    viz_UniProt(color_list=uniprot_colors, style='SSN', output_file=uniprot_out_file)
    vp_u = get_visual_properties(df=node_table, column='UniProt Annotation Status', property_name='Fill Color')
    legend_out_file = uniprot_out_file[:-4] + '_legend.png'
    legend(visual_properties=vp_u, file=legend_out_file)

    # Visualize taxonomy
    taxonomy_colors = define_colors_for_discrete_mapping(node_table, 'Kingdom')
    taxonomy_out_file = os.path.join(os.getcwd(), root_dir, 'results', 'images', 'Image_4.png')
    viz_taxonomy(df=node_table, color_list=taxonomy_colors, style='SSN', output_file=taxonomy_out_file)
    vp_t = get_visual_properties(df=node_table, column='Kingdom', property_name='Fill Color')
    legend_out_file = taxonomy_out_file[:-4] + '_legend.png'
    legend(visual_properties=vp_t, file=legend_out_file)

    # Visualize clusters with solved structure
    pdb_out_file = os.path.join(os.getcwd(), root_dir, 'results', 'images', 'Image_5.png')
    viz_pdb(df=node_table, style='SSN', output_file=pdb_out_file)

    # Set multiindex for the new table
    table_multiindex(df=node_table, index_columns=['clusterONE', 'shared name'])

    clusters = get_number_of_clusters(node_table)

    swiss_table = swiss_prot_df(node_table)

    swiss_num = swiss_prot_number(df=swiss_table)

    trembl_num = trembl_number(number_of_clusters=clusters, swiss_num=swiss_num)

    swiss_seq_num, swiss_clust_num, swiss_acc_list = swissprot_seq_num(df=swiss_table, root_dir=root_dir)

    seq_num = seqence_number(node_table)

    ec_numbers, swiss_clusters = ec(node_table)

    kingdom_count, phylum_count = count_taxonomy(node_table)

    pfam = count_pfam(node_table)

    interpro_family = count_interpro_family(node_table)
    interpro_domain = count_interpro_domain(node_table)

    pdb_df = get_pdb_table(node_table)
    pdb_out_file = os.path.join(os.getcwd(), root_dir, 'results', 'pdb.txt')
    pdb_clust_list = get_pdb_codes_per_cluster(df=pdb_df, outfile=pdb_out_file)

    number_of_pdb_clusters, number_of_pdb_seq = count_pdb(cluster_list=pdb_clust_list, df=pdb_df)

    # Write general SSN information
    general_info_output_dir = os.path.join(os.getcwd(), root_dir, 'results/')

    write_general_info(file=general_info_output_dir, clust_num=clusters, swiss_num=swiss_num,
                       swiss_cluster_numbers=swiss_clust_num, swiss_seq_num=swiss_seq_num, seq_num=seq_num,
                       trembl_num=trembl_num, swiss_clusters=swiss_clusters, ec_number=ec_numbers,
                       kingdom_count=kingdom_count, phylum_count=phylum_count, pfam=pfam,
                       interpro_domain=interpro_domain, interpro_family=interpro_family,
                       prot_struct_num=number_of_pdb_seq, clust_struct_num=number_of_pdb_clusters)

    print(f'General information about the full SSN written to {general_info_output_dir}general_info.txt\n')

if __name__ == '__main__':
    main()