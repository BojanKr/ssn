from py2cytoscape import cyrest
import os

# Create folders
results_folder = os.path.join(os.getcwd(), 'results/')
data_folder = os.path.join(os.getcwd(), 'data/')
file = os.listdir(os.path.join(os.getcwd(), 'network/'))[0]
print(file)
for file in os.listdir(os.path.join(os.getcwd(), 'network/')):
    project_folder = os.path.join(os.getcwd(), file+'/')
    if not os.path.exists(project_folder):
        os.makedirs(project_folder)
    else:
        new_name = project_folder[:-1]+'_copy_'+str(1)+'/'
        if not os.path.exists(new_name):
            os.rename(project_folder, new_name)
        else:
            directory_list = []
            for directory in os.listdir(os.getcwd()):
                if directory.startswith(file):
                    directory_list.append(directory)
            os.rename(project_folder, project_folder[:-1]+'_copy_'+str(len(directory_list)))
            print(f'Folder {file} already exists. It has been renamed to {file}_copy_{str(len(directory_list))}')
        os.makedirs(project_folder)
    project_results_folder = project_folder + 'results/'
    project_data_folder = project_folder + 'data/'
    os.makedirs(project_results_folder)
    os.makedirs(project_results_folder + 'Clusters/')
    os.makedirs(project_results_folder + 'images/')
    os.makedirs(project_results_folder + 'mast/')
    os.makedirs(project_results_folder + 'meme/')
    os.makedirs(project_results_folder + 'meme/motifs/')
    os.makedirs(project_data_folder)
    os.makedirs(project_data_folder + 'Clusters/')
    os.makedirs(project_data_folder + 'Clusters/tables/')
    os.makedirs(project_data_folder + 'sequences/')


# Create cyrest client
cytoscape = cyrest.cyclient()

# Load network to Cytoscape
print('Loading network to Cytoscape')
for file in os.listdir(os.path.join(os.getcwd(), 'network/')):
    network_path = os.path.join(os.getcwd(), 'network', file)
    cytoscape.network.load_file(network_path)
print('Network loaded successfully')

# calculate layout
print('Applying favourite layout')
cytoscape.layout.apply_preferred()
print('Layout applied successfully')

# Default visual style
cytoscape.vizmap.apply(styles='default')

# Run the first rule (run_clusterONE)
#print(f'Running clusterONE algorithm and extracting general information about the full {file} SSN.\n')
os.system('snakemake -s rules/run_clusterONE.rules --cores 1')
#print('Step complete.\n')
#print(f'Results are written to {project_results_folder}.\n')
#print(f'General information about SSN is writen to file general_info.txt\n')

# Run the second rule (to get information per cluster)
#print(f'Extracting information for each defined cluster from the {file} SSN.\n')
os.system('snakemake -s rules/cluster_info.rules --cores 1')
#print('Step complete.\n')
#print(f'Results are written to {project_results_folder}Clusters/\n')

# Run the third rule (to calculate seqeunce motifs)
#print(f'Running MEME from MEME-Suite to calculate sequence motifs for defined clusters from the {file} SSN.\n')
os.system('snakemake -s rules/sequence_motifs.rules --cores 1')
#print('Step complete.\n')
#print(f'Results written to {project_results_folder}meme/\n')

# Run the found rule (search calculated fingerprints against the dataset)
#print(f'Running MAST from MEME-Suite to search previously calculated sequence motifs against the dataset used to calculate the {file} SSN.\n')
os.system('snakemake -s rules/run_mast.rules --cores 1')
#print('Step complete.\n')
#print(f'Results written to {project_results_folder}mast/\n')
