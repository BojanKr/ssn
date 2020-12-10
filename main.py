from py2cytoscape import cyrest
import os

# Create folders
for file in os.listdir(os.path.join(os.getcwd(), 'network/')):
    if not file.endswith('xgmml'):
        continue
    else:
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
    if not file.endswith('xgmml'):
        continue
    else:
        network_path = os.path.join(os.getcwd(), 'network', file)
        try:
            cytoscape.network.load_file(network_path)
            print('Network loaded successfully')
        except:
            continue

# Create view and apply layout
views = cytoscape.view.list()
if not views['views']:
    print('Creating view')
    cytoscape.view.create()
    print('Applying favourite layout')
    cytoscape.layout.apply_preferred()
else:
    print('Applying favourite layout')
    cytoscape.layout.apply_preferred()

# Default visual style
cytoscape.vizmap.apply(styles='default')

# Run the first rule (run_clusterONE)
print(f'Running clusterONE algorithm and extracting general information about the full SSN.\n')
os.system('snakemake -s rules/run_clusterONE.rules --cores 1')
print('Step complete.\n')


# Run the second rule (to get information per cluster)
print(f'Extracting information for each defined cluster from the SSN.\n')
os.system('snakemake -s rules/cluster_info.rules --cores 1')
print('Step complete.\n')

# Run the third rule (to calculate seqeunce motifs)
print(f'Running MEME from MEME-Suite to calculate sequence motifs for defined clusters from the SSN.\n')
os.system('snakemake -s rules/sequence_motifs.rules --cores 1')
print('Step complete.\n')

# Run the found rule (search calculated fingerprints against the dataset)
print(f'Running MAST from MEME-Suite to search previously calculated sequence motifs against the dataset used to calculate the SSN.\n')
os.system('snakemake -s rules/run_mast.rules --cores 1')
print('Step complete.\n')
