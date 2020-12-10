import os

print('Runing MEME from MEME-Suite to calculate seqeunce motfs for each of the defiend clusters.')

def get_root_dir():
    for file in os.listdir(os.path.join(os.getcwd(), 'network/')):
        if not file.endswith('xgmml'):
            continue
        else:
            dir = os.path.join(os.getcwd(), file + '/')

            return dir

root_dir = get_root_dir()

directory = os.path.join(os.getcwd(), root_dir, 'data', 'sequences/')
path_to_meme = os.path.join(os.getcwd(), 'tools', 'meme-5.2.0', 'src', 'meme')

# Provide the nu,ber of motifs to be searched and their width
number_of_motifs = input('How many motifs should be calculated? ')
min_width = input('Min length of discovered motifs (in number of residues): ')
max_width = input('Max length of discovered motifs (in number of residues): ')

# Run meme for each cluster
for cluster in os.listdir(directory):
    if cluster.endswith('.txt'):
        output_directory = os.path.join(os.getcwd(), root_dir, 'results', 'meme', 'motifs/') + cluster[:-4] + '/'
        meme = f'{path_to_meme} {directory+cluster} -o {output_directory} -protein -minw {min_width} -maxw {max_width} -nmotifs {number_of_motifs}'

        print('Running MEME from MEME-Suite wih the following command:\n')
        print(meme)
        os.system(meme)

with open(os.path.join(os.getcwd(), root_dir, 'results', 'done.txt'), 'w') as f:
    f.write('test')