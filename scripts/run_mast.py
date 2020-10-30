import os

root_dir = os.listdir(os.path.join(os.getcwd(), 'network/'))[0]

def select_cluster_numbers():

    """ Creates a list of cluster numbers that should be searched against
        a sequence dataset of interest """

    cluster_list = []

    while True:
        cluster_number = input('Specify clusters whose motifs should be searched against the sequence dataset (press Enter when done): ')
        if cluster_number == '':
            break
        else:
            cluster_list.append(cluster_number)

    return cluster_list


def main():

    mast_dir = os.path.join(os.getcwd(), 'tools', 'meme-5.2.0', 'src', 'mast')
    dataset = os.path.join(os.getcwd(), 'dataset', 'dataset.txt')
    motif_directory = os.path.join(os.getcwd(), root_dir, 'results', 'meme', 'motifs/')
    output_directory = os.path.join(os.getcwd(), root_dir, 'results', 'mast/')
    clusters = select_cluster_numbers()

    for cluster in clusters:
        motif_dir = motif_directory + 'Cluster_' + str(cluster) + '/' + 'meme.txt'
        num_of_motifs = input('How many motifs to scan against the dataset? ')
        # Create a list of motifs (specified by number as assigned by meme) that comprise the fingerprint
        # which will be searched against the seqeunce dataset
        fingerprint = []
        for i in range(int(num_of_motifs)):
            motif = input('Enter motif numbers (as asssigned by MEME): ')
            fingerprint.append('-mi ' + motif)
        output_dir = output_directory + 'Cluster_' + str(cluster) + '/'
        full_fingerprint = ' '.join(fingerprint)
        print(full_fingerprint)
        mast = f'{mast_dir} -o {output_dir} {full_fingerprint} {motif_dir} {dataset}'
        print(mast)

        with open(os.path.join(os.getcwd(), root_dir, 'data', 'used_clusters.txt'), 'a') as f:
            f.write(str(cluster))

        os.system(mast)


if __name__ == '__main__':
    main()