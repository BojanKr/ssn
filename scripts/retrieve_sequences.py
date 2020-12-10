from timeit import default_timer as timer
import requests
import pandas as pd
import os


def get_root_dir():
    for file in os.listdir(os.path.join(os.getcwd(), 'network/')):
        if not file.endswith('xgmml'):
            continue
        else:
            dir = os.path.join(os.getcwd(), file + '/')

            return dir


root_dir = get_root_dir()


directory = os.path.join(os.getcwd(), root_dir, 'data', 'Clusters/')
output_directory = os.path.join(os.getcwd(), root_dir, 'data', 'sequences/')
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

def write_seqeunces(url, directory, file):

    """ Fetches seqeunces from uniprot and writes to a file """

    with requests.get(url, stream=True) as r:
        r.raise_for_status()

        with open(directory+file, 'ab') as fn:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    fn.write(chunk)


def main():

    print('Fetching sequences from UniProt.\n')
    print(f'Sequences written to {output_directory}\n')

    for file in os.listdir(directory):
        if file.endswith('.txt'):
            # Read file into df
            file_df = pd.read_csv(directory+file, header=None)
            # Create a list of accession numbers
            acc_list = file_df[0].tolist()
            # Create a query list (adding accession: on from of every accession number
            query_list = [f'accession:{acc}' for acc in acc_list]
            if len(query_list) < 150:
                # Create a one-string query
                full_query = "+OR+".join(query_list)
                # Create url
                url = f"https://www.uniprot.org/uniprot/?query={full_query}&format=fasta"
                write_seqeunces(url=url, directory=output_directory, file=file)
            else:
                # Cut the large list into pieces
                pieces = [query_list[x:x+150] for x in range(0, len(query_list), 150)]
                for piece in pieces:
                    query = "+OR+".join(piece)
                    url = f"https://www.uniprot.org/uniprot/?query={query}&format=fasta"
                    write_seqeunces(url=url, directory=output_directory, file=file)

    with open(os.path.join(os.getcwd(), root_dir, 'results', 'motifs.txt'), 'w') as f:
        f.write('test')


if __name__ == '__main__':
    main()