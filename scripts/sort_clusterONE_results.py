import pandas as pd
import os

def main():

    print('Cleaning and reformatting ClusterONE results.\n')

    file = os.listdir(os.path.join(os.getcwd(), 'network'))[0]

    # Read clusterONE results
    cone_df = pd.read_csv(os.path.join(os.getcwd(), file, 'results', 'clusterONE_results.csv'))

    # Keep only cluster numbers and accession numbers of members of clusters
    cone_df = cone_df[['Cluster', 'Members']]

    # Write sorted clusterONE results to a file
    out_file = os.path.join(os.getcwd(), file, 'results', 'sorted_clusterONE_results.csv')
    with open(out_file, 'w') as f:
        f.write('shared name' + ',' + 'clusterONE' + '\n')
    # Sort results
    for index, row in cone_df.iterrows():
        for acc in row['Members'].split():
            with open(out_file, 'a') as f:
                f.write(acc + ',' + str(row['Cluster']) + '\n')

    print(f'Cleaned and sorted clusterONE results written to {out_file}')

if __name__ == '__main__':
    main()

