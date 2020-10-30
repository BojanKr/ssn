import os

root_dir = os.listdir(os.path.join(os.getcwd(), 'network/'))[0]
data_dir = os.path.join(os.getcwd(), root_dir, 'data', 'ready_edge_table.csv')
results_dir = os.path.join(os.getcwd(), root_dir, 'results', 'clusterONE_results.csv')
script_dir = os.path.join(os.getcwd(), 'scripts', 'clusterONE.py')

shell = "java -jar tools/cluster_one-1.0.jar " + data_dir + " --min-density auto\
 --input-format edge_list --output-format csv --max-overlap 0.8 --merge-method single\
  --penalty 2 --min-size 10 --similarity match > " + results_dir

print('Running clusterONE with the following command:\n')
print(shell)
print('\n')
print(f'To change parameters edit the shell command in {script_dir}')

os.system(shell)