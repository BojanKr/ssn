# ssn

## Description
This repository is a workflow allowing automated analyses of
Sequence Similarity Networks (SSNs) of diverese groups of proteins.

## Installation

Here are the requirements for the workflow to work.

1. You need to Install:
- Cytoscape, is used to read and integrate the networks with any type of attribute data. The version used here is 3.8.0. You can download it here: https://cytoscape.org/
- cyREST, which allows to execute workflow under Cytoscape.
You can install it using Cytoscape:
Open Cytoscape > Go to Apps > Apps Manager > Install Apps > Search for cyREST > Click install
It might be already installed. Then you will see written in the list: cyREST (Installed)

2. Create a folder called /tools under /ssn-master

3. Then under the newly created folder /tools install:
- MEME-Suite, is a motif-based sequence analysis tools.
Download the latest version here: http://meme-suite.org/doc/download.html
Move the unzip version of "meme-5.2.0" under the folder /tools
Follow the instructions for a quick install here: http://meme-suite.org/doc/install.html?man_type=web
- Cluster-one, is a graph clustering algorithm, allowing here to generate
overlapping clusters.
You can download the ".jar" under the section "Download"
here: https://www.paccanarolab.org/cluster-one/
Save the ".jar" file under the /tools folder


## Usage

1. Open Cytoscape, it must be running when you run the workflow.

2. Place the SSN you want to analyse in format '.xgmml' in the folder '/network'

3. Place the seqeunce dataset used to calculate the SSN in folder '/dataset' and name it dataset.txt (the dataaset must be in FASTA format)

4. To execute the workflow use: 
```python 
python main.py
```

5. The workflow will now go through a series of steps which are split into four snakefiles.

After the first snakefile (run_clusterONE.rules) executes, the SSN will be clustered into isofunctional clusters using the ClusterONE algorithm (1). Clusterring information and general information about the SSN will be written to a file (see under Output)

After the second snakefile executes (cluster_info.rules), the defined isofunctional clusters will be separated from the full SSN as subclusters and information for each defined cluster will be written to a file (see under Output)

The third snakefile (run_meme.rules) will run the MEME tool (2) from the MEME-Suite (3) and search for sequence motifs in each cluster defined in previous two rules.
During this step the user will be prompted to enter the number of motifs that should be found and their minimum and maximum width. 
The output of this step are folders that contain sequence logos and text files containing motifs calculated by MEME (see more under Output).

The fourth snakefile (run_mast.rules) will run the MAST tool (4) from the MEME-SUite (3) and based on its results it will make predictions about the possible functions of 
some uncharacterized parts of the SSN. Outputs of this step are text files with MAST search hits and a predictions files which lists all of the 
clusters whose function can be predicted based on the MAST search (see more under Output). The user will be prompted in this step to enter number of a 
defined cluster whose fingerprint should be search against the dataset used to make the SSN. The user will also have to provide the number of motifs 
that will comprise the said fingerprint and names of those motifs (names are numbers assigned by the MEME tool so this information can be found in the
output of the previous step in the workflow).


# NB
You must activate the environment ssn_env.yaml using conda in order for the
workflow to function correctly under the right environment. 
Read more about environment using conda here: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#
```
conda activate ssn_env.yaml
```

## Output

1. When the forkflow starts, a new (root) folder is created named after the file of the network analysed.
You will find two folders here, one with the results and one with the data produced and collected during the execution of the workflow.

2. Output files to look for are writen to /root folder/results/
Output of the first snakefile: 
- /results/general_info.txt (general information about the SSN)
- /results/pdb.txt (a list of all proteins with solved structure found in the SSN)
- /results/clusterONE_results.csv (results of applying the ClusterONE algorithm)
- /results/images/ (visualized informaton about the full SSN)
	- Image 1 shows the topology of an unclsutered SSN
	- Image 2 shows isofunctional clusters found after applying ClusterONE algorithm to SSN
	- Image 3 shows all proteins from SwissProt found in the SSN
	- Image 4 shows taxonomical composition of the SSN on the Kingdom level
	- Image 5 shows position of every protein with sloved structure in the SSN
	
Output of the second snakefile:
- /results/Clusters/ (information for each defined cluster is written here)

Output of the third snakefile:
- /results/meme/ (sequence motifs and their logos calculated by MEME are written here)

Output of the fourth snakefile:
- /results/mast/ (in folders named after clusters are results of MAST search and predictions made based on this search)
	- inside a folder named after the selected cluster are following files:
		mast.txt, mast.html, mast.xm (produced by the MAST tool, these files all contain hits from the MAST search)
		predictions.txt (predictions made based on the MAST search)
		predictions.png (predictions made based on the MAST search visualize on the "pivot" network)		


## Example

In order for you to test that you have a proper setup of the workflow,
you can find an example of a network already in the '/network' folder called:
"47305_case_study_full_ssn.xgmml"

## Refereneces:
1. Nepusz T, Yu H, Paccanaro A. Detecting overlapping protein complexes in protein-protein interaction networks. Nat Methods. 2012 May 18;9(5):471–2. 
2. Bailey L Timothy EC. Fitting a mixture model by expectation maximization to discover motifs in biopolymers. Proc Int Conf Intell Syst Mol Biol. 1994;2:28–36. 
3. Bailey TL, Boden M, Buske FA, Frith M, Grant CE, Clementi L, et al. MEME Suite: Tools for motif discovery and searching. Nucleic Acids Res. 2009;37:202–8.
4. Bailey TL, Gribskov M. Combining evidence using p-values: Application to sequence homology searches. Bioinformatics. 1998;14(1):48–54. 
