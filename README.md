# Eeg_network_pipeline
Functions for network discovery and comparison in EEG data.

The aim of the repo is to collect the Matlab scripts/functions used for network-based analyses of EEG data in the **Sound and Speech Perception Research Group** at TTK, Budapest (PI: Istvan Winkler). We assume that data is already preprocessed (e.g. re-referenced, filtered for muscle and other artefacts, bandpass filtered to ranges of interest, reconstructed into source-space). Main analysis steps include: 

1. Connectivity measurements;

2. Edge pruning and other network discovery/construction methods;

3. Clustering / modularity / other network structure discovery methods;

4. Network similarity measures / network comparisons.

Functions are grouped mostly according to the above classification.

We pulled a few functions from the great [Network Community Toolbox](http://commdetect.weebly.com/): 
<br>`/modularity/zrand.m`
<br>`/modularity/consensus_similarity`

Functions under `/modularity` are mainly wrappers for calling the `genlouvain` and `iterated_genlouvain` functions from the great [Genlouvain toolbox](https://github.com/GenLouvain/GenLouvain)
The folder `/python_leidenalg` is a very simple interface for using the great [leidenalg python package](https://github.com/vtraag/leidenalg) (note that the leiden algorithm has been included in [igraph](https://igraph.org/redirect.html)).
