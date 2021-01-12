# Eeg_network_pipeline
Functions for functional connectivity estimation, network discovery and comparison in EEG data.

We collect in this repo the Matlab scripts/functions used for network-based analyses of EEG data in the **Sound and Speech Perception Research Group** at TTK, Budapest (PI: Istvan Winkler). Development is ongoing, changes might occur at any time and at any level.

### Dependencies / compatilibity:
#### Matlab
- Development is primarily with Matlab 2017a or higher, Octave compatibility is rarely - if ever - tested. 
- Some functions (mostly the surrogate data generator wrappers found in `/edgePruning`) use `parfor` for across-subject loops, requiring Parallel Computing Toolbox for optimal performance. Note that without the Parallel Computing Toolbox `parfor` loops work essentially as for loops but with a different index order. 
- Some functions rely on Statistics and Machine Learning Toolbox (occasional calls to `normcdf`, `normpdf`, `kstest`, etc.). It is a (very) long-term goal to not rely on the Stat. Toolbox.
- Occasionally we might rely on Signal Processing Toolbox, hopefully only for stopgap measures. 
#### Toolboxes for modularity
- `/python_leidenalg` is a very simple interface for using the great [leidenalg python package](https://github.com/vtraag/leidenalg) (note that the leiden algorithm has been included in [igraph](https://igraph.org/redirect.html))
- Functions under `/modularity` are mainly wrappers for calling the `genlouvain` and `iterated_genlouvain` functions from the great [Genlouvain toolbox](https://github.com/GenLouvain/GenLouvain)

### Functions from other collections
- We pulled functions (ones we rely on) from the great [Network Community Toolbox](http://commdetect.weebly.com/):
<br>`/modularity/consensus_iterative`, `/modularity/consensus_similarity`, `/modularity/multislice_static_unsigned`, `/modularity/zrand.m`
<br> We felt free to do so since there is no mention of a specific license / any restrictions on their website - we treat these as if under MIT license. Please cite the [Network Community Toolbox](http://commdetect.weebly.com/) and the corresponding papers (see on their website and in the helps of the functions) whenever using their work. 
- We also pulled functions (ones we rely on) from the great [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/):
<br>`/measures/betweenness_wei.m`, `/measures/charpath.m`, `/measures/clustering_coef_wu.m`, `/measures/distance_wei.m`, `/measures/edge_betweenness_wei.m`, `/measures/efficiency_wei.m`, `/measures/module_degree_zscore.m`, `/measures/participation_coef.m`, `/measures/rich_club_wu.m`, `/measures/strengths_und.m`, `/measures/weight_conversion.m`, 
<br>`/modularity/null_model_und_sign_mod`, `/modularity/randmio_und_signed_mod`
<br> As with the Network Community Toolbox, we felt free to do so since there is no mention of a specific license / any restrictions on their website - we treat these functions as if under MIT license. Please cite the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/) and the corresponding papers (see on their website and in the helps of the functions) whenever using their work. 

### EEG data format
We assume that data is already preprocessed (e.g. re-referenced, filtered for muscle and other artefacts, bandpass filtered to ranges of interest, source-reconstructed, averaged into ROIs based on a parcellation). We expect preprocessed data in 3D/4D arrays, with dimensions `[ROIs/channels X samples X epochs (X conditions)]`, in one .mat file per subject.

### Main analysis steps and corresponding functions

1. __Connectivity estimation (`/connectivity`):__
<br> Functions under `/connectivity` are used for connectivity estimation on epoch-level data (matrices with dimensions `[ROIs/channels X samples]`). 
<br> Currently we only support undirected metrics: PLV, iPLV, PLI, wPLI, ampCorr, orthAmpCorr. Their functions have the same name.
<br> Further functions:
<br>  `connectivitySpeedtest`   - Helper function to estimate the speed of per-epoch connectivity estimations. Especially helpful before surrogate data testing. 
<br>  `connectivityWrapperReal` - Wrapper for connectivity estimation on a whole dataset. 

2. __Edge pruning and other network discovery/construction methods:__

3. Clustering / modularity / other network structure discovery methods;

4. Network similarity measures / network comparisons.

Functions are grouped mostly according to the above classification.

We pulled a few functions from the great [Network Community Toolbox](http://commdetect.weebly.com/): 
<br>`/modularity/zrand.m`
<br>`/modularity/consensus_similarity`


