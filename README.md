# Hippocampus_cortex_gradient_youth

Code and data for the manuscript entitled "Hippocampus reorganizes its function and geometry along a dual long-axis for cognitive maturation in youth" by Zeng et al.

## File descriptions

1. 

## Installation



## Downloading data



## Original data

The original and preprocessed HCP-D data, after meeting eligibility requirements, can be accessed here: https://humanconnectome.org/study/hcp-lifespan-development.

## Dependencies

Part of the dependency packages (especially for those with modifications) have been stored in the `Dependencies/` folder to ensure version compatibility. 

1. Hippocampal gradient computation and visualization: Matlab toolbox [cifti-matlab](https://github.com/Washington-University/cifti-matlab)(v2.1.0), [gifti](https://github.com/gllmflndn/gifti), [BrainSpace tool](https://github.com/MICA-MNI/BrainSpace)(v0.1.10) with a few modifications in `plotter.m` used by `plot_hemispheres.m` for better brain visualization, [gramm](https://github.com/piermorel/gramm), [SurfStat](https://math.mcgill.ca/keith/surfstat/) with a few modifications to add some useful functions.
2. Hippocampal geometric eigenmodes computation: Partly adapted from [BrainEigenmodes](https://github.com/NSBLab/BrainEigenmodes/tree/main) by Pang et al., need [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation), [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall), [Connectome Workbench](https://www.humanconnectome.org/software/get-connectome-workbench), and [Gmsh](https://gmsh.info/) software, and the python environment is presented in `Dependencies/Python/BrainEigMod_package_list.txt`.
3. Transcriptomic association and gene enrichment analysis of the hippocampal gradients: Partly adapted from [Hippocampus_AP_Axis](https://github.com/illdopejake/Hippocampus_AP_Axis) by Vogel et al. The python environment is presented in `Dependencies/Python/Hipp_Gene_package_list.txt`. [Metascape web tool](www.metascape.org). [Developmental-specific expression analysis (SEA) tool](http://genetics.wustl.edu/jdlab/cseatool-2/)
4. Analysis of developmental effects: R packages [mgcv](https://rdocumentation.org/packages/mgcv/versions/1.8-42)(v1.8-42) and [gratia](https://rdocumentation.org/packages/gratia/versions/0.8.1)(v0.8.1).

## Compatibility

The codes have been tested on Python 3.8, MATLAB R2018a, and R 4.2.3.

## Demo



## Citation

If you use our code in your research, please cite us as follows:

[PREPRINT] D. Zeng, Q. Li, D. Li, ..., Y. He, X. Zuo, S. Li, Dual long-axis reorganization of hippocampus in youth, bioRxiv (2023) (DOI: [10.1101/2023.11.03.565423](https://www.biorxiv.org/content/10.1101/2023.11.03.565423v1.article-metrics))

## Further details

Please contact debin_z@126.com if you need any further details.
