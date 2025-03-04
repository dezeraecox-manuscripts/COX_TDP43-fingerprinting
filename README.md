[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14960397.svg)](https://doi.org/10.5281/zenodo.14960397)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14965410.svg)](https://doi.org/10.5281/zenodo.14965410)


# COX_TDP43-fingerprinting

This repository contains the analysis code associated with the MND<sub>TDP</sub> single-molecul pulldown project, led by Dr Dezerae Cox. This manuscript has been submitted for publication under the title ***Fingerprinting disease-derived protein aggregates reveals unique signature of Motor Neuron Disease***.

This manuscript has been submitted as a preprint via BioRxiv [here](biorxiv/link). A link to the final version will be provided upon publication.

## Prerequisites

This analysis assumes a standard installation of Python3 (=> 3.8). For specific package requirements, see the ```environment.yml``` file, or  create a new conda environment containing all packages by running ```conda create -f environment.yml```. 

## Raw data

The proteomics ```.RAW``` files have been deposited via the [PRIDE][1]<sup>[1]</sup> partner repository to the [ProteomeXchange Consortium][2]<sup>[2]</sup> under the dataset identifier PXD061429. In addition, preprocessed data (including identification and quantitation for proteomics, and property summaries for microscopy data; hereon termed raw data) have also been uploaded as an open-access [Zenodo dataset](https://doi.org/10.5281/zenodo.14960397). These data can be collected automatically using the ```raw_data.py``` script in each of the respective ```src``` analysis folders.


## Workflow

Raw data were preprocessed using either [MaxQuant][3]<sup>[3]</sup> or [smma][4]<sup>[4]</sup>. Individual scripts to produce each visualisation from the resulting summary data are presented within the ```src``` folder, alongside some utilities scripts for statistical testing, imputing or visualisation.

## References

[1]: https://www.ebi.ac.uk/pride/archive/

1. Perez-Riverol, Yasset, Attila Csordas, Jingwen Bai, Manuel Bernal-Llinares, Suresh Hewapathirana, Deepti J Kundu, Avinash Inuganti, et al. “The PRIDE Database and Related Tools and Resources in 2019: Improving Support for Quantification Data.” Nucleic Acids Research 47, no. D1 (January 8, 2019): D442–50. https://doi.org/10.1093/nar/gky1106.

[2]: http://proteomecentral.proteomexchange.org

2. Deutsch, Eric W., Attila Csordas, Zhi Sun, Andrew Jarnuczak, Yasset Perez-Riverol, Tobias Ternent, David S. Campbell, et al. “The ProteomeXchange Consortium in 2017: Supporting the Cultural Change in Proteomics Public Data Deposition.” Nucleic Acids Research 45, no. Database issue (January 4, 2017): D1100–1106. https://doi.org/10.1093/nar/gkw936.

[3]: https://www.maxquant.org/

3. Tyanova, Stefka, Tikira Temu, and Juergen Cox. “The MaxQuant Computational Platform for Mass Spectrometry-Based Shotgun Proteomics.” Nature Protocols 11, no. 12 (December 2016): 2301–19. https://doi.org/10.1038/nprot.2016.136.

[4]: https://github.com/dezeraecox/smma

4. “SMMA: Single Molecule Microscopy Analysis (SMMA) Module Containing Analysis Functions for Single-Molecule Data”. Last accessed February 24, 2025. https://github.com/dezeraecox/smma.
