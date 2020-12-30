
Here you may find a description of the folder contents and a summary of the integrative modeling procedure used to localize a model of TFIIK in the larger core-Mediator Preinitiation complex.

1. data

The input data consists of the following files:

- topology file - specifies the system's representation; which molecules are to be modeled, which pdb and sequence files may be accessed to represent each domain, and at which level may each component be represented. The listed molecules correspond only to those present in the EMDB 3850 structure (i.e. core Mediator, TFIIA, TFIIB, SPT15/TBP, TFIIF, TFIIE, TFIIH, and RNA Polymerase II on a DNA template) and the TFIIK trimer proteins (Kin28, Ccl1, and Tfb3). 

    topology.txt

- sequence files - sequence files were all obtained from uniprot and formatted slightly for compatibility in IMP

    emd_3850.fasta.txt   (which was separated into the following files for convenience)
        dna.fasta.txt
        gtfs.fasta.txt
        head.fasta.txt
        middle.fasta.txt
        pol_ii.fasta.txt
        tfiih.fasta.txt

    tfiik.fasta.txt

- pdb models

    emd_3850.pdb (equivalent to PDB 5OQM and which was similarly separated into the following files for convenience)
        dna.pdb
        gtfs.pdb
        head.pdb
        middle.pdb
        pol_ii.pdb
        tfiih.pdb

    tfiik.pdb

- em models

    emd_3850.map.mrc (simply renamed the em density retrieved from EMDB 3850)
    emd_3850_tfiik_seg.map.mrc (regions obtained after subjecting emd_3850.map.mrc to map segmentation at threshold 0.05)

- gaussian mixture model files - produced by using IMP's create_gmm.py functionality
    emd_3850_tfiik_seg.map.mrc.gmm.50.mrc (used for visual inspection)
    emd_3850_tfiik_seg.map.mrc.gmm.50.txt (file containing functions used during modeling)

- crosslinking datasets - taken from Robinson 2016 (Cell) and Schilbach 2017 (Nature); selected unique and unambiguous crosslinks from both datasets pertaining only to the proteins specified in the topology file. For convenience, each large dataset was also separated to list only those crosslinks pertaining to TFIIK subunits and may be used in modeling.

    robinson_mpic_xls.csv 
    schilbach_mpic_xls.csv

    robinson_tfiik_xls.csv
    schilbach_tfiik_xls.csv

2. results

- two main folders may be found:
    - clustering: contains results for kmeans clustering of two independent simulations, with representative models and corresponding localization probability density envelope mrc files, along with pdf results showing both a heatmap and a dendrogram representation of the clustering of the top 500 models for each run under subfolder figures.
    - precision: contains raw results for calculation of per-residue root mean squared fluctuation in the top 500 models used for clustering as well as plots of each queried protein's rmsf values under the subfolder figures.

3. chimera_session

- a chimera session file may be opened to visualize the results from the modeling procedure.

4. scripts

    modeling.py - this is the modeling protocol file. It accesses the data folder llocated at "../data" and in the default settings will utilize the segmented and gaussian mixture model envelope derived for TFIIK, the PDB files derived from EMDB 3850, the PDB file for the TFIIK model, the sequence datasets from Uniprot, and the full crosslinking dataset from Robinson et al 2016. Alternatively, the modeling script may be edited to include both datasets or consider the smaller datasets of crosslinks pertaining only to TFIIK for similar modeling results, at reduced computational cost.

    clustering.py - this is the clustering protocol applied to the resulting dataset and is based on a script derived from the RNA Polymerase II tutorial on the IMP website (version 2.12.0). This script will consider the top 500 models of an individual run and cluster them according to the TFIIK trimer subunits.

    prcision_rmsf.py - this is the precision protocol applied to the resulting clusters obtained from clustering.py and is based on a script derived from the RNA Polymerase II tutorial on the IMP website (version 2.12.0). This script will obtain the per-residue root mean square fluctuation values of the TFIIK trimer subunits in a specified cluster.
