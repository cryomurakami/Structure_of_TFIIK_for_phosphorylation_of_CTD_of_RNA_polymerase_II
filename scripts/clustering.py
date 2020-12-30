import matplotlib as mpl
mpl.use('Agg')

import IMP
import IMP.pmi
import IMP.pmi.macros
import sys,os

# most common settings
num_top_models = 500
run_number = 1

merge_directories = ["../run"+str(run_number)+"/"]
prefiltervalue = 2000
out_dir = "kmeans_r%i_m%i/" %(run_number,num_top_models)

#################################
# should not have to change below
##################################

model=IMP.Model()

# initialize the macro
mc=IMP.pmi.macros.AnalysisReplicaExchange0(model,
                                           merge_directories=merge_directories)

# fields that have to be extracted for the stat file
feature_list=["ISDCrossLinkMS_Distance_intrarb",
              "ISDCrossLinkMS_Distance_interrb",
              "ISDCrossLinkMS_Data_Score",
              "GaussianEMRestraint_None",
              "SimplifiedModel_Linker_Score_None",
              "ISDCrossLinkMS_Psi",
              "ISDCrossLinkMS_Sigma"]

# Dictionary of densities to be calculated
# the key is the name of the file and the value is the selection
# example: {"med17-CTD":[(200,300,"med17")],"med17-CTD.med14":[(200,300,"med17"),"med14"]   }
density_names = {"Ccl1":["Ccl1"],
                 "Kin28":["Kin28"],
                 "Tfb3":["Tfb3"],}

# list of component names needed to calculate the RMSD for the clustering
rmsd_names = {"Ccl1":"Ccl1",
              "Kin28":"Kin28",
              "Tfb3":"Tfb3"}

# components used for structural alignment
align_names = None # (None because EM provides reference frame)

mc.clustering(prefiltervalue=prefiltervalue,                   # prefilter the models by score
              number_of_best_scoring_models=num_top_models,    # number of models to be clustered
              alignment_components=None,                       # list of proteins you want to use for structural alignment
              rmsd_calculation_components=rmsd_names,          # list of proteins used to calculated the rmsd

              distance_matrix_file="distance.rawmatrix.pkl",   # save the distance matrix
              outputdir=out_dir,                               # location for clustering results
              feature_keys=feature_list,                       # extract these fields from the stat file
              load_distance_matrix_file=False,                 # skip the matrix calculation and read the precalculated matrix
              display_plot=True,                               # display the heat map plot of the distance matrix
              exit_after_display=False,                        # exit after having displayed the distance matrix plot
              get_every=1,                                     # skip structures for faster computation
              number_of_clusters=num_clusters,                 # number of clusters to be used by kmeans algorithm
              voxel_size=3.0,                                  # voxel size of the mrc files
              density_custom_ranges=density_names)             # setup the list of densities to be calculated
