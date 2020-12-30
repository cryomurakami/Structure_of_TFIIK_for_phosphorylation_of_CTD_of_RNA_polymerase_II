
"""
#############################################
##  Integrative Modeling Script
#############################################
#
# Short modeling script combining EM and Crosslinking data
#           to localize TFIIK on core-Mediator-PIC
#
# Adapted for TFIIK core-Mediator-PIC complex by Jose Gorbea
#
"""
import os
import sys

import IMP
import RMF
import IMP.atom
import IMP.core
import IMP.algebra
import IMP.container
import IMP.rmf

import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.em
import IMP.pmi.restraints.stereochemistry
import ihm.cross_linkers


#-------------------------------------------------------------
# Define Input Files
#-------------------------------------------------------------
datadirectory = "../data/"
topology_file = datadirectory + "topology.txt"
target_gmm_file = datadirectory + "emd_3850_tfiik_seg.map.mrc.gmm.50.txt"

#-------------------------------------------------------------
# Set MC Sampling Parameters
#-------------------------------------------------------------
num_frames = 20000
if '--test' in sys.argv:
    num_frames = 100
num_mc_steps = 10

#-------------------------------------------------------------
# Create movers
#-------------------------------------------------------------

# rigid body movement params
rb_max_trans = 4.00
rb_max_rot = 0.3
# flexible bead movement
bead_max_trans = 4.00
# super rigid body movement
max_srb_trans = 4.0
max_srb_rot = 0.3

#-------------------------------------------------------------
# Build the Model Representation
#-------------------------------------------------------------

# Initialize model
m = IMP.Model()

# Create list of components from topology file
topology = IMP.pmi.topology.TopologyReader(topology_file,
                                           pdb_dir=datadirectory,
                                           fasta_dir=datadirectory
                                           )
domains = topology.get_components()

bs = IMP.pmi.macros.BuildSystem(m)
bs.add_state(topology)

representation, dof = bs.execute_macro(max_rb_trans=rb_max_trans,
                                       max_rb_rot=rb_max_rot,
                                       max_bead_trans=bead_max_trans,
                                       max_srb_trans=max_srb_trans,
                                       max_srb_rot=max_srb_rot)

#-------------------------------------------------------------
# Define Degrees of Freedom
#-------------------------------------------------------------

# Select and gather all particles to fix in place
#  Namely, all subunits except TFIIK components, 
#  but fix Tfb3 region associated to core complex

fixed_particles=[]
for prot in [
             "Rpb1","Rpb2","Rpb3","Rpb4","Rpb5","Rpb6","Rpb7","Rpb8","Rpb9","Rpb10","Rpb11","Rpb12",
             "Med6","Med8","Med11","Med17","Med18","Med20","Med22",
             "Med4","Med7","Med9","Med10","Med14","Med19","Med21","Med31",
             "Spt15","Toa1","Toa2","Sua7","Tfg1","Tfg2","Tfa1","Tfa2",
             "Rad3","Tfb1","Tfb2","Tfb4","Tfb5","Ssl1","Ssl2",
             "NDNA","TDNA"
             ]:
    fixed_particles+=IMP.atom.Selection(representation,molecule=prot).get_selected_particles()
fixed_particles+=IMP.atom.Selection(representation, molecule="Tfb3", residue_indexes=range(8,145)).get_selected_particles()


# Fix corresponding movers using dof
fixed_beads,fixed_rbs=dof.disable_movers(fixed_particles,
                                         [IMP.core.RigidBodyMover,
                                          IMP.pmi.TransformMover])


# Shuffling for random initial conformations 
IMP.pmi.tools.shuffle_configuration(representation,
                                    excluded_rigid_bodies=fixed_rbs,
                                    max_translation=20,
                                    verbose=False,
                                    cutoff=5.0,
                                    niterations=100)

# Add default mover parameters to simulation
outputobjects = []  # reporter objects (for stat files)
sampleobjects = []  # sampling objects

#-------------------------------------------------------------
# Define Scoring Function Components
#-------------------------------------------------------------

# Here we are defining a number of restraints on our system.
# For all of them we call add_to_model() so they are incorporated into scoring
# We also add them to the outputobjects list, so they are reported in stat
# files

#-------------------------------------------------------------
# Connectivity Restraints
#-------------------------------------------------------------

moldict = bs.get_molecules()[0]
for molname in moldict:
    for mol in moldict[molname]:
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            mol, scale=2.0, label=molname)
        cr.add_to_model()
        cr.set_label(molname)

        outputobjects.append(cr)
        sampleobjects.append(cr)

#-------------------------------------------------------------
# Excluded Volume (evaluated at 20 residue resolution)
#-------------------------------------------------------------

ev1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=representation,
                                                              resolution=10)
ev1.set_label('Excluded Volume')
ev1.add_to_model()
outputobjects.append(ev1)

#-------------------------------------------------------------
# XL-MS Datasets:
#-------------------------------------------------------------

# To use this restraint we have to first define the data format
# Here assuming that it's a CSV file with column names that may need to change

# Other options include the linker length and the slope for nudging
# components together

kw = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
# kw.set_unique_id_key("id")
kw.set_protein1_key("prot1")
kw.set_protein2_key("prot2")
kw.set_residue1_key("res1")
kw.set_residue2_key("res2")
# kw.set_id_score_key(None)

# XL-MS Dataset 1

Robinson = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
Robinson.create_set_from_file(datadirectory + 'robinson_mpic_xls.csv')
# Robinson.create_set_from_file(datadirectory + 'robinson_tfiik_xls.csv')

x_Robinson = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=representation,
    label="Robinson",

    database=Robinson,
    length=30,
    resolution=1.0,

    slope=0.02,
    weight=1.0,
    )

x_Robinson.add_to_model()
sampleobjects.append(x_Robinson)
outputobjects.append(x_Robinson)
dof.get_nuisances_from_restraint(x_Robinson)

# # XL-MS Dataset 2

# Schilbach = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
# Schilbach.create_set_from_file(datadirectory + 'schilbach_mpic_xls.csv')
# Schilbach.create_set_from_file(datadirectory + 'schilbach_tfiik_xls.csv')

# x_Schilbach = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
#     root_hier=representation,
#     label="Schilbach",
#
#     database=Schilbach,
#     length=30,
#     resolution=1.0,
#
#     slope=0.02,
#     weight=1.0,
#     )
#
# x_Schilbach.add_to_model()
# sampleobjects.append(x_Schilbach)
# outputobjects.append(x_Schilbach)
# dof.get_nuisances_from_restraint(x_Schilbach)


#----------------------------------------------------------------------
# Cryo-EM Model
#----------------------------------------------------------------------

"""
# EM Model Description

Original Cryo-EM Model Contains:
    
    core Mediator (Head and Middle modules, with a truncated Med14 stalk)
    12 subunit RNA Pol II
    TBP, TFIIB, TFIIA, TFIIE, TFIIF
    complete TFIIH

    DOI 10.1038/nature24282
    EMDB 3850

Model Processed to facilitate localization of TFIIK:

    Segmented Cryo-EM Model on UCSF-Chimera.
    Extracted low-resolution TFIIK Density based on descriptions on
    work by Robinson et al (2016; /10.1016/j.cell.2016.08.050) and Schilbach et al
    (2017; 10.1038/nature24282).


"""

# Electron Microscopy Restraint

em_components = IMP.pmi.tools.get_densities(representation)

gemt = IMP.pmi.restraints.em.GaussianEMRestraint(em_components,
                                                  target_gmm_file,
                                                  slope=0.000001,
                                                  weight=100.0)
gemt.add_to_model()
outputobjects.append(gemt)

#--------------------------
# Monte-Carlo Sampling
#--------------------------

# This object defines all components to be sampled as well as the sampling protocol

mc1 = IMP.pmi.macros.ReplicaExchange0(m,
                                      root_hier=representation,
                                      monte_carlo_sample_objects=dof.get_movers(),
                                      output_objects=outputobjects,
                                      crosslink_restraints=sampleobjects,
                                      monte_carlo_temperature=1.0,
                                      
                                      simulated_annealing=True,
                                      simulated_annealing_minimum_temperature=1.0,
                                      simulated_annealing_maximum_temperature=1.5,
                                      simulated_annealing_minimum_temperature_nframes=200,
                                      simulated_annealing_maximum_temperature_nframes=20,
                                      
                                      replica_exchange_minimum_temperature=1.0,
                                      replica_exchange_maximum_temperature=2.5,
                                      
                                      number_of_best_scoring_models=100,
                                      monte_carlo_steps=num_mc_steps,
                                      number_of_frames=num_frames,
                                      global_output_directory="output")

# Start Sampling
mc1.execute_macro()
