import lomap
import sys
import networkx as nx
import numpy as np
import csv
import sys

##################################################################################################
# This example file walks through how to read in molecules with Lomap using an interactive
# clustering mode. To select defaults during clustering, simply hit enter or type the displayed
# options. The user selected clusters will then be sent for design optimization. 
##################################################################################################

#-------------------------------------------------------#
# Generate similarity scores.
#-------------------------------------------------------#
# Read molecules from test directory.
db_mol = lomap.DBMolecules('../test/radial/', output=True, radial=True)
    
# Generate the strict and loose symmetric similarity score matrices.
strict, loose = db_mol.build_matrices()
    
# Convert the similarity matrix to numpy array
sim_np = strict.to_numpy_2D_array()

# Clean data if Lomap produces rare error. If score is NaN, replace with 0.0
n_arr = np.where(np.isnan(sim_np), 0.0, sim_np)

#-------------------------------------------------------#
# Clustering.
#-------------------------------------------------------#
# Generate distance array.
X = 1 - n_arr

# Create ID_list from dbmols prior to clustering. To be incorporated.
mol_names = ['ID']
for i in range(n_arr.shape[1]):
    fname = db_mol[i].getName()
    fsplit = fname.split(".", 1)[0]
    mol_names.append(fsplit)
# Convert the name list to an np array.
mol_arr_pre = np.array(list(mol_names))
# Cleave the 'ID' entry
ID_list = mol_arr_pre[1:]


# Perform clustering.
labels, cluster_kwarg = lomap.cluster_interactive(X, ID_list)

# Generate sub-arrays of clusters. Stored in dictionaries.
sub_arr, sub_ID = lomap.sub_arrays(labels, n_arr, ID_list)

#-------------------------------------------------------#
# Optimization.
#-------------------------------------------------------#
# Example reference ligands.
ref_ligs = ['ejm_31']

# Send the user selected clusters for optimization.
lomap.clusters2optimize(sub_arr, sub_ID, clusters2optim = cluster_kwarg, ref_ligs=ref_ligs)
