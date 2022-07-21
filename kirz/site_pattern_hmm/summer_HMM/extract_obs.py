import sys
import gzip
import numpy as np
import pandas as pd
import time


# Define a function to split a genotype matrix into non-overlapping windows.
def genotype_matrix_windows(
        variant_positions,
        polarized_genotype_matrix,
        window_size=50_000,
        sequence_length=20_000_000,
):
    # Intialize a dicctionary with the start and stop position for each window.
    windows = {}
    index = 1
    for window_start in range(0, int(sequence_length), int(window_size)):
        windows[index] = [window_start, (window_start + window_size)]
        index += 1
    # Locate what window each variant is in.
    index = 0
    pos = variant_positions[index]
    for key in windows:
        start, stop = windows[key]
        while start <= pos < stop:
            windows[key].append(index)
            index += 1
            if index < len(variant_positions):
                pos = variant_positions[index]
            else:
                break
    return windows


# Extracts the observed sequence (binned)
# Variable positions corresponds to the first index in the genotype matrix
# Different from the intro positions, which are just the start and stop positions of all introgressed segements
# Input: var_pos - filepath to the array of all variable positions (filename rep_id_{REP}_var_pos.csv.gz)
#        pol_gen_mat - filepath to the polarized genotype matrix (filename rep_id_{REP}_polarized_geno_mat.csv.gz)
# Output: the observed sequence and a dictionary containing information about it
def extract_O(variable_positions, polarized_genotype_matrix):

    var_pos = variable_positions
    pol_geno_mat = polarized_genotype_matrix


    # TODO: OLD, HARDCODED INSTRUCTIONS
    # Load the mutated tree sequence.
    # rep_id_1_mts = tskit.load('./cs282_sim_data/rep_id_1_mut_tree_seq.ts')
    # Load the variable positions.
    # rep_id_1_var_pos = np.loadtxt('../cs282_sim_data/rep_id_1_var_pos.csv.gz', delimiter=',')
    # Load the genotype matrix.
    # rep_id_1_polarized_geno_mat = np.loadtxt('../cs282_sim_data/rep_id_1_geno_mat.csv.gz', dtype=int,delimiter=',')
    # Load the introgressed region dataframe.
    # rep_id_1_intro_pos_df = pd.read_csv('../cs282_sim_data/rep_id_1_intro_pos.csv.gz', float_precision='round_trip')
    # Inspect the tree-sequence summary.
    # rep_id_1_mts


    # TODO: NON-HARDCODED
    var_pos = np.loadtxt(var_pos, delimiter=',')
    # Load the genotype matrix.
    pol_geno_mat = np.loadtxt(pol_geno_mat, dtype=int, delimiter=',')
    # Load the introgressed region dataframe.
    true_intro_pos = pd.read_csv('../cs282_sim_data/rep_id_1_intro_pos.csv.gz', float_precision='round_trip')




    # Typically 50k long
    # print(rep_id_1_var_pos)
    # print(rep_id_1_polarized_geno_mat)

    # Indexed from 1 - 400
    Windows = genotype_matrix_windows(var_pos, pol_geno_mat, window_size=500)

    #print(Windows[1])
    #print(len(Windows[1])-2)
    #print('This should be smaller than 50,000 {0}'.format(rep_id_1_var_pos[1352]))
    #print('This should be larger than 50,000 {0}'.format(rep_id_1_var_pos[1353]))
    #print(Windows[1][2])
    #print(Windows[1][-1])
    #wind_1_idx = np.asarray(Windows[400][2:], dtype=np.int32)
    #wind_1_geno_mat = rep_id_1_polarized_geno_mat[wind_1_idx,:]
    #print(wind_1_idx)
    #print(wind_1_idx.shape)
    #print(wind_1_geno_mat)
    #print(wind_1_geno_mat.shape)



    # Intialize observed sequence.
    obs_seq = []
    # Define what C, a pattern consistent with introgression, would look like.
    c_pattern = np.array([0, 0, 1, 1])
    # Intialize the start time.
    start = time.time()
    # Iterate through all the windows by key.
    for key in Windows:
        # Extract the values for the window key.
        window_vals = Windows[key]
        # Print the tracker for me.
        print('there are {0} variants in window {1}'.format(len(window_vals), key))
        # If there are variants in that window. Does this mean a window with a single 'C' in it gets left out?
        if len(window_vals) > 2:
            # Extract variable positions in that window. [2:] excludes start pos and end pos
            variants = np.asarray(window_vals[2:], dtype=np.int32)
            # Subset the genotype matrix for that window.
            window_geno_mat = pol_geno_mat[variants,:]
            print(window_geno_mat)
            # Define what C matrix would look like given an arbitrary number of variants.
            c_mat = np.tile(c_pattern, (window_geno_mat.shape[0], 1))
            # If the C matrix is equal to the windowed matrix declare it consistent.
            if np.array_equal(c_mat, window_geno_mat):
                print('C')
                obs_seq.append('C')
            # Else declare the window non-cosistent.
            else:
                print('N')
                obs_seq.append('N')
        # If there are no variants in the window declare in non-consistent.
        else:
            print('N')
            obs_seq.append('N')

    # Intialize the end time.
    end = time.time()
    ## Convert the observation sequence list to an array.
    obs_seq_array = np.asarray(obs_seq)
    print('there are {0} many consistent observations'.format(np.count_nonzero(obs_seq_array == 'C')))
    print('the consistent observations occur in window(s) {0}'.format(np.where(obs_seq_array == 'C')))
    print('the run time for generating one observed sequence is {0} seconds'.format((end-start)/float(60)))



    # throw vp and pgm into genotype_matrix_windows

    O = obs_seq_array
    Bin_dict = {}
    return O, Bin_dict


extract_O(sys.argv[1], sys.argv[2])
