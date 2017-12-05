# import psycopg2

import sys
import numpy as np
import _pickle as pickle
import os
import argparse


# TODO let calculations run en masse, and print to stats file accordingly. make choosing the cluster directory extensible
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-k_base', type=int, required=True, help='base k value for running any calculations.')
    parser.add_argument('-k_inc', help='max k value increment. example if k_from is 6, and inc is 2, '
                                       'you should expect 6,7')
    parser.add_argument('-env', choices=['nick', 'npr', 'pbs', 'lsb456'], required=True,
                        help='The environment/machine you\'re working on. Helpful for filepaths.')
    parser.add_argument('-m', choices=['matrix', 'stats'], help='Which programs you would like to run')
    parser.add_argument()
    if not len(sys.argv) > 1:
        print("no arguments specified. Refer to -h or --help.")
        exit(0)
    args = parser.parse_args()
    environment = args.env
    mode = args.m
    base_k = args.k_base
    to_k = base_k + 1
    if args.k_inc:
        to_k = args.k_inc + base_k
    base_path = ""
    matrix_output_path = ""
    if environment == "nick":
        base_path = "/Users/nickpredey/Documents/Networks/PickleFilesMinHash/"
        matrix_output_path = "NOT YET IMPLEMENTED"
        results_output_path = base_path + "test_matrix_40.csv"
    elif environment == "npr":
        base_path = "/home/catherine/Networks_Nick_NPR/"
        results_output_path = base_path + "test_matrix_40.csv"
        matrix_output_path = "/data/matrix_k6.txt"
        # matrix_output_path = "/media/catherine/ExtraDrive1/Network_Matrices/" #nprito, I believe
    elif environment == "pbs":
        base_path = "/home/catherine/Networks_Nick/"
        matrix_output_path = "/media/catherine/My Book/Network_Matrices/"
        results_output_path = base_path + "test_matrix_40.csv"
    elif environment == "lsb456":
        base_path = "/home/lsb456/Networks_nick/"
        results_output_path = base_path + "test_matrix_40.csv"
        matrix_output_path = "/media/CP_MyBook/Pickle_Matrices/"

    if not os.path.exists(base_path):
        print("The path %s  does not exist on this machine. Are you in the right environment?" % basePath)

    pickle_file_path = base_path + "all_sequences.p"
    accession_pickle_file_path = base_path + "accession_dict.p"
    # results_output_path = basePath + "test_matrix.csv"
    clusters_input_path = base_path + "USearch_AA_Centroids/clusters_40"  # This is only on NPR for right now

    if mode == "matrix":
        sequence_dict = pickle.load(open(pickle_file_path, "rb"))
        print_to_matrices(sequence_dict, matrix_output_path, base_k, to_k)
    elif mode == "stats":
        if not os.path.exists(clusters_input_path):
            print("cluster input path %s does not exist. exiting.." % clusters_input_path)
            exit(0)
        accession_dict = pickle.load(open(accession_pickle_file_path, "rb"))
        print("Reading matrix...")
        matrix_k = np.loadtxt(matrix_output_path)
        print("Matrix reading finished")
        calc_minHash_stats(pickle_file_path, accession_dict, matrix_k, clusters_input_path, results_output_path)
    else:
        print("No mode specified.")

main()
