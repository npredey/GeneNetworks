# import psycopg2
import csv

import sys
from datasketch import MinHash
import numpy as np
import _pickle as pickle
import pprint
import os
from Bio import SeqIO
import itertools
from sortedcontainers import SortedSet
import re
import argparse


def get_minhash_pickle_filename(k, bp):
    return bp + "accession_jaccard_k" + str(k) + ".p"


def count_kmers(read, k):
    counts = {}
    num_kmers = len(read) - k + 1
    for i in range(num_kmers):
        kmer = read[i:i + k]
        if kmer not in counts:
            counts[kmer] = 0
        counts[kmer] += 1
    return counts


def get_jaccard_index(sequence1, sequence2, k):
    seq1_minHash, seq2_minHash = MinHash(), MinHash()
    seq1_kmers = count_kmers(sequence1, k)
    seq2_kmers = count_kmers(sequence2, k)
    seq1_keys = list(seq1_kmers.keys())
    seq2_keys = list(seq2_kmers.keys())
    for key in seq1_keys:
        seq1_minHash.update(key.encode('utf8'))
    for key in seq2_keys:
        seq2_minHash.update(key.encode('utf8'))
    return seq1_minHash.jaccard(seq2_minHash)


def get_actual_jaccard_index(kmers1, kmers2):
    actual_jaccard = float(len(kmers1.intersection(kmers2))) / float(len(kmers1.union(kmers2)))
    return actual_jaccard


def create_matrices(from_i, to_i, matrix_length):
    matrix_dictionary = dict()
    for i in range(from_i, to_i):
        jaccard_matrix = np.empty(shape=(matrix_length, matrix_length), dtype=np.float32)
        for j in range(1, matrix_length):
            jaccard_matrix[j][0] = j
            jaccard_matrix[0][j] = j
        matrix_dictionary[i] = jaccard_matrix
    return matrix_dictionary


def print_sequences_to_pickle():
    conn_string = "host='localhost' dbname='networks' user='nickpredey' password='postgres'"
    conn = psycopg2.connect(conn_string)
    cursor = conn.cursor()
    cursor.execute("SELECT sequence_id, accession_number, seq FROM sequences")
    records = cursor.fetchall()
    numSequences = len(records)
    sequence_dictionary = dict()
    for record in records:
        sequence_id = record[0]
        seq_and_accession = (record[1], record[2])
        sequence_dictionary[sequence_id] = seq_and_accession
    filepath = "/Users/nickpredey/Documents/Networks/PickleFilesMinHash/all_sequences.p"
    pickle.dump(sequence_dictionary, open(filepath, 'wb'))


def print_minhash_to_pickle(from_k, to_k, sequence_dict):
    for kk in range(from_k, to_k):
        print(kk)
        accession_number_dict = dict()
        for key, value in sequence_dict.items():
            minHash = MinHash()
            kmer = count_kmers(value[1], kk)
            for kmer_key in kmer.keys():
                minHash.update(kmer_key.encode('utf8'))
            accession_number_dict[key] = minHash
        filepath = get_minhash_pickle_filename(kk, basePath)
        pickle.dump(accession_number_dict, open(filepath, 'wb'))


# This takes in the directory of cluster files and returns all of the absolute file
# paths to be used later in the program
def get_cluster_filenames_from_directory(clusterFileDirectory):
    file_paths = []
    for folder, subs, files in os.walk(clusterFileDirectory):
        for filename in files:
            file_paths.append(os.path.abspath(os.path.join(folder, filename)))
    return file_paths


def get_and_print_accessions_to_pickle(sequence_dictionary):
    accession_dict = dict()
    for key, value in sequence_dictionary.items():
        accession = value[0]
        sequence = value[1]
        accession_dict[accession] = (key, sequence)
    pprint.pprint(accession_dict)
    pickle.dump(accession_dict, open(accession_pickle_file_path, "wb"))


def get_lookup_row_and_column(accession_key1, accession_key2):
    lookup_row = accession_key1
    lookup_col = accession_key2
    if accession_key2 > accession_key1:
        lookup_row = accession_key2
        lookup_col = accession_key1
    return lookup_row, lookup_col


def compare_accessions(gene_list, accession_dict, matrix_k):
    in_cluster_minHash = SortedSet()
    out_cluster_minHash = SortedSet()
    for pair in itertools.combinations(gene_list, r=2):
        accession1 = accession_dict[pair[0]][0]
        accession2 = accession_dict[pair[1]][0]
        lookup_row, lookup_col = get_lookup_row_and_column(accession1, accession2)
        in_cluster_minHash.add(matrix_k[lookup_row][lookup_col])
        for key, value in accession_dict.items():
            if accession1 != key:
                lookup_row, lookup_col = get_lookup_row_and_column(accession1, key)
                out_cluster_minHash.add(matrix_k[lookup_row][lookup_col])
            if accession2 != key:
                lookup_row, lookup_col = get_lookup_row_and_column(accession2, key)
                out_cluster_minHash.add(matrix_k[lookup_row][lookup_col])
    return in_cluster_minHash[0], in_cluster_minHash[len(in_cluster_minHash) - 1], \
           out_cluster_minHash[0], out_cluster_minHash[len(out_cluster_minHash) - 1]


def calc_minHash_stats(pickle_file_path, accession_dict, matrix_k, cluster_file_path):
    sequence_dict = pickle.load(open(pickle_file_path, "rb"))
    # accession_dict = pickle.load(open(accession_pickle_file_path, "rb"))
    # matrix_k = np.loadtxt(matrix_output_path + "matrix_k6.txt")
    cluster_file_paths = get_cluster_filenames_from_directory(cluster_file_path)
    file_stats = list()
    for path in cluster_file_paths:
        split_path = str(path).split("/")
        cluster_name = re.findall("\d+", split_path[len(split_path) - 1])
        num_records = 0
        gene_list = list()
        for record in SeqIO.parse(path, "fasta"):
            gene_list.append(record.id)
            num_records = num_records + 1
        if num_records > 2:
            print("Checking cluster #" + str(cluster_name))
            in_max, in_min, out_max, out_min = compare_accessions(gene_list, accession_dict, matrix_k)
            file_stats.append((cluster_name, num_records, in_max, in_min, out_max, out_min))
    with open(results_output_path, newline='') as csvfile:
        results_writer = csv.writer(csvfile, delimiter=' ', quoting=csv.QUOTE_NONE)
        for result in file_stats:
            results_writer.writerow = result


def print_to_matrices(sequence_dict, matrix_output_path, from_k, to_k):
    # print_minhash_to_pickle(3, 11, sequence_dict)
    print(len(sequence_dict))
    numSequences = len(sequence_dict)
    matrix_length = numSequences + 1
    # matrix_dictionary = create_matrices(2, 11, matrix_length)

    if not os.path.exists(matrix_output_path):
        command = input("External drive %s doesn't exist..continue? (y/n)" % matrix_output_path)
        if command == 'n':
            exit("exiting early because the external drive doesn't exist")
    else:
        print("The external directory exists.")
    print("Running from k = " + str(from_k) + " to k = " + str(to_k))

    for k in range(from_k, to_k):
        jaccard_matrix = np.empty(shape=(matrix_length, matrix_length), dtype=np.float32)
        # pickle.dump(jaccard_matrix, open(matrix_output_path + "test_matrix.p", "wb"), protocol=3)
        # np.savetxt(matrix_output_path + "test_matrix.txt", jaccard_matrix)
        for i in range(1, matrix_length):
            jaccard_matrix[i][0] = i
            jaccard_matrix[0][i] = i
        pickle_file = get_minhash_pickle_filename(k, basePath)
        sequence_minhash = pickle.load(open(pickle_file, "rb"))
        np.fill_diagonal(jaccard_matrix, 1)
        for i in range(2, matrix_length - 1):
            print(i)
            for j in range(1, i + 1):
                minhash_col = sequence_minhash[jaccard_matrix[0][j]]
                minhash_row = sequence_minhash[jaccard_matrix[i][0]]
                # minhash_value = minhash_row.jaccard(minhash_col)
                # print("Minhash value: " + str(minhash_value))
                jaccard_matrix[i][j] = minhash_row.jaccard(minhash_col)
                # pprint.pprint(jaccard_matrix[i])
        np.savetxt(matrix_output_path + "matrix_k" + str(k) + ".txt", jaccard_matrix)
        del jaccard_matrix  # protect against possible MemoryError--doesn't do anything
        # pickle.dump(jaccard_matrix, open(matrix_output_path + "matrix_k" + str(k) + ".p", "wb"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-k_base', type=int, required=True, help='base k value for running any calculations.')
    parser.add_argument('-k_inc', help='max k value increment. example if k_from is 6, and inc is 2, '
                                       'you should expect 6,7')
    parser.add_argument('-env', choices=['nick', 'npr', 'pbs', 'lsb456'], required=True,
                        help='The environment/machine you\'re working on. Helpful for filepaths.')
    parser.add_argument('-m', choices=['matrix', 'stats'], help='Which programs you would like to run')

    # matrix_output_path = "/media/CP_MyBook/Pickle_Matrices/"
    matrix_output_path = "/media/catherine/ExtraDrive1/Network_Matrices/"
    # matrix_output_path = "/media/catherine/My Book/Network_Matrices/"
    pickle_file_path = basePath + "all_sequences.p"
    accession_pickle_file_path = basePath + "accession_dict.p"
    results_output_path = basePath + "test_matrix.csv"
    clusters_input_path = basePath + "USearch_AA_Centroids/clusters_35"

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
    elif environment == "npr":
        base_path = "/home/catherine/Networks_Nick_NPR/"
        matrix_output_path = "/media/catherine/ExtraDrive1/Network_Matrices/"
    elif environment == "pbs":
        base_path = "/home/catherine/Networks_Nick/"
        matrix_output_path = "/media/catherine/My Book/Network_Matrices/"
    elif environment == "lsb456":
        base_path = "/home/lsb456/Networks_nick/"
        matrix_output_path = "/media/CP_MyBook/Pickle_Matrices/"

    if not os.path.exists(base_path):
        print(f"The path \"{base_path}\" does not exist on this machine. Are you in the right environment?")

    pickle_file_path = basePath + "all_sequences.p"
    accession_pickle_file_path = basePath + "accession_dict.p"
    results_output_path = basePath + "test_matrix.csv"
    clusters_input_path = basePath + "USearch_AA_Centroids/clusters_35"  # This is only on NPR for right now

    if mode == "matrix":
        sequence_dict = pickle.load(open(pickle_file_path, "rb"))
        print_to_matrices(sequence_dict, matrix_output_path, from_k, to_k)
    elif mode == "stats":
        if not os.path.exists(clusters_input_path):
            print("cluster input path %s does not exist. exiting.." % clusters_input_path)
            exit(0)
        accession_dict = pickle.load(open(accession_pickle_file_path, "rb"))
        print("Reading matrix...")
        matrix_k = np.loadtxt(matrix_output_path + "matrix_k" + str(from_k) + ".txt")
        print("Matrix reading finished")
        calc_minHash_stats(pickle_file_path, accession_dict, matrix_k, clusters_input_path)
    else:
        print("No mode specified.")


main()
