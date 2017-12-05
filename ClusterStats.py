# import psycopg2
import csv

from Bio import SeqIO
import itertools
from sortedcontainers import SortedSet
import re
import argparse
from util import *


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


# gene_list is a list of
def compare_accessions(cluster_gene_list, accession_dict, matrix_k):
    in_cluster_minHash = SortedSet()
    out_cluster_minHash = SortedSet()
    cluster_gene_list_as_int = set()
    for gene in cluster_gene_list:
        cluster_gene_list_as_int.add(accession_dict[gene][0])
    for pair in itertools.combinations(cluster_gene_list, r=2):
        accession1 = int(accession_dict[pair[0]][0])
        # print("accession1: " + str(accession1))
        accession2 = int(accession_dict[pair[1]][0])
        # print("accession2: " + str(accession2))
        lookup_row, lookup_col = getLookupRowAndColumn(accession1, accession2)
        in_cluster_minHash.add(matrix_k[lookup_row][lookup_col])
        # this iteration calculates for outside of the cluster
        for key, value in accession_dict.items():
            outside_cluster_index = int(value[0])
            if outside_cluster_index in cluster_gene_list_as_int:
                continue
            if accession1 != outside_cluster_index:
                lookup_row, lookup_col = getLookupRowAndColumn(accession1, outside_cluster_index)
                out_cluster_minHash.add(matrix_k[lookup_row][lookup_col])
            if accession2 != outside_cluster_index:
                lookup_row, lookup_col = getLookupRowAndColumn(accession2, outside_cluster_index)
                out_cluster_minHash.add(matrix_k[lookup_row][lookup_col])
    return in_cluster_minHash[len(in_cluster_minHash) - 1], in_cluster_minHash[0], \
           out_cluster_minHash[len(out_cluster_minHash) - 1], out_cluster_minHash[0]


def calc_minHash_stats(pickle_file_path, accession_dict, matrix_k, cluster_file_path, results_output_path):
    sequence_dict = pickle.load(open(pickle_file_path, "rb"))
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
    with open(results_output_path, 'w', newline='') as csvfile:
        results_writer = csv.writer(csvfile, delimiter=' ', quoting=csv.QUOTE_NONE)
        for result in file_stats:
            results_writer.writerow(result)


# TODO let calculations run en masse, and print to stats file accordingly. make choosing the cluster directory extensible
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-env', choices=['nick', 'npr', 'pbs', 'lsb456'], required=True,
                        help='The environment/machine you\'re working on. Helpful for filepaths.')
    #  parser.add_argument('-m', choices=['stats'], help='Which programs you would like to run')
    parser.add_argument('-clust_id', help='The cluster id (threshold) to run.')
    parser.add_argument('-k', help='The k value for the matrix.')
    if not len(sys.argv) > 1:
        print("no arguments specified. Refer to -h or --help.")
        exit(0)
    args = parser.parse_args()
    environment = args.env
    #  mode = args.m

    base_path = ""
    matrix_output_path = ""
    kmer_matrix_value = 6
    cluster_id_value = 40
    if args.k:
        kmer_matrix_value = args.k
    else:
        print('No k value specified. Defaulting to {}'.format(kmer_matrix_value))
    if args.clust_id:
        cluster_id_value = args.clust_id
    else:
        print('No cluster id value specified. Defaulting to {}'.format(cluster_id_value))
    if environment == "nick":
        base_path = "/Users/nickpredey/Documents/Networks/"
        matrix_output_path = "NOT YET IMPLEMENTED"
        results_output_path = base_path + "test_matrix_40.csv"
    elif environment == "npr":
        base_path = "/home/catherine/Networks_Nick_NPR/"
        results_output_path = '{0}test_matrix_{1}.csv'.format(base_path, cluster_id_value)
        matrix_output_path = "/data/matrix_k{0}.txt".format(kmer_matrix_value)
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
        print("The path {}  does not exist on this machine. Are you in the right environment?".format(base_path))

    pickle_file_path = base_path + "all_sequences.p"
    accession_pickle_file_path = base_path + "accession_dict.p"
    # results_output_path = basePath + "test_matrix.csv"
    clusters_input_path = "{0}USearch_AA_Centroids/clusters_{1}".format(base_path, cluster_id_value)

    if not os.path.exists(clusters_input_path):
        print("cluster input path {} does not exist. exiting..".format(clusters_input_path))
        exit(0)
    accession_dict = pickle.load(open(accession_pickle_file_path, "rb"))
    print("Reading matrix...")
    matrix_k = np.loadtxt(matrix_output_path)
    print("Matrix reading finished")

    print(
        "Running with parameters: ID Threshold= {0}\nMatrix k value = {1}".format(cluster_id_value, kmer_matrix_value)
    )
    calc_minHash_stats(pickle_file_path, accession_dict, matrix_k, clusters_input_path, results_output_path)


main()
