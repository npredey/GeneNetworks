# import psycopg2
import argparse
import itertools
import re
import sys

from Bio import SeqIO
from sortedcontainers import SortedSet

from util.CommandLineUtil import *
from util.util import *


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
        lookup_row, lookup_col = get_lookup_row_and_column(accession1, accession2)
        in_cluster_minHash.add(matrix_k[lookup_row][lookup_col])
        # this iteration calculates for outside of the cluster
        for key, value in accession_dict.items():
            outside_cluster_index = int(value[0])
            if outside_cluster_index not in cluster_gene_list_as_int:
                lookup_row, lookup_col = get_lookup_row_and_column(accession1, outside_cluster_index)
                out_cluster_minHash.add(matrix_k[lookup_row][lookup_col])

                lookup_row, lookup_col = get_lookup_row_and_column(accession2, outside_cluster_index)
                out_cluster_minHash.add(matrix_k[lookup_row][lookup_col])
    return in_cluster_minHash[len(in_cluster_minHash) - 1], in_cluster_minHash[0], \
           out_cluster_minHash[len(out_cluster_minHash) - 1], out_cluster_minHash[0]


def calc_minHash_stats(accession_dict, matrix_k, cluster_io_paths):
    header = 'cluster_id num_genes_in_cluster max_in min_in max_out min_out'
    for io_paths in cluster_io_paths:
        cluster_input_dir = io_paths[0]
        cluster_result_output_path = io_paths[1]
        print_header = True
        if not os.path.exists(cluster_input_dir):
            print("cluster input path {} does not exist. exiting..".format(cluster_input_dir))
            continue
        cluster_file_paths = get_cluster_filenames_from_directory(cluster_input_dir)
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
                with open(cluster_result_output_path, 'a', newline='') as csvfile:
                    output = "{} {} {} {} {} {}".format(str(cluster_name[0]), str(num_records), str(in_max),
                                                        str(in_min),
                                                        str(out_max), str(out_min))
                    if print_header:
                        csvfile.write('{0}\n'.format(header))
                        print_header = False
                    print(output)
                    csvfile.write(output + '\n')


def parse_cluster_stats(cluster_stats_file_path):
    cluster_name = 0
    cluster_size = 1
    max_in = 2
    min_in = 3
    max_out = 4
    min_out = 5

    cluster_file_list = list()
    if os.path.isdir(cluster_stats_file_path):
        cluster_file_list = get_cluster_filenames_from_directory(cluster_stats_file_path)
    else:
        cluster_file_list.append(cluster_stats_file_path)
    for file in cluster_file_list:
        output_file_path = str(file).split(".")[0] + "_bad_cluster_list.csv"
        bad_clusters = list()
        with open(file, "rb") as infile:
            is_header = True
            for line in infile:
                line = str(line)
                if is_header:
                    if 'cluster_id' not in line:
                        print("File is not in the correct format.")
                        break
                    else:
                        is_header = False
                        outfile = open(output_file_path, 'w+')
                else:
                    split_line = line.split(" ")
                    clust_name = split_line[cluster_name]
                    if float(split_line[max_out]) > float(split_line[max_in]):
                        outfile.write("{}\n".format(clust_name))
            if not is_header:
                outfile.close()


# TODO let calculations run en masse, and print to stats file accordingly. make choosing the cluster directory extensible
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-env', choices=['nick', 'npr', 'pbs', 'lsb456'], required=True,
                        help='The environment/machine you\'re working on. Helpful for filepaths.')
    #  parser.add_argument('-m', choices=['stats'], help='Which programs you would like to run')
    parser.add_argument('-k',
                        help='The k value for the matrix.')
    parser.add_argument(
        '-cluster_vals',
        nargs="*",  # expects ≥ 0 arguments
        type=int,
        default=[35, 40, 50, 60, 70, 80, 90],
        dest='thresholds'  # default to all cluster directories
    )
    parser.add_argument('-output_dir',
                        help='Directory name to hold each cluster output by matrix k value.'
                             'it will be created if it does not already exist. Should contain'
                             'information about the k value of the matrix',
                        dest='output_directory')
    parser.add_argument('-t', help='Flag to force an empty matrix for testing purposes.', default=False)
    parser.add_argument('-p', help='Parse a directory or file of cluster stats', type=str)
    if not len(sys.argv) > 1:
        print("no arguments specified. Refer to -h or --help.")
        exit(0)
    args = parser.parse_args()
    environment = args.env
    output_directory = args.output_directory
    #  mode = args.m
    if args.p:
        parse_cluster_stats(args.p)
        exit(0)
    thresholds = args.thresholds
    kmer_matrix_value = 6
    if args.k:
        kmer_matrix_value = args.k
    else:
        print('No k value specified. Defaulting to {}'.format(kmer_matrix_value))

    base_path, cluster_io_paths, matrix_output_path = parse_environment(environment,
                                                                        kmer_matrix_value,
                                                                        thresholds,
                                                                        output_directory)

    if not os.path.exists(base_path):
        print("The path {}  does not exist on this machine. Are you in the right environment?".format(base_path))

    if not os.path.exists(output_directory):
        print("Output directory {0} does not exist. It will be created.".format(output_directory))
        os.makedirs(output_directory)

    pickle_file_path = base_path + "all_sequences.p"
    accession_pickle_file_path = base_path + "accession_dict.p"

    accession_dict = pickle.load(open(accession_pickle_file_path, "rb"))

    print("Reading matrix...")
    if args.t:
        matrix_k = np.zeros(shape=(10 ** 5, 10 ** 5))
    else:
        matrix_k = np.loadtxt(matrix_output_path)
    print("Matrix reading finished")

    print(
        'Running with parameters: ID Thresholds= {0}'
        '\nMatrix k value = {1}'
        '\nOutput directory = {2}'.format(thresholds, kmer_matrix_value, output_directory)
    )
    calc_minHash_stats(accession_dict, matrix_k, cluster_io_paths)


main()
