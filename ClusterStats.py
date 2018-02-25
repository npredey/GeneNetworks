# import psycopg2
import argparse
import itertools
import re
import sys
import csv

from Bio import SeqIO
from sortedcontainers import SortedSet

from util.CommandLineUtil import *
from util.util import *


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
        if os.path.exists(cluster_input_dir):
            cluster_file_paths = get_filenames_from_directory(cluster_input_dir)
            output_string = ''
            for path in cluster_file_paths:
                split_path = str(path).split("/")
                cluster_name = re.findall("\d+", split_path[len(split_path) - 1])[0]
                num_records = 0
                gene_list = list()
                for record in SeqIO.parse(path, "fasta"):
                    gene_list.append(record.id)
                    num_records = num_records + 1
                if num_records > 2:
                    print("Checking cluster #{}".format(cluster_name))
                    in_max, in_min, out_max, out_min = compare_accessions(gene_list, accession_dict, matrix_k)
                    comparison_result = "{} {} {} {} {} {}".format(cluster_name, num_records, in_max,
                                                                   in_min, out_max, out_min)
                    if print_header:
                        output_string += '{}\n'.format(header)
                        print_header = False
                    print(comparison_result)
                    output_string += comparison_result + '\n'
            with open(cluster_result_output_path, 'a', newline='') as csvfile:
                csvfile.write(output_string + '\n')
        else:
            print("cluster input path {} does not exist. terminating this run.".format(cluster_input_dir))


def check_cluster_sequence_length(cluster_directories, clust_name, id_threshold):
    cluster_dir_to_check = ''
    for directory in cluster_directories:
        if str(id_threshold) in directory:
            cluster_dir_to_check = directory
            break
    cluster_file_to_check = cluster_dir_to_check + '/clusters' + clust_name

    sequences_in_cluster = list()
    for record in SeqIO.parse(cluster_file_to_check, 'fasta'):
        sequences_in_cluster.append(len(record.seq))

    check_result = ''
    length_difference = max(sequences_in_cluster) - min(sequences_in_cluster)
    if length_difference < min(sequences_in_cluster) / 2:
        check_result = "Length_difference=" + str(length_difference)
    return check_result


# parse a directory containing cluster files from USEARCH.
# returns a tuple of
def parse_clusters(cluster_files_directory, output_file):
    output_header = "threshold,num_clusters,average_cluster_size,min_size,max_size,num_>2"
    print(cluster_files_directory)
    print_header = True
    output_string = ''
    for base_dir, subs, files in os.walk(cluster_files_directory):
        for sub in subs:
            threshold = sub.split('_')[-1]
            num_clusters = 0
            average_cluster_size = 0
            min_size = 100000
            max_size = 0
            num_greater_2 = 0
            num_genes = 0

            for cluster_folder, subs1, _files in os.walk(build_filepath(base_dir, sub)):
                print(threshold)
                for file in _files:
                    cluster_file = build_filepath(base_dir + '/' + sub, file)
                    records = list(SeqIO.parse(cluster_file, 'fasta'))
                    records_length = len(records)
                    num_genes += records_length
                    num_clusters += 1
                    if records_length > max_size:
                        max_size = records_length
                    elif records_length < min_size:
                        min_size = records_length
                    if records_length > 2:
                        num_greater_2 += 1

            # TODO median cluster size might be more insightful than the mean cluster size
            average_cluster_size = int(float(num_genes) / num_clusters)
            output_string += "{},{},{},{},{},{}\n".format(threshold, num_clusters, average_cluster_size, min_size,
                                                          max_size, num_greater_2)
    with open(output_file, 'w') as output:
        output.write(output_header + '\n')
        output.write(output_string)


def get_cluster_stats(cluster_file_list):
    cluster_name_index = 0
    cluster_size_index = 1
    max_in_index = 2
    min_in_index = 3
    max_out_index = 4
    min_out_index = 5
    for file in cluster_file_list:
        output_file_path = str(file).split(".")[0] + "_bad_cluster_list.csv"
        id_threshold = int(''.join(filter(str.isdigit,
                                          file.split('/')[-1])))

        k_value = int(''.join(filter(str.isdigit,
                                     file.split('/')[-2])))

        if 'bad' in str(file):
            continue
        with open(file, "r") as infile:
            csv_reader = csv.reader(infile, delimiter=' ')
            is_header = True
            num_bad_clusters = 0
            num_good_clusters = 0
            bad_but_correctable_cluster = 0
            for row in csv_reader:
                # print(row)
                if is_header:
                    # if 'cluster_id' not in csv_reader.read():
                    #     print('File {} is not in the correct format.\n {} is the problem line'.format(file, line))
                    #     break
                    # else:
                    is_header = False
                    outfile = open(output_file_path, 'w+')
                    outfile.write('max_in min_in difference\n')
                else:
                    # split_line = line.split(" ")
                    if len(row) < 5 or not any(char.isdigit() for char in row[0]):
                        continue
                    clust_name = row[cluster_name_index]
                    bad_cluster = False
                    output_line = ''
                    if float(row[max_out_index]) > float(row[max_in_index]):
                        output_line += clust_name + ' '
                        bad_cluster = True
                        # outfile.write("{}\n".format(clust_name.decode('utf-8')))

                    if float(row[max_out_index]) > float(row[min_in_index]):
                        output_line += clust_name + ' '
                        bad_cluster = True

                    if bad_cluster:
                        num_bad_clusters += 1
                        check_result = check_cluster_sequence_length(directory_names, clust_name, id_threshold)
                        if check_result == '':
                            bad_but_correctable_cluster += 1
                            continue
                        outfile.write("{} {}\n".format(output_line, check_result))
                        bad_cluster = False
                    else:
                        num_good_clusters += 1


def parse_cluster_stats(cluster_stats_file_path, clusters_dir, graph_output_file):
    graph_output_header = "k, threshold, percent_bad, correctable"
    graph_output_string = ''
    cluster_name_index = 0
    cluster_size_index = 1
    max_in_index = 2
    min_in_index = 3
    max_out_index = 4
    min_out_index = 5

    num_bigger_max_in = 0
    num_bigger_min_in = 0

    cluster_file_list = list()
    if os.path.isdir(cluster_stats_file_path):
        cluster_file_list = get_filenames_from_directory(cluster_stats_file_path)
    else:
        cluster_file_list.append(cluster_stats_file_path)

    if not os.path.exists(clusters_dir):
        print("Cluster directory {} does not exist. Exiting...".format(clusters_dir))
        return
    directory_names = [val[0] for val in os.walk(clusters_dir)]
    print_graph_header = True
    for file in cluster_file_list:
        output_file_path = str(file).split(".")[0] + "_bad_cluster_list.csv"
        id_threshold = int(''.join(filter(str.isdigit,
                                          file.split('/')[-1])))

        k_value = int(''.join(filter(str.isdigit,
                                     file.split('/')[-2])))

        if 'bad' in str(file):
            continue
        with open(file, "r") as infile:
            csv_reader = csv.reader(infile, delimiter=' ')
            is_header = True
            num_bad_clusters = 0
            num_good_clusters = 0
            bad_but_correctable_cluster = 0
            for row in csv_reader:
                # print(row)
                if is_header:
                    # if 'cluster_id' not in csv_reader.read():
                    #     print('File {} is not in the correct format.\n {} is the problem line'.format(file, line))
                    #     break
                    # else:
                    is_header = False
                    outfile = open(output_file_path, 'w+')
                    outfile.write('max_in min_in difference\n')
                else:
                    # split_line = line.split(" ")
                    if len(row) < 5 or not any(char.isdigit() for char in row[0]):
                        continue
                    clust_name = row[cluster_name_index]
                    bad_cluster = False
                    output_line = ''
                    if float(row[max_out_index]) > float(row[max_in_index]):
                        output_line += clust_name + ' '
                        bad_cluster = True
                        # outfile.write("{}\n".format(clust_name.decode('utf-8')))

                    if float(row[max_out_index]) > float(row[min_in_index]):
                        output_line += clust_name + ' '
                        bad_cluster = True

                    if bad_cluster:
                        num_bad_clusters += 1
                        check_result = check_cluster_sequence_length(directory_names, clust_name, id_threshold)
                        if check_result == '':
                            bad_but_correctable_cluster += 1
                            continue
                        outfile.write("{} {}\n".format(output_line, check_result))
                        bad_cluster = False
                    else:
                        num_good_clusters += 1

            percent_bad = -1
            if not is_header:
                total_clusters = num_good_clusters + num_bad_clusters
                percent_bad = float(num_bad_clusters) / total_clusters
                outfile.write("Num bad clusters: {}\nNum Good clusters: {}\ntotal: {}\n Percent bad{} ".format(
                    num_bad_clusters, num_good_clusters, total_clusters, percent_bad))
                outfile.close()

                bad_but_correctable_cluster = float(bad_but_correctable_cluster) / total_clusters
                graph_output_string += "{},{},{},{}\n".format(k_value, id_threshold, percent_bad,
                                                              bad_but_correctable_cluster)

            # with open(graph_output_file, 'w') as graph_output:
            #     if print_graph_header:
            #         graph_output.write(graph_output_header + '\n')
            #         print_graph_header = False
            #     bad_but_correctable_cluster = float(bad_but_correctable_cluster)/total_clusters
            #     graph_output.write("{},{},{},{}\n".format(k_value, id_threshold, percent_bad, bad_but_correctable_cluster))
    with open(graph_output_file, 'a') as graph_output:
        graph_output.write(graph_output_header + '\n')
        graph_output.write(graph_output_string)


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
    parser.add_argument('-diroutput',
                        help='Directory name to hold each cluster output by matrix k value.'
                             'it will be created if it does not already exist. Should contain'
                             'information about the k value of the matrix. For example, Matrix_k7_stats.'
                             ' In general, output_dir will be the output directory of any program running.',
                        dest='output_directory', type=str)
    parser.add_argument('-foutput', help='This parameter describes the output file for any one of the programs.',
                        dest='output_file')
    parser.add_argument('-t', help='Flag to force an empty matrix for testing purposes.', default=False)
    parser.add_argument('-p', help='Parse a directory or file of cluster stats', type=str)
    parser.add_argument('-pc', help='Parse a directory of clusters (from USEARCH or similar program)', type=str)
    parser.add_argument('-clusters_dir', help='Cluster files directory for parsing the cluster stats files.', type=str)
    if not len(sys.argv) > 1:
        print("no arguments specified. Refer to -h or --help.")
        exit(0)
    args = parser.parse_args()
    environment = args.env
    output_directory = args.output_directory
    output_file = args.output_file
    clusters_dir = args.clusters_dir
    if args.p:
        parse_cluster_stats(args.p, clusters_dir, output_file)
    elif args.pc:
        parse_clusters(args.pc, output_file)

    elif args.k and args.thresholds:
        kmer_matrix_value = args.k
        thresholds = args.thresholds
        base_path, cluster_io_paths, matrix_output_path = parse_environment(environment,
                                                                            kmer_matrix_value,
                                                                            thresholds,
                                                                            output_directory)

        if args.output_directory and not os.path.exists(output_directory):
            print("Output directory {0} does not exist. It will be created.".format(output_directory))
            os.makedirs(output_directory)
        else:
            print("Must specify an output directory for running matrix stats. Refer to -h.")

        pickle_file_path = base_path + "all_sequences.p"
        accession_pickle_file_path = base_path + "accession_dict.p"

        accession_dict = pickle.load(open(accession_pickle_file_path, "rb"))

        if args.t:
            print("Reading matrix...")
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
