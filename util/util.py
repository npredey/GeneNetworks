import _pickle as pickle
import os

import numpy as np
from datasketch import MinHash


def build_filepath(base_path, file):
    return base_path + '/' + file


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
def get_filenames_from_directory(file_directory):
    file_paths = []
    for folder, subs, files in os.walk(file_directory):
        for filename in files:
            if 'DS_Store' not in filename:
                file_paths.append(os.path.abspath(os.path.join(folder, filename)))
    return file_paths
