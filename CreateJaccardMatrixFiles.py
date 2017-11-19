#import psycopg2
import csv
from datasketch import MinHash
import numpy as np
import _pickle as pickle
import pprint
import os
from Bio import SeqIO
import itertools
from sortedcontainers import SortedSet
import re

basePath = "/Users/nickpredey/Documents/Networks/PickleFilesMinHash/"
#basePath = "/home/lsb456/Networks_nick/"
#basePath = "/home/catherine/Networks_Nick/"
#matrix_output_path = "/media/CP_MyBook/Pickle_Matrices/"
matrix_output_path = "/media/catherine/My Book/Network_Matrices/"
pickle_file_path = basePath + "all_sequences.p"
accession_pickle_file_path = basePath + "accession_dict.p"
results_output_path = basePath + "test_matrix.csv"
from_k = 7
to_k = 9


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


def getJaccardIndex(sequence1, sequence2, k):
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


def getActualJaccardIndex(kmers1, kmers2):
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
    cursor.execute("SELECT sequence_id, accession_number, seq from sequences")
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
def getClusterFilenamesFromDirectory(clusterFileDirectory):
    file_paths = []
    for folder, subs, files in os.walk(clusterFileDirectory):
        for filename in files:
            file_paths.append(os.path.abspath(os.path.join(folder, filename)))
    return file_paths

def getAndPrintAccessionsToPickle(sequence_dictionary):
    accession_dict = dict()
    for key, value in sequence_dictionary.items():
        accession = value[0]
        sequence = value[1]
        accession_dict[accession] = (key, sequence)
    pprint.pprint(accession_dict)
    pickle.dump(accession_dict, open(accession_pickle_file_path, "wb"))

def getLookupRowAndColumn(accession_key1, accession_key2):
    lookup_row = accession_key1
    lookup_col = accession_key2
    if accession_key2 > accession_key1:
        lookup_row = accession_key2
        lookup_col = accession_key1
    return lookup_row, lookup_col


def compareAccessions(gene_list, accession_dict, matrix_k):
    in_cluster_minHash = SortedSet()
    out_cluster_minHash = SortedSet()
    for pair in itertools.combinations(gene_list, r=2):
        accession1 = accession_dict[pair[0]][0]
        accession2 = accession_dict[pair[1]][0]
        lookup_row, lookup_col = getLookupRowAndColumn(accession1, accession2)
        in_cluster_minHash.add(matrix_k[lookup_row][lookup_col])
        for key, value in accession_dict.items():
            if accession1 != key:
                lookup_row, lookup_col = getLookupRowAndColumn(accession1, key)
                out_cluster_minHash.add(matrix_k[lookup_row][lookup_col])
            if accession2 != key:
                lookup_row, lookup_col = getLookupRowAndColumn(accession2, key)
                out_cluster_minHash.add(matrix_k[lookup_row][lookup_col])
    return in_cluster_minHash[0], in_cluster_minHash[len(in_cluster_minHash)-1], \
           out_cluster_minHash[0], out_cluster_minHash[len(out_cluster_minHash)-1]
sequence_dict = pickle.load(open(pickle_file_path, "rb"))
accession_dict = pickle.load(open(accession_pickle_file_path, "rb"))
matrix_k = np.loadtxt(matrix_output_path + "matrix_k6.txt")
cluster_file_paths = getClusterFilenamesFromDirectory("/Users/nickpredey/Documents/USearch_AA_PT35")
file_stats = list()
for path in cluster_file_paths:
    split_path = str(path).split("/")
    cluster_name = re.findall("\d+", split_path[len(split_path)-1])
    num_records = 0
    gene_list = list()
    for record in SeqIO.parse(path, "fasta"):
        gene_list.append(record.id)
        num_records = num_records + 1
    if num_records > 2:
        print("Checking cluster #" + str(cluster_name))
        in_max, in_min, out_max, out_min = compareAccessions(gene_list, accession_dict, matrix_k)
        file_stats.append((cluster_name, num_records, in_max, in_min, out_max, out_min))
with open(results_output_path, newline='') as csvfile:
    results_writer = csv.writer(csvfile, delimiter=' ', quoting=csv.QUOTE_NONE)
    for result in file_stats:
        results_writer.writerow = result
exit()
#print_minhash_to_pickle(3, 11, sequence_dict)
print(len(sequence_dict))
numSequences = len(sequence_dict)
matrix_length = numSequences + 1
#matrix_dictionary = create_matrices(2, 11, matrix_length)

if not os.path.exists(matrix_output_path):
    command = input("External drive %s doesn't exist..continue? (y/n)" % matrix_output_path)
    if command == 'n':
        exit("exiting early because the external drive doesn't exist")
else:
    print("The external directory exists.")
print("Running from k = " + str(from_k) + " to k = " + str(to_k))

for k in range(from_k, to_k):
    jaccard_matrix = np.empty(shape=(matrix_length, matrix_length), dtype=np.float32)
    #pickle.dump(jaccard_matrix, open(matrix_output_path + "test_matrix.p", "wb"), protocol=3)
    #np.savetxt(matrix_output_path + "test_matrix.txt", jaccard_matrix)
    for i in range(1,  matrix_length):
        jaccard_matrix[i][0] = i
        jaccard_matrix[0][i] = i
    pickle_file = get_minhash_pickle_filename(k, basePath)
    sequence_minhash = pickle.load(open(pickle_file, "rb"))
    np.fill_diagonal(jaccard_matrix, 1)
    for i in range(2, matrix_length-1):
        print(i)
        for j in range(1, i+1):
            minhash_col = sequence_minhash[jaccard_matrix[0][j]]
            minhash_row = sequence_minhash[jaccard_matrix[i][0]]
           # minhash_value = minhash_row.jaccard(minhash_col)
           # print("Minhash value: " + str(minhash_value))
            jaccard_matrix[i][j] = minhash_row.jaccard(minhash_col)
       #pprint.pprint(jaccard_matrix[i])
    np.savetxt(matrix_output_path + "matrix_k" + str(k) + ".txt", jaccard_matrix)
    del jaccard_matrix #protect against possible MemoryError--doesn't do anything
    #pickle.dump(jaccard_matrix, open(matrix_output_path + "matrix_k" + str(k) + ".p", "wb"))
exit()

accession_number_to_int_keys = list(accession_number_dict.keys())
for i in range(1, matrix_length):
    jaccard_matrix[i][0] = accession_number_to_int_keys[sequence_index]
    jaccard_matrix[0][i] = accession_number_to_int_keys[sequence_index]
    sequence_index = sequence_index + 1
for i in range(1, matrix_length):
    print(i)
    sequence_1 = accession_number_dict[jaccard_matrix[i][0]][1]
    for j in range(1, matrix_length):
        print(j)
        sequence_2 = accession_number_dict[jaccard_matrix[0][i]][1]
        jaccard_matrix[i][j] = getJaccardIndex(sequence_1, sequence_2, k)

rowtowrite = list(jaccard_matrix[0])
print(rowtowrite)
exit()
with open('/Users/nickpredey/Documents/Networks/jaccard.csv', newline='') as csvfile:
    jaccard_writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE)
    #row

