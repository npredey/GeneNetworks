import _pickle as pickle
import os
import numpy as np


def get_and_print_accessions_to_pickle(sequence_dictionary):
    accession_dict = dict()
    for key, value in sequence_dictionary.items():
        accession = value[0]
        sequence = value[1]
        accession_dict[accession] = (key, sequence)
    pprint.pprint(accession_dict)
    pickle.dump(accession_dict, open(accession_pickle_file_path, "wb"))


def output_path_exist(matrix_output_path, from_k, to_k):
    if not os.path.exists(matrix_output_path):
        command = input("External drive {} doesn't exist..continue? (y/n)".format(matrix_output_path))
        if command == 'n':
            exit("exiting early because the external drive doesn't exist")
    else:
        print("The external directory exists.")
    print("Running from k = " + str(from_k) + " to k = " + str(to_k))


def print_to_matrices(sequence_dict, matrix_output_path, from_k, to_k, base_path):
    # print_minhash_to_pickle(3, 11, sequence_dict)
    print(len(sequence_dict))
    num_sequences = len(sequence_dict)
    matrix_length = num_sequences + 1
    # matrix_dictionary = create_matrices(2, 11, matrix_length)
    output_path_exist(matrix_output_path, from_k, to_k)

    for k in range(from_k, to_k):
        jaccard_matrix = np.empty(shape=(matrix_length, matrix_length), dtype=np.float32)
        # pickle.dump(jaccard_matrix, open(matrix_output_path + "test_matrix.p", "wb"), protocol=3)
        # np.savetxt(matrix_output_path + "test_matrix.txt", jaccard_matrix)
        for i in range(1, matrix_length):
            jaccard_matrix[i][0] = i
            jaccard_matrix[0][i] = i

        pickle_file = get_minhash_pickle_filename(k, base_path)
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
