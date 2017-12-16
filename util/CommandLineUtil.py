

def parse_environment(environment, matrix_k_value, thresholds, cluster_output_directory):
    all_cluster_paths = list()
    cluster_io_paths = list()
    if environment == "nick":
        base_path = "/Users/nickpredey/Documents/Networks/"
        matrix_output_path = "/Users/nickpredey/Documents/Networks/test_matrix.txt"
        results_output_path = base_path + "test_matrix_{0}.csv"
    elif environment == "npr":
        base_path = "/home/catherine/Networks_Nick_NPR/"
        results_output_path = cluster_output_directory + 'cluster_{0}_results.csv'
        matrix_input_path = "/data/matrix_k{0}.txt".format(matrix_k_value)
        clusters_input_path = base_path + "USearch_AA_Centroids/clusters_{}"
        # matrix_output_path = "/media/catherine/ExtraDrive1/Network_Matrices/" #nprito, I believe
    elif environment == "pbs":
        base_path = "/home/catherine/Networks_Nick/"
        matrix_output_path = "/media/catherine/My Book/Network_Matrices/"
        results_output_path = base_path + "test_matrix_40.csv"
    elif environment == "lsb456":
        base_path = "/home/lsb456/Networks_nick/"
        results_output_path = base_path + "test_matrix_40.csv"
        matrix_output_path = "/media/CP_MyBook/Pickle_Matrices/"

    for threshold in thresholds:
        current_input = clusters_input_path.format(threshold)
        current_output = results_output_path.format(threshold)
        cluster_io_paths.append((current_input, current_output))
    return base_path, cluster_io_paths, matrix_input_path