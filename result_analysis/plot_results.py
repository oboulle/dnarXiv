
import os
import sys
import matplotlib.pyplot as plt

def compute_array(array, repeat_nbr):
    computed_array = []
    
    for i in range(int(len(array)/repeat_nbr)):
        sub_array = array[repeat_nbr*i:repeat_nbr*(i+1)]
        while None in sub_array: sub_array.remove(None)
        if sub_array == []: computed_array.append(None)
        else: computed_array.append(sum(sub_array) / len(sub_array))
    return computed_array


def compute_results(results_path, repeat_nbr):
    
    if not os.path.isfile(results_path):
        print("error : file not found : "+results_path)
        exit(1)
    
    source_length_array = []
    n_mol_array = []
    precision_array = []
    reading_time_array = []
    
    results_file = open(results_path)
    line = results_file.readline()
    while line != "": #read the file
        source_length = results_file.readline().replace("source_length ", "").replace("\n","")
        n_mol = results_file.readline().replace("n_mol ", "").replace("\n","")
        precision = results_file.readline().replace("precision ", "").replace("\n","")
        reading_time = results_file.readline().replace("reading_time ", "").replace("\n","")
        
        if source_length != "None":
            source_length_array.append(int(source_length))
        else:
            source_length_array.append(None)
        if n_mol != "None":
            n_mol_array.append(int(n_mol))
        else:
            n_mol_array.append(None)
        if precision != "None":
            precision_array.append(float(precision))
        else: #precision = None -> no result has been found -> 0% precision
            precision_array.append(0)
        if reading_time != "None":
            reading_time_array.append(int(reading_time))
        else:
            reading_time_array.append(None)
        
        line = results_file.readline()
        line = results_file.readline()
    
    computed_source_length_array = compute_array(source_length_array, repeat_nbr)
    computed_n_mol_array = compute_array(n_mol_array,repeat_nbr)
    computed_precision_array = compute_array(precision_array, repeat_nbr)
    computed_reading_time_array = compute_array(reading_time_array, repeat_nbr)
    return computed_source_length_array, computed_n_mol_array, computed_precision_array, computed_reading_time_array


def plot_arrays(array_x, array_y, label_x, label_y):

    plt.plot(array_x, array_y)
    
    plt.xlabel(label_x)
    plt.ylabel(label_y)
    plt.show()

# =================== main ======================= #
if __name__ == '__main__':
    
    if len(sys.argv) != 2:
        print("usage : plot_results.py results_path")
        sys.exit(1)

    results_path = sys.argv[1]
    repeat_nbr = 5 #number of time the scripts have been repeated with the same parameters
    source_length, n_mol, precision, reading_time = compute_results(results_path, repeat_nbr)
    
    plot_arrays(n_mol, reading_time, "n_mol", "reading_time (s)")
    plot_arrays(n_mol, precision, "n_mol", "precision")

    
    print("\tcompleted !")
