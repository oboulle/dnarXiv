
import os
import sys
import matplotlib.pyplot as plt


def compute_results(results_path):
    
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
        line = results_file.readline()
        line = results_file.readline()
            
        if source_length != "None" and n_mol != "None":
            source_length_array.append(int(source_length))
            n_mol_array.append(int(n_mol))
        else:
            continue #skip the line
        
        if precision != "None":
            precision_array.append(float(precision))
        else: #precision = None -> no result has been found -> 0% precision
            precision_array.append(0)
        
        if reading_time != "None":
            reading_time_array.append(int(reading_time))
        else:
            reading_time_array.append(0)

    return source_length_array, n_mol_array, precision_array, reading_time_array


def plot_arrays(array_x, array_y, label_x, label_y):

    plt.plot(array_x, array_y)
    
    plt.xlabel(label_x)
    plt.ylabel(label_y)

    if label_y == "precision": plt.ylim(0,100)
    plt.show()

# =================== main ======================= #
if __name__ == '__main__':
    
    if len(sys.argv) != 2:
        print("usage : plot_results.py results_path")
        sys.exit(1)
    
    results_path = sys.argv[1]
    source_length, n_mol, precision, reading_time = compute_results(results_path)
    
    plot_arrays(n_mol, reading_time, "n_mol", "reading_time (s)")
    plot_arrays(n_mol, precision, "n_mol", "precision")
    
    print("\tcompleted !")
