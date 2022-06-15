import sys
import numpy as np
import random
from PIL import Image


"""
convert images to dna sequences
"""


base_4_to_dna_dict = {0: "A", 1: "C", 2: "G", 3: "T"}

base_16_to_dna_dict = {0: "AA", 1: "AC", 2: "AG", 3: "AT",
                       4: "CA", 5: "CC", 6: "CG", 7: "CT", 
                       8: "GA", 9: "GC", 10: "GG", 11: "GT", 
                       12: "TA", 13: "TC", 14: "TG", 15: "TT"}


#___Encoding

def encode_png(doc_path: str, dna_fragment_length: int, mode: str) -> (str, dict):
    """
    encode a png file to a sequence in the given mode
    'L' (grey shades), 'RGB', ...
    """ 
    mat = png_to_matrix(doc_path, mode)
    sequence, encoding_data = mat_to_sequence(mat, mode)
        
    # adjust size of dna sequence to get a round number of fragments
    rest = len(sequence) % (dna_fragment_length)
    if rest > 0: sequence += "A" * (dna_fragment_length - rest)
    
    pixel_nbr = int(encoding_data["height"])*int(encoding_data["width"])
    
    s_sequence = shuffle_sequence(sequence, pixel_nbr, mode)

    return s_sequence, encoding_data


def png_to_matrix(image_path: str, mode: str) -> np.array:
    """
    convert a png file to a matrix in the given mode
    'L' (grey shades), 'RGB', ...
    """ 
    img = Image.open(image_path).convert(mode)
    mat = np.array(img)
    return mat


def shuffle_sequence(sequence: str, pixel_nbr: int, mode: str) -> str:
    """
    shuffle the sequence representing the pixels of the matrix, so when fragments are lost, some pixels are randomly lost in all the matrix instead of a continuous line"
    """
    index_list = [k for k in range(pixel_nbr)]
    s_index_list = random.Random(0).sample(index_list, pixel_nbr) #the seed is set for the shuffle, it will be used for the unshuffle at the decoding
    
    if mode == "L":
        pixel_size_in_dna = 2
    else: #RGB
        pixel_size_in_dna = 3
        
    s_sequence = ['A' for k in range(len(sequence))]
    for i in range(pixel_nbr):
        new_index = s_index_list[i]
        s_sequence[pixel_size_in_dna*i:pixel_size_in_dna*(i+1)] = sequence[pixel_size_in_dna*new_index:pixel_size_in_dna*(new_index+1)]
    
    return "".join(s_sequence)


def unshuffle_sequence(sequence: str, pixel_nbr: int, pixel_size_in_dna: int) -> str:
    """
    remove the shuffle by replacing each dna_pixel at its original position
    create an index list and shuffle it with the same seed of the original shuffling
    use the new positions of the shuffled index to find the original positions of the dna pixels in the sequence
    """
    index_list = [k for k in range(pixel_nbr)] #init sorted list of index
    s_index_list = random.Random(0).sample(index_list, pixel_nbr) #shuffle the indexes with the same seed of the shuffling of the sequence

    us_sequence = ['_' for k in range(pixel_nbr*pixel_size_in_dna)] #init empty result sequence
    for i in range(pixel_nbr):
        old_index = s_index_list[i]
        us_sequence[pixel_size_in_dna*old_index:pixel_size_in_dna*(old_index+1)] = sequence[pixel_size_in_dna*i:pixel_size_in_dna*(i+1)] # s_index_list[4] = 56 -> the 4th element of sequence is 56th in the unshuffled sequence
    
    return "".join(us_sequence)


def pixel_to_dna_sequence(pxl_value, mode: str) -> str:
    """
    if mode = L, pxl value is an int 0-255
    if mode = RGB, pxl value is a list [0-255, 0-255, 0-255]
    """
    if mode == "L":
        rounded_pxl_value = int(pxl_value/16) # round the shade from 0-255 to 0-15
        kmer_value = base_16_to_dna_dict[rounded_pxl_value] #pair of nucleotide representing the value
    else: #RGB
        rounded_pxl_value_list = [ int(k/64) for k in pxl_value ] # round the color from 0-255 to 0-3
        kmer_value_list = [ base_4_to_dna_dict[k] for k in rounded_pxl_value_list ] #list of nucleotide representing the value
        kmer_value = "".join(kmer_value_list)
    return kmer_value


def mat_to_sequence(mat: np.array, mode: str) -> (str, dict):
    """
    convert a matrix of grey shades into a dna sequence
    """
    mat_height = len(mat)
    mat_width = len(mat[0])
    
    encoding_data = {"color" : mode, "height" : str(mat_height), "width" : str(mat_width)}
    kmer_list = []
    #create a list of all pixels converted to 2-kmers 
    for i in range(mat_height):
        for j in range(mat_width):
            pxl_value = mat[i][j] #pixel value representing a shade of grey (0-255)
            kmer_value = pixel_to_dna_sequence(pxl_value, mode) #pair of nucleotide representing the value
            kmer_list.append(kmer_value)
    
    sequence = "".join(kmer_list)
    
    return sequence, encoding_data


#___Decoding

def decode_png(sequence: str, output_path: str, encoding_data: dict) -> None:
    """
    convert a sequence to a png image
    """    
    img_mode = encoding_data["color"]
    img_height = int(encoding_data["height"])
    img_width = int(encoding_data["width"])
    
    if img_mode == "L":
        us_sequence = unshuffle_sequence(sequence, img_height*img_width, 2)
        mat = sequence_to_grey_mat(us_sequence, img_height, img_width)
    else: #RGB
        us_sequence = unshuffle_sequence(sequence, img_height*img_width, 3)
        mat = sequence_to_rgb_mat(us_sequence, img_height, img_width)
    matrix_to_png(mat, img_mode, output_path)


def sequence_to_grey_mat(sequence: str, mat_height: int, mat_width: int) -> list:
    """
    convert a sequence to a matrix representing a grey shaded image
    the given sequence is supposed to have the same height and the width *2 of the original image
    """
    mat = [ [] for k in range(mat_height) ]
    sequence_rows = [sequence[i:i+mat_width*2] for i in range(0, len(sequence), mat_width*2)]

    for i in range(mat_height):
        seq_row = sequence_rows[i]
        for j in range(mat_width):
            kmer = seq_row[2*j:2*(j+1)]
            pxl_value = -1 #count as a missing pixel if the kmer is not recognized
            for key, value in base_16_to_dna_dict.items():
                if value == kmer:
                    rounded_pxl_value = key
                    pxl_value = 17*rounded_pxl_value #from 0-15 to 0-255
                    break
            mat[i].append(pxl_value)

    def check_adjacent_pixel(i, j, sum, nbr_adj_pixel):
        # if pixel in matrix and not missing
        if i >= 0 and i < mat_height and j >= 0 and j < mat_width and mat[i][j] != -1:
            sum += mat[i][j] # sum of values of adjacent pixels ++
            nbr_adj_pixel += 1 # nbr of adjacent pixels ++
        return sum, nbr_adj_pixel
    
    #process the missing pixels
    filtered_mat = [row[:] for row in mat] #deep copy of mat
    for i in range(mat_height):
        for j in range(mat_width):
            if mat[i][j] == -1: # is a missing pixel
                sum, nbr_adj_pixel = 0, 0
                sum, nbr_adj_pixel = check_adjacent_pixel(i-1, j, sum, nbr_adj_pixel) # up pixel
                sum, nbr_adj_pixel = check_adjacent_pixel(i+1, j, sum, nbr_adj_pixel) # bottom pixel
                sum, nbr_adj_pixel = check_adjacent_pixel(i, j-1, sum, nbr_adj_pixel) # left pixel
                sum, nbr_adj_pixel = check_adjacent_pixel(i, j+1, sum, nbr_adj_pixel) # right pixel
                
                if nbr_adj_pixel > 0:
                    filtered_mat[i][j] = round(float(sum/nbr_adj_pixel)) #the missing pixels are replaced by the mean of adjacent pixels
                else:
                    filtered_mat[i][j] = 0

    return filtered_mat


def sequence_to_rgb_mat(sequence: str, mat_height: int, mat_width: int) -> list:
    """
    convert a sequence to a matrix representing a rgb image
    the given sequence is supposed to have the same height and the width *3 of the original image
    """
    mat = [ [] for k in range(mat_height) ]
    sequence_rows = [sequence[i:i+mat_width*3] for i in range(0, len(sequence), mat_width*3)]
    
    for i in range(mat_height):
        seq_row = sequence_rows[i]
        for j in range(mat_width):
            kmer = seq_row[3*j:3*(j+1)]
            pxl_value_list = [-1,-1,-1] #count as a missing pixel if the kmer is not recognized
            for k in range(3):
                for key, value in base_4_to_dna_dict.items():
                    if value == kmer[k]:
                        rounded_pxl_value = key
                        pxl_value = 85*rounded_pxl_value #from 0-3 to 0-255
                        pxl_value_list[k] = pxl_value
                        break
            mat[i].append(pxl_value_list)
            
    def check_adjacent_pixel(i, j, sum, nbr_adj_pixel):
        # if pixel in matrix and not missing
        if i >= 0 and i < mat_height and j >= 0 and j < mat_width and mat[i][j][k] != -1:
            sum += mat[i][j][k] # sum of values of adjacent pixels ++
            nbr_adj_pixel += 1 # nbr of adjacent pixels ++
        return sum, nbr_adj_pixel
    
     #process the missing pixels
    filtered_mat = [[line[:] for line in row] for row in mat] #deep copy of mat
    for i in range(mat_height):
        for j in range(mat_width):
            for k in range(3):
                if mat[i][j][k] == -1: # is a missing pixel
                    sum, nbr_adj_pixel = 0, 0
                    sum, nbr_adj_pixel = check_adjacent_pixel(i-1, j, sum, nbr_adj_pixel) # up pixel
                    sum, nbr_adj_pixel = check_adjacent_pixel(i+1, j, sum, nbr_adj_pixel) # bottom pixel
                    sum, nbr_adj_pixel = check_adjacent_pixel(i, j-1, sum, nbr_adj_pixel) # left pixel
                    sum, nbr_adj_pixel = check_adjacent_pixel(i, j+1, sum, nbr_adj_pixel) # right pixel
    
                    if nbr_adj_pixel > 0:
                        filtered_mat[i][j][k] = round(float(sum/nbr_adj_pixel)) #the missing pixels are replaced by the mean of adjacent pixels
                    else:
                        filtered_mat[i][j][k] = 0
    return filtered_mat   
  
    
def matrix_to_png(mat: np.array, mode: str, image_path: str) -> None:
    """
    save the png image from the matrix
    """
    np_mat = np.asarray(mat, dtype=np.uint8)
    img = Image.fromarray(np_mat, mode)
    img.save(image_path)
  

def resize_img(image_path: str, result_path: str, height: int, width: int) -> None:
    """
    resize the image
    """
    img = Image.open(image_path)
    img = img.resize((height,width),Image.ANTIALIAS)
    img.save(result_path)         


if __name__ == "__main__":
    """mat = png_to_matrix("img_small.png", 'RGB')
    seq, meta = rgb_mat_to_sequence(mat)
    s_seq = shuffle_sequence(seq, 3)
    u_seq = unshuffle_sequence(s_seq, len(seq), 3)
    mat2 = sequence_to_rgb_mat(seq, 118, 106)
    for l in mat2:
        for k in l:
           # print(k)
            pass
    matrix_to_png(mat2, 'RGB', "img_decoded.png")
    """
    sequence = sys.argv[1]
    decode_png(sequence, "test.png", "RGB;118;106")

