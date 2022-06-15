#scoring constants
score_match = 10
score_gap = -7
score_mismatch = -5


def local_align(x: str,y: str) -> (float, list):
    '''
    Smith-Waterman algorithm
    Fill in the matrix with alignment scores and obtain the best score and position
    worst possible alignment : score=0
    perfect alignment : score=score_match
    :param x: the long sequence
    :param y: the short sequence to align on x
    :return: the score and the coordinates [start_index;stop_index] of the best local alignment found
    '''
    
    #filling of the scoring matrix
    matrix = [ [ 0 for j in range(len(x)+1) ] for i in range(len(y)+1) ] #init empty matrix 

    best_score = 0
    best_score_index = [0, 0]
    
    for i in range(1,len(y)+1):
        for j in range(1,len(x)+1):
            matrix[i][j] = max(
                matrix[i][j-1] + score_gap,
                matrix[i-1][j] + score_gap,
                matrix[i-1][j-1] + (score_match if x[j-1] == y[i-1] else score_mismatch),
                0
                )

            if matrix[i][j] >= best_score:
                best_score = matrix[i][j]
                best_score_index = [i,j]
                
                
    #reconstruction of the aligned sequence
    alignment = ''

    #coordinates of the best score in the matrix
    i = best_score_index[0]
    j = best_score_index[1]

    while matrix[i][j] > 0:
        diag = matrix[i-1][j-1]
        up = matrix[i-1][j]
        left = matrix[i][j-1]
        if matrix[i][j] == diag + (score_match if x[j-1] == y[i-1] else score_mismatch):
            i = i - 1
            j = j - 1
            alignment += x[j]
        elif matrix[i][j] == left + score_gap:
            j = j - 1
            alignment += '_'
        elif matrix[i][j] == up + score_gap:
            i = i - 1
        else:
            #when matrix[i][j] = 0
            break
    #print(alignment[::-1]) #reconstructed alignment, unused
    adjusted_score = best_score/len(y) #score is adjusted to depend on the length of the sequence to align
    
    #at the end of the while loop, j is the starting index of the aligned sequence, best_score_index[1] is the stop index
    return adjusted_score, [j, best_score_index[1]]

# =================== main ======================= #
if __name__ == '__main__':
    x = 'CTGGATAAATAATGACCGCTTCAGGTAGCATCTACGTATTTCGGTGGTCAGTACACTGAA'#AGTAGGAATCATGCTCGATCGTTAACCTCCAATAGACCTGAGGTAAACTCACGATTCTCCTAGAAATGCTAACATACTCCTCGCGATAAGCTATACCCA'
    y = 'ACTTCGTTCAGTTACGTATTGC'
    print('Input sequences are: ')
    print(x)
    print(y)
    print(local_align(x,y))
