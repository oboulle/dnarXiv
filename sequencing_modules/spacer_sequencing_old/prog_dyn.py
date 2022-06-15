
# programmation dynamique
def compareFrag(f1,f2,T):
    f1 = '*'+f1
    f2 = '*'+f2
    l1 = len(f1)
    l2 = len(f2)
    M = [[0 for i in range(l2)] for j in range(l1)]
    for i in range(l1):
        M[i][0] = i
    for j in range(l2):
        M[0][j] = j
    for i in range(1,l1):
        m = 1000
        for j in range(1,l2):
            if f1[i]==f2[j] :
                sub = 0
            else:
                sub = 1
            M[i][j] = min(M[i-1][j-1]+sub,M[i-1][j]+1,M[i][j-1]+1)
            m = min(M[i][j],m)
        if m > T:
            return m
    return M[l1-1][l2-1]

    a1 = ""
    a2 = ""
    am = ""
    i = l1-1
    j = l2-1
    while i > 0 and j > 0:
        if f1[i]==f2[j]:
            sub = 0
        else:
            sub = 1
        
        if M[i-1][j-1]+sub <= M[i][j-1]+1 and M[i-1][j-1]+sub <= M[i-1][j]+1 :
            a1 = f1[i] + a1
            a2 = f2[j] + a2
            if sub == 0 :
                am = '|'+am
            else:
                am = ' '+am
            i -= 1
            j -= 1
        else:
            am = ' ' + am
            if M[i][j-1] > M[i-1][j]:
                a1 = f1[i] + a1
                a2 = '-' + a2
                i -= 1
            else:
                a1 = '-' + a1
                a2 = f2[j] + a2
                j -= 1

    return M[l1-1][l2-1],a1,a2,am
