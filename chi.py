import numpy as np

Q_in = '''
111101011001100101000001
011110101100110010100001
001111010110011001010001
000111101011001100101001
000011110101100110010101
000001111010110011001011
100000111101011001100101
010000011110101100110011
101000001111010110011001
010100000111101011001101
001010000011110101100111
100101000001111010110011
'''.replace('\n','')


G_in = '''
100000000000011100100111
010000000000101010010111
001000000000110001001111
000100000000111011100100
000010000000111101010010
000001000000111110001001
000000100000100111011100
000000010000010111101010
000000001000001111110001
000000000100100100111011
000000000010010010111101
000000000001001001111110
'''.replace('\n', '')

docstr_to_mat = lambda s: np.array([np.array([int(s[i*24 + j]) for j in range(24)]) for i in range(12)])

Q = docstr_to_mat(Q_in)
G = docstr_to_mat(G_in)

def get_perms(arr,seen,perm,perms,i=0):
    if i == len(perm):
        perms.append(perm)
        return
    for n in arr:
        if n not in seen:
            seen.add(n)
            perm[i] = (arr[i],n)
            get_perms(arr,seen,perm[:],perms,i+1)
            seen.remove(n)


def generate_perms():
    perms = []
    F = [1,2,3,4,5,18,21]
    C = [6,15,17,20,23]
    D = [7,8,11,12,13,14,19,22]
    fixed = [(a,a) for a in F] + [(a,b) for a,b in zip([9,10,16,24],[10,9,24,16])]

    cperms, dperms = [], []
    get_perms(C, set(), [None]*len(C), cperms)
    get_perms(D, set(), [None]*len(D), dperms)

    for cperm in cperms:
        for dperm in dperms:
            perm = [0] * 24
            for from_coord, to_coord in fixed + cperm + dperm:
                perm[to_coord - 1] = from_coord - 1
            perms.append(perm)
    return perms

def find_chi(G_columns, Q_rows):
    for perm in generate_perms():
        G_columns_perm = G_columns[perm]           # permutes each basis *column* of G
        mat_prod = Q_rows.dot(G_columns_perm) % 2  # matrix product in F^2
        if np.array_equal(mat_prod, np.zeros((12,12))):
            print(perm)
         

find_chi(G.transpose(),Q)
print('possibly not found')










