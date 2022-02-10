

rows = [
'111101011001100101000001',
'011110101100110010100001',
'001111010110011001010001',
'000111101011001100101001',
'000011110101100110010101',
'000001111010110011001011',
'100000111101011001100101',
'010000011110101100110011',
'101000001111010110011001',
'010100000111101011001101',
'001010000011110101100111',
'100101000001111010110011'
]

rows = [[int(n) for n in list(row)] for row in rows]

xor = lambda a,b: [(m+n) % 2 for m,n in zip(a,b)]

def perm(v):
    res = [None] * len(v)
    for i in range(23):
        res[(i+1) % 23] = v[i]
    res[23] = v[22]
    return res

target = perm([1,0,0,1,0,1,0,0,0,0,0,1,1,1,1,0,1,0,1,1,0,0,1,1])

def choose(r,b):
    res = []
    for i in range(12):
        if b & 1:
            res.append((i,r[i]))
        b >>= 1
    return res

combs = [choose(rows,i) for i in range(1 << (12))]

for comb in combs:
    xord = [0] * 24
    idxs = []
    for i,c in comb:
        xord = xor(xord, c)
        idxs.append(i)
    if all(xord[i] == target[i] for i in range(24)):
        print(comb)
        print(idxs)
