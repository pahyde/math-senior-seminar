import numpy as np
from itertools import permutations
import pprint as pp

Q_in = '''
111101011001100101000001011110101100110010100001001111010110011001010001000111101011001100101001
000011110101100110010101000001111010110011001011100000111101011001100101010000011110101100110011
101000001111010110011001010100000111101011001101001010000011110101100111100101000001111010110011
'''.replace('\n','')

G_in = '''
100000000000011100100111010000000000101010010111001000000000110001001111000100000000111011100100
000010000000111101010010000001000000111110001001000000100000100111011100000000010000010111101010
000000001000001111110001000000000100100100111011000000000010010010111101000000000001001001111110
'''.replace('\n', '')

def generate_perms():
    perms = []
    F,C,D = [1,2,3,4,5,18,21], [6,15,17,20,23], [7,8,11,12,13,14,19,22]
    fixed = [(a,a) for a in F] + [(a,b) for a,b in zip([9,10,16,24],[10,9,24,16])]
    cperms = [[(C[i], perm[i]) for i in range(len(perm))] for perm in permutations(C)]
    dperms = [[(D[i], perm[i]) for i in range(len(perm))] for perm in permutations(D)]
    for cperm in cperms:
        for dperm in dperms:
            perm = [0] * 24
            for from_coord, to_coord in fixed + cperm + dperm:
                perm[to_coord - 1] = from_coord - 1
            perms.append(perm)
    return perms

docstr_to_mat = lambda s: np.array([np.array([int(s[i*24 + j]) for j in range(24)]) for i in range(12)])
Q = docstr_to_mat(Q_in)              # row vector basis for R
G = docstr_to_mat(G_in).transpose()  # col vector basis for B


G = G.transpose()


def get_words(gen_mat):
    words = []
    for mask in range(1 << 12):
        word = np.array([0] * 24)
        for i in range(12):
            if mask & 1:
                word = (word + gen_mat[i]) % 2
            mask >>= 1
        if word.dot(word) == 8:
            words.append(word)
    return words

B = get_words(G)
R = get_words(Q)

weight = lambda word: word.dot(word)
C = [i-1 for i in [6,9,10,15,16,17,20,21,23,24]]
D = [i-1 for i in [7,8,11,12,13,14,19,22]]
def find_ones(positions, code, sset):
    return [i+1 for word in code for i,b in enumerate(word) if i in sset and b == 1 if np.array_equal(word[[0,1,2,3,4,17,20]], np.array(positions)) and weight(word[sset]) == 1]



class Partition:
    def __init__(self, indices, size, target_bits = None, target_weight = None):
        self.idxs = indices
        self.size = size
        self._target_bits = target_bits
        self._target_weight = target_weight

    def indices(self):
        return self.idxs
    
    def target_bits(self, bits):
        return Partition(self.idxs, self.size, target_bits = bits, target_weight = self._target_weight)

    def target_weight(self, weight):
        return Partition(self.idxs, self.size, target_bits = self._target_bits, target_weight = weight)

    def satisfies_target_bits(self, word):
        if self._target_bits == None:
            return True
        return np.array_equal(word[self.idxs], np.array(self._target_bits))

    def satisfies_target_weight(self, word):
        if self._target_weight == None:
            return True
        return weight(word[self.idxs]) == self._target_weight



class Golay_Model:
    def __init__(self, generator_matrix, codewords = None):
        if codewords == None:
            codewords = self._codewords(generator_matrix)
        self.mat       = generator_matrix
        self.codewords = codewords

    def get_words(self):
        return [list(word) for word in self.codewords]
    
    def get_generator_mat(self):
        return self.mat

    def permute(self, permutation):
        return Golay_Model(self.mat[:,permutation], [word[permutation] for word in self.codewords])

    def is_equivalent_to(self, golay_model):
        mat_prod = self.mat.dot(golay_model.get_generator_mat().transpose()) % 2
        return np.array_equal(mat_prod, np.zeros((12,12)))

    def filter(self, partition):
        filtered_words = [] 
        for word in self.codewords:
            if not partition.satisfies_target_bits(word):
                continue
            if not partition.satisfies_target_weight(word):
                continue
            filtered_words.append(word)
        return Golay_Model(self.mat, filtered_words)

    def __len__(self):
        return len(self.codewords)

    def _codewords(self, gen_mat):
        words = []
        for mask in range(1 << 12):
            word = np.array([0] * 24)
            for i in range(12):
                if mask & 1:
                    word = (word + gen_mat[i]) % 2
                mask >>= 1
            words.append(word)
        return words

class Infinite_Set:
    def __init__(self, elements = set(), is_infinite = True):
        self.elements = elements
        self.is_infinite = is_infinite

    def __or__(self, other):
        if self.is_infinite:
            return Infinite_Set()
        return Infinite_Set(self.elements | other)

    def __and__(self, other):
        if self.is_infinite:
            return Infinite_Set(other, False)
        return Infinite_Set(self.elements & other, False)

    def __str__(self):
        if self.is_infinte:
            return 'inf'
        return str(self.elements)
    
    def __repr__(self):
        if self.is_infinite:
            return 'inf'
        return repr(self.elements)  

    def __eq__(self, other):
        return self.elements == other.elements

    def __len__(self):
        if self.is_infinite:
            return np.inf
        return len(self.elements)

    def get_elements(self):
        return self.elements
        


indices = [i for i in range(24)]
# We need Golay codewords with weight 8 to form a 5 - (24,8,1) block design
R = Golay_Model(Q).filter(Partition(indices, 24).target_weight(8))
B = Golay_Model(G).filter(Partition(indices, 24).target_weight(8))

# Defining the relavent partitions
fixed = Partition([i-1 for i in [1,2,3,4,5,18,21]], 24)
C     = Partition([i-1 for i in [6,9,10,15,16,17,20,21,23,24]], 24)
D     = Partition([i-1 for i in [7,8,11,12,13,14,19,22]], 24)

# Function used to set fixed bits
# returns possible tuples of the from (x1,x2,x3,...,xl)
# with exactly 4 bits 1 and the rest 0
def get_four_tuples(length):
    four_tuples = []
    for i in range(1 << length):
        bits = []
        for _ in range(length):
            bits.append(1 if i & 1 else 0)
            i >>= 1
        if sum(1 for b in bits if b == 1) == 4:
            four_tuples.append(bits)
    return four_tuples
            

# from_to_map. maps individual indices in B to a set of candidate indices in R
# initially each index in B can map to anything.
from_to_map = {i: Infinite_Set() for i in range(24)}

# All fixed indices must map back to themselves
for i in fixed.indices():
    from_to_map[i] &= set([i])

# We conjecture that 23 -> 15 and 20 -> 20
from_to_map[23] &= {15}
from_to_map[20] &= {20}

# We seek corresponding codeword in B and R
# which have exacty 4 one bits in fixed possitions, all other 0
for four_tuple in get_four_tuples(7):
    # For each partition C and D, we look for
    # codewords that have a single bit in that partition. 
    # We insist these words also satisfy the four tuple requirement mentioned above
    for partition in [C,D]:
        B_words = B.filter(fixed.target_bits(four_tuple)).filter(partition.target_weight(1)).get_words()
        R_words = R.filter(fixed.target_bits(four_tuple)).filter(partition.target_weight(1)).get_words()

        # Proceed if B and R actually have codewords meeting these requirements
        if len(B_words) != 0: 
            # indices in words of B that meet this requirements must map to those of R
            B_from_idxs = [i for word in B_words for i in partition.indices() if word[i] == 1]
            R_to_idxs   = [i for word in R_words for i in partition.indices() if word[i] == 1]
            # iteratively narrow the possible R indices which a give B index can map to
            # accomplished by set intersection
            for from_idx in B_from_idxs:
                from_to_map[from_idx] &= set(R_to_idxs)

# collect pairs function. Explanation below.
def collect_pairs(from_to):
    reduced = []
    for i, (from_idx, to_idxs) in enumerate(from_to.items()):
        if len(to_idxs) == 1:
            reduced.append([from_idx, list(to_idxs.get_elements())[0]])
        else:
            for j in range(i+1,24):
                if to_idxs == from_to[j]:
                    reduced.append([(i,j), tuple(to_idxs.get_elements())])
    return reduced
                
# Upon inspection, each index in B must map to either one or two indices in R
# If there is a index if B, b1, which maps to a pair in R, (r1,r2),
# then there is another b2 which maps to (r1, r2). 
# Thus we are either absolutly certain about a mapping b -> r 
# or we have two choices for a situation such as {b1, b2} -> {r1, r2}

# collect pairs such as b1 -> (r1, r2) and b2 -> (r1, r2) into (b1,b2) -> (r1,r2)
reduced_from_to_map = collect_pairs(from_to_map)

# deterimined maps: b -> r. undetermined maps: (b1, b2) -> (r1, r2).
determined   = [mapping for mapping in reduced_from_to_map if not isinstance(mapping[0], tuple)]
undetermined = [mapping for mapping in reduced_from_to_map if isinstance(mapping[0], tuple)]

def possible_mappings(from_to_pairs):
    mappings = []
    for mask in range(1 << len(from_to_pairs)):
        mapping = []
        for (b1,b2), (r1,r2) in from_to_pairs:
            if mask & 1:
                mapping.append([b1,r1])
                mapping.append([b2,r2])
            else:
                mapping.append([b1,r2])
                mapping.append([b2,r1])
            mask >>= 1
        mappings.append(mapping)
    return mappings


def make_permutation_array(idx_map):
    perm_array = [None] * len(idx_map)
    for from_idx, to_idx in idx_map:
        perm_array[to_idx] = from_idx
    return perm_array

equivalence_maps = []

for pair_selection in possible_mappings(undetermined):
    candidate_permutation = make_permutation_array(determined + pair_selection)
    image_B = B.permute(candidate_permutation)
    if R.is_equivalent_to(image_B):
        equivalence_maps.append(candidate_permutation)


def cycle_permutation(permutation_array):
    seen = set()
    perm_map = {from_idx+1: to_idx+1 for to_idx, from_idx in enumerate(permutation_array)}
    cycles = []
    for i in range(1,24+1):
        if perm_map[i] == i:
            continue
        if i not in seen:
            seen.add(i)
            cycle = [i]
            while perm_map[i] not in seen:
                i = perm_map[i]
                seen.add(i)
                cycle.append(i)
            cycles.append(cycle)
    return ' '.join('({})'.format(', '.join(map(str,cycle))) for cycle in cycles)

print('equivalence maps chi satisfying out constraints')
for equivalence_map in equivalence_maps:
    print(cycle_permutation(equivalence_map))

sigma = np.array([22,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23])
chi = np.array(equivalence_maps[1])
chi_inv = [0] * 24
for ti, fi in enumerate(chi):
    chi_inv[fi] = ti
chi_inv = np.array(chi_inv)
tau = np.array([i for i in range(24)])[chi][sigma][chi_inv]
tau = tuple(tau)

chip = np.array(equivalence_maps[0])
chi_invp = [0] * 24
for ti, fi in enumerate(chi):
    chi_invp[fi] = ti
chi_invp = np.array(chi_invp)
taup = np.array([i for i in range(24)])[chip][sigma][chi_invp]
taup = tuple(taup)

print('tau')
print(cycle_permutation(tau))

rho = tuple([j + 3*i for i in range(8) for j in [2,0,1]])

        
def find_unique(unique, rho, tau, perm, tau_prime):
    stack = [perm]
    stage = [0]
    sp = 0
    i = 0
    while True:
        if sp == -1:
            return
        if stack[sp] == taup or len(unique) >= 244823040:
            break
        if (stage[sp] == 0 and stack[sp] in unique) or stage[sp] == 3:
            stack.pop()
            sp -= 1
            continue
        if stage[sp] == 0:
            unique.add(stack[sp])
            stage[sp] += 1
            continue
        if stage[sp] == 1:
            tau_perm = tuple([stack[sp][t] for t in tau])
            stack.append(tau_perm)
            stage[sp] += 1
            stage.append(0)
            sp += 1
        if stage[sp] == 2:
            rho_perm = tuple([stack[sp][r] for r in rho])
            stack.append(rho_perm)
            stage[sp] += 1
            stage.append(0)
            sp += 1

def find_uniquer(unique, rho, tau, perm, tau_prime, depth=1):
    print(len(unique))
    if depth > 158:
        return
    if perm in unique:
        return
    if perm == tau_prime:
        return 
    unique.add(perm)
    depth += 1
    find_uniquer(unique, rho, tau, tuple([perm[r] for r in rho]), tau_prime, depth)
    find_uniquer(unique, rho, tau, tuple([perm[t] for t in tau]), tau_prime, depth)

print('finding size')

test = [0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1]
test = np.array(test)

print(list(test[chi]))
print(list(test))
#unique = set()
#find_uniquer(unique, rho, tau, tuple(i for i in range(24)), taup)
#print(len(unique))

    
    

    


'''
print(len(B))
collect = [set() for _ in range(24)]
for i in range(7):
    for j in range(i+1,7):
        for k in range(j+1,7):
            pos = [1]*7
            pos[i] = 0
            pos[j] = 0
            pos[k] = 0
            for sset in [C,D]:
                bwords = find_ones(pos, B, sset)
                rwords = find_ones(pos, R, sset) 
                for b_idx in bwords:
                    if not len(collect[b_idx-1]):
                        collect[b_idx-1] |= set(rwords)
                    else:
                        collect[b_idx-1] &= set(rwords)

perm_map = [(i+1,collect[i]) for i in range(len(collect))]
unset_idxs = [idxs for idxs in perm_map if len(idxs[1]) > 1]
pairs = [(set(f for f,t in unset_idxs if t == to_idxs), to_idxs) for from_idx, to_idxs in unset_idxs]
print(len(unset_idxs))
print(pairs)
                        ''' 


'''

perms = generate_perms()
for chi in perms:
    G_perm   = G[chi]                                # permutes each basis *column* of G
    mat_prod = Q.dot(G_perm) % 2                     # matrix product in F^2
    if np.array_equal(mat_prod, np.zeros((12,12))):  # chi is an equivalence map if Q(G^T) = 0
        print("succes: ")
        print(chi)

'''











