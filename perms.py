
def eta(arr):
    return [arr[i] for i in [4,0,1,2,7,3,6,5]]

def epsilon(arr):
    return [arr[i] for i in [3,0,1,2,7,4,5,6]]


def find(curr, seen):
    key = ','.join(map(str,curr))
    if key in seen:
        return
    seen.add(key)
    find(eta(curr), seen)
    find(epsilon(curr), seen)

seen = set()
find([1,2,3,4,5,6,7,8], seen)

print(len(seen)) 


