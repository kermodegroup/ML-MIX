import ase
import ase.io
from ase.build import bulk

def build_struct(elements, n=10, rank=0,path='./',structure='fcc',a=3.0):
    struct = bulk('H', structure, a=a, cubic=True)
    struct *= n

    elements_str = '-'.join(elements)

    # randomise atoms to different elements
    nelements = len(elements)
    for i, atom in enumerate(struct):
        element = elements[i % nelements]
        atom.symbol = element

    struct.rattle(0.1)
    if rank == 0:
        ase.io.write(f'{path}/{elements_str}.xyz', struct, parallel=False)

    return struct




