import ase.io, spglib, numpy as np
from ase import Atoms
p = ase.io.read('llzo_mp-942733.cif')
print('primitive', len(p), p.get_chemical_formula(), 'vol', p.get_volume())
cell = (p.cell[:], p.get_scaled_positions(), p.numbers)
print('spacegroup', spglib.get_spacegroup(cell, symprec=1e-3))
lat, pos, num = spglib.standardize_cell(cell, to_primitive=False, no_idealize=False, symprec=1e-3)
c = Atoms(numbers=num, scaled_positions=pos, cell=lat, pbc=True)
print('conventional', len(c), c.get_chemical_formula(), np.round(c.cell.cellpar(), 4))
ase.io.write('llzo_conv.cif', c)
