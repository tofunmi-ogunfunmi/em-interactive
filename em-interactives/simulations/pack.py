import numpy as np, ase.io, gzip, base64, json
c = ase.io.read('llzo_conv.cif'); a = c.cell[0,0]; cz = c.cell[2,2]
DEC = 8.0   # decades of dynamic range in the log encoding
META = {'nbd': (1.0, 128, 2), 'mid': (6.0, 128, 2), 'atom': (20.0, 192, 3)}
LAM = 0.019687  # Å at 300 kV
presets = []
for name in ['nbd', 'mid', 'atom']:
    z = np.load(f'{name}.npz'); d = z['data'].transpose(0, 1, 3, 2)   # -> [iy, ix, ky, kx]
    mrad, crop, b = META[name]; dk = float(z['dk']); s = crop // 2; c0 = (s - (b - 1) / 2) / b
    M = float(d.max()); r = d / M
    v = np.where(r > 10**-DEC, np.clip(np.round(255 * (np.log10(np.maximum(r, 1e-30)) / DEC + 1)), 1, 255), 0).astype(np.uint8)
    blob = gzip.compress(v.tobytes(), 9)
    # physics check on the position-averaged pattern: allowed (hk0) need h,k even in I4_1/acd
    avg = d.mean(axis=(0, 1)); n = avg.shape[0]
    def at(h, k):
        x = c0 + h / a / dk; y = c0 + k / a / dk; xi, yi = int(round(x)), int(round(y))
        return avg[yi-1:yi+2, xi-1:xi+2].max() / avg.max() if 1 <= xi < n-1 and 1 <= yi < n-1 else float('nan')
    print(f'{name}: alpha={mrad} mrad ({mrad/1000/LAM:.3f} 1/Å), binned dk={dk:.4f}, c0={c0:.2f}, gz={len(blob)/1e6:.2f} MB')
    print('   (2,0,0)=%.2e (1,0,0)=%.2e (2,2,0)=%.2e (1,1,0)=%.2e (4,0,0)=%.2e' % (at(2,0), at(1,0), at(2,2), at(1,1), at(4,0)))
    presets.append({'id': name, 'mrad': mrad, 'alpha': mrad / 1000 / LAM, 'dk': dk, 'c0': c0, 'M': M, 'decades': DEC,
                    'b64': base64.b64encode(blob).decode()})
# projected atom columns in one unit cell (Å), x right, y down (same as scan)
atoms = {}
for el in ['La', 'Zr', 'O', 'Li']:
    xy = c.positions[np.array(c.get_chemical_symbols()) == el][:, :2] % a
    u = np.unique(np.round(xy / 0.05).astype(int), axis=0) * 0.05
    atoms[el] = [round(float(q), 3) for q in u.ravel()]
    print(el, 'columns', len(u))
out = {'a': a, 'c': cz, 'thickness': 8 * cz, 'ns': 24, 'n': 64, 'atoms': atoms, 'presets': presets,
       'source': 'Li7La3Zr2O12, tetragonal I4_1/acd, Materials Project mp-942733; abTEM multislice, 300 kV, Lobato potentials, no thermal diffuse scattering'}
js = '/* Simulated 4D-STEM data for LLZO [001] (abTEM). Generated file; see em-sims notes. */\nwindow.LLZO_4DSTEM=' + json.dumps(out, separators=(',', ':')) + ';\n'
open('llzo-4dstem-data.js', 'w').write(js); print('data file MB', len(js) / 1e6)
