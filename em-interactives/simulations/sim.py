import time, sys, ase.io, numpy as np, abtem
c = ase.io.read('llzo_conv.cif'); a, cz = c.cell[0,0], c.cell[2,2]
NS = 24          # probe positions per unit-cell edge
OUT = 64         # stored pattern size (pixels)
# name, lateral cells, gpts, semiangle (mrad), crop (px, centered), bin
PRESETS = [('nbd', 8, 512, 1.0, 128, 2), ('mid', 4, 256, 6.0, 128, 2), ('atom', 4, 384, 20.0, 192, 3)]
which = sys.argv[1:] or [p[0] for p in PRESETS]
for name, nxy, gpts, alpha, crop, b in PRESETS:
    if name not in which: continue
    atoms = c * (nxy, nxy, 8)
    pot = abtem.Potential(atoms, gpts=gpts, slice_thickness=cz/6, projection='infinite', parametrization='lobato').build().compute(progress_bar=False)
    probe = abtem.Probe(energy=300e3, semiangle_cutoff=alpha, defocus=0)
    data = np.zeros((NS, NS, OUT, OUT), np.float32); t0 = time.time()
    for iy in range(NS):
        y = iy * a / NS
        scan = abtem.LineScan(start=(0, y), end=(a, y), gpts=NS, endpoint=False)
        m = probe.scan(pot, scan=scan, detectors=abtem.PixelatedDetector(max_angle=None)).compute(progress_bar=False)
        arr = np.asarray(m.array)                          # (NS, gpts, gpts), fft-shifted
        n = arr.shape[-1]; h = n // 2; s = crop // 2
        if iy == 0:
            pk = np.unravel_index(arr[0].argmax(), arr[0].shape); print(name, 'max at', pk, 'centre', h, 'k-sampling', m.sampling, flush=True)
        cr = arr[:, h - s:h + s, h - s:h + s]
        data[:, iy] = cr.reshape(NS, OUT, b, OUT, b).sum(axis=(2, 4))   # data[ix, iy]
        print(name, f'row {iy+1}/{NS}  {time.time()-t0:.0f}s', flush=True)
    dk = m.sampling[-1] * b
    np.savez_compressed(f'{name}.npz', data=data.transpose(1, 0, 2, 3), dk=dk, alpha=alpha, a=a, step=a/NS)  # stored [iy, ix, ky, kx]
    print(name, 'done; binned dk', dk, flush=True)
