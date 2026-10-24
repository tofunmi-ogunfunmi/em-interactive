# Interactive electron microscopy

Browser simulations of three electron-microscopy techniques used to study solid-state battery materials:

- **FIB-SEM tomography**: mill a solid electrolyte slice by slice, image each face and build a 3D reconstruction of lithium filaments, voids and particles.
- **4D-STEM**: scan a probe over a polycrystalline sample, record a diffraction pattern at every position, and form virtual bright-field, dark-field and orientation images. Includes **multislice simulations of the LLZO garnet electrolyte** computed with abTEM.
- **EELS**: move a probe across a Li metal | interphase | LLZO interface, fit and subtract the background under each edge, and acquire elemental maps.

**Live page:** https://tofunmi-ogunfunmi.github.io/em-interactives/

By [Tofunmi Ogunfunmi](https://tofunmi-ogunfunmi.github.io/), UC Berkeley.

## Contents

| Path | What it is |
| --- | --- |
| `index.html` | The page |
| `js/em-sims.js` | All simulations and viewers (plain JavaScript, no dependencies) |
| `js/llzo-4dstem-data.js` | Precomputed abTEM 4D-STEM data for LLZO (~6 MB, loaded only when that demo is on screen) |
| `simulations/` | Python scripts that generate the LLZO data, and the input structure |

## Running locally

Open `index.html` in a browser, keeping the `js` folder next to it. Or serve the folder:

```
python3 -m http.server
```

## Reproducing the LLZO data

```
pip install abtem spglib
cd simulations
python build_cell.py   # primitive mp-942733 cell -> conventional tetragonal cell (llzo_conv.cif)
python sim.py          # multislice 4D-STEM at 1, 6 and 20 mrad (about 5 min on 2 CPU cores)
python pack.py         # crop, log-encode and compress -> llzo-4dstem-data.js
cp llzo-4dstem-data.js ../js/
```

Settings: tetragonal Li7La3Zr2O12 (I4₁/acd) along [001], 8 unit cells (10.2 nm) thick, 300 kV, Lobato potentials,
2.1 Å slices, 24 × 24 probe positions per unit cell, no aberrations, no thermal diffuse scattering.

## Models and limitations

These are teaching tools, not quantitative simulations.

- The FIB-SEM, toy 4D-STEM and EELS demos use simplified models: an analytic 3D volume with approximate SEM contrast; a single-slice phase-object diffraction model; and schematic EELS edge shapes with Poisson plasmon scattering.
- The LLZO data are full multislice calculations, but use a DFT-relaxed cell (about 1% larger than experiment) and omit thermal diffuse scattering, so dark-field contrast is lower than in experiments.

## Credits

- Inspired by the interactive 4D-STEM demo from the [Colin Ophus Lab (COLab), Stanford](https://colab.stanford.edu/interactive) and the [py4DSTEM](https://github.com/py4dstem/py4DSTEM) toolkit. Code here is written independently.
- LLZO simulations: J. Madsen and T. Susi, "The abTEM code: transmission electron microscopy from first principles," *Open Research Europe* 1, 24 (2021). [doi:10.12688/openreseurope.13015.1](https://doi.org/10.12688/openreseurope.13015.1)
- LLZO structure: [Materials Project mp-942733](https://next-gen.materialsproject.org/materials/mp-942733) (CC BY 4.0). A. Jain et al., *APL Materials* 1, 011002 (2013).

## License

Code: MIT (see `LICENSE`). Simulated data: CC BY 4.0, please cite this repository and the sources above.
