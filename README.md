# MHD viztools Fortran — visualization tools for planetary magnetosphere MHD simulations

[日本語](README.ja.md)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24030429.svg)](https://doi.org/10.5281/zenodo.24030429)

A set of Fortran tools to combine, analyse and visualize the output of global MHD simulations of planetary
magnetospheres. **They are a modernized version of the visualization tools developed and used by the late Prof. Tatsuki Ogino
for his magnetospheric MHD codes** (slice plots, polar ionosphere plots, 3-D display with VRML, etc.), rewritten to be
easy to use on present-day computers. They were developed and tested with the output of the Jupiter code (jupiter7),
and work with Saturn and Earth simulations that use the same output format.

## Features

* Fortran 2008 + OpenMP. External libraries are needed only for HDF5 output
* Everything is set in namelists (no recompilation when changing grids, files, quantities or slices)
* Huge data (standard grid 3002×2002×1002 × 8 variables, about 190 GB) are handled with a few GB of memory, because only
  the planes, rows and variables needed are read directly from the files
* Normalized quantities are converted back to physical units (cm⁻³, km/s, nPa, nT, eV, µA/m², kV, ...). Normalization
  constants of Jupiter, Saturn and Earth are built in, and the distance unit (RJ / RS / RE) follows the planet
* Derived quantities (temperature p/ρ etc.), unit conversions and differences are computed in double precision
* 3-D scenes are written as glTF 2.0 (.glb) and X3D, viewable in a browser (the included viewer.html), Blender and ParaView

## Tools

| tool | purpose |
|---|---|
| `combine3d` | combine per-MPI-rank output into one 3-D data file (`m.data`) |
| `vorticity` | vorticity components parallel and perpendicular to B |
| `bin2hdf` | convert to HDF5 + XDMF (or raw binary + XDMF) for ParaView / Python |
| `slice_plot` | colour plots of any slice (xy / xz / yz at any position) with vectors (PostScript) |
| `polar_map` | mapping to the polar ionosphere: field-aligned currents, convection and potential, open-closed boundary, energy flux |
| `mag3d` | 3-D visualization: field lines (connected to the planet, solar wind, loops), slices, stacked slices, selectable slices, isosurfaces, arrows, planet images, colour bars |

The output formats of the original tools (`m.data`, `.vort`, `.h5`/`.xmf`) are kept, so existing data and analysis scripts still work.
Problems found in the original tools (mapping and potential solver of the polar plots, etc.) and how they were fixed are described in the [user manual](doc/manual_en.md).

## Requirements

* Fortran compiler: gfortran, Intel ifx/ifort, NEC nfort, Fujitsu frtpx, NVIDIA nvfortran (selected in the Makefile)
* OpenMP (optional)
* HDF5 with the Fortran library (only for HDF5 output of `bin2hdf`)
* ghostscript to convert plots (`scripts/ps2png.sh`), python3 + numpy for the tests

## Usage

```sh
make                          # gfortran + OpenMP; for Intel: make COMPILER=intel
make test                     # regression tests (synthetic data, analytic solutions)

cp namelist/mag3d.nml .       # copy an example and edit file names, grid and contents
OMP_NUM_THREADS=16 ./bin/mag3d mag3d.nml         # -> .glb
```

* Each tool runs as `./bin/<tool> <namelist file>`. Example settings are in `namelist/` (with Japanese comments in `namelist/ja/`)
* Batch job example (Slurm): `scripts/job_example.sh`
* PostScript plots are converted to PNG / PDF by `scripts/ps2png.sh file.ps [pdf]`
* 3-D display: open `viewer/viewer.html` in a browser and choose or drag & drop .glb files
  (layer visibility, opacity, selectable-slice sliders, views, colour bars, PNG export)

Physical units and the planet (the `&units` group of each namelist file):

```
&units
  planet = 'jupiter',          ! 'jupiter' / 'saturn' / 'earth' / 'custom'
  temp_unit = 'eV'
/
```

See the [user manual](doc/manual_en.md) for details ([Japanese manual](doc/manual_ja.md)).

## Directories

```
src/        Fortran sources
namelist/   example settings of each tool (namelist/ja/: with Japanese comments)
viewer/     glTF viewer (three.js), X3D / VRML viewer (X_ITE)
textures/   planet images (Earth: Natural Earth; Jupiter, Saturn: schematic)
scripts/    PostScript conversion, batch job example
test/       regression tests
doc/        user manuals (manual_en.md, manual_ja.md)
```

## History and acknowledgements

The basic parts of these visualization tools were developed and published on the web by Prof. Tatsuki Ogino of Nagoya University.
FUKAZAWA, who studied under Prof. Ogino, has modernized them for present-day computers and publishes them together with
the improvements and additional tools he has developed since then (not published before). We are deeply grateful to
Prof. Ogino, who laid the foundation of these tools.

* Original visualization tools: Tatsuki Ogino
* Improvements, additional tools and modernization: Keiichiro FUKAZAWA (Research Institute for Humanity and Nature, ORCID [0000-0002-4862-6411](https://orcid.org/0000-0002-4862-6411))
* Earth image: Natural Earth (public domain, https://www.naturalearthdata.com/)
* The Jupiter and Saturn images are schematic band patterns (generated by `textures/make_textures.py`). To use real
  images, replace them with equirectangular images of the same names
* Viewers: three.js (MIT License) and X_ITE (MIT License), loaded from a CDN

## Citation

If you use this software in your research, please cite it using the information in [CITATION.cff](CITATION.cff) (Zenodo DOI: [10.5281/zenodo.24030429](https://doi.org/10.5281/zenodo.24030429)):

> Ogino, T., & FUKAZAWA, K. (2026). MHD viztools Fortran: visualization tools for planetary magnetosphere MHD simulations (Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24030429

## License

[MIT License](LICENSE)
