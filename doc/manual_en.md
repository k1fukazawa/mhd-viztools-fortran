# MHD viztools Fortran user manual

Fortran tools to combine, analyse and visualize the output of global MHD simulations of planetary magnetospheres.
They are a modernized version of the late Prof. Tatsuki Ogino's visualization tools. Below, "legacy code" / "legacy tools" refer to these.

* Fortran 2008 + OpenMP. External libraries are needed only for HDF5 output
* Grid sizes, paths, file numbers and plot contents are all given in namelists (no recompilation)
* Output file formats are identical to the legacy tools (`m.data`, `.vort`, `.h5`/`.xmf`); existing Python tools and ParaView settings keep working
* `make test` checks the tools
* **Physical units**: the visualization and analysis tools multiply normalized quantities by the normalization constants and plot / write physical units (default; set in the `&units` group, see below)

## Contents

```
viztools/
  README.md          overview (GitHub top page) / README.ja.md (Japanese)
  doc/manual_en.md   this manual (Japanese: doc/manual_ja.md)
  LICENSE            MIT License
  CITATION.cff       citation information
  Makefile
  src/
    mhd_common.f90   common: stream access to sequential files, file names, timing
    mhdunits.f90     common: normalized -> physical units (&units, constants of each planet)
    combine3d.f90    combination of per-rank files
    vorticity.f90    vorticity parallel / perpendicular to B
    bin2hdf.F90      HDF5/XDMF conversion (MHD and vorticity)
    psplot.f90       PostScript plotting library
    slice_plot.f90   slice plots
    polar_map.f90    mapping to the polar ionosphere (FAC, convection, open-closed boundary)
    mag3d.f90        3-D visualization (field lines, slices, isosurfaces, vectors -> glTF / X3D)
    mhdfield.f90     partial reading, interpolation and field-line tracing of 3-D data (for mag3d)
    scene3d.f90      3-D scene and writers for glTF (.glb) / X3D / PNG
    isosurf.f90      isosurface extraction (marching tetrahedra)
  textures/          planet images (earth.jpg: Natural Earth, jupiter/saturn.jpg: schematic) and generator script
  viewer/
    viewer.html      glTF (.glb) viewer (three.js, open in a browser)
    x3d_view.html    viewer for X3D and legacy VRML (.wrl) (X_ITE)
  namelist/          example settings of each tool (namelist/ja/: with Japanese comments)
  scripts/           ps2png.sh (PS -> PNG/PDF), job_example.sh (Slurm batch job example)
  test/              regression tests (synthetic data, analytic solutions)
```

## Build

```sh
make                                   # gfortran + OpenMP (bin2hdf with raw output only)
make COMPILER=intel                    # ifx; for ifort: make COMPILER=intel FC=ifort
make COMPILER=nec                      # nfort (SX-Aurora)
make COMPILER=fujitsu                  # frtpx
make OMP=no                            # without OpenMP
make HDF5=yes HDF5_DIR=/path/to/hdf5 COMPILER=intel   # enable HDF5 output
make test                              # regression tests (gfortran, python3 + numpy)
```

If the HDF5 include / library directories are not `$HDF5_DIR/include`, `$HDF5_DIR/lib`, give
`H5INC=-I...` and `H5LIB="-L... -lhdf5_fortran -lhdf5 -lz"`.
Executables are created in `bin/`. The `-mcmodel=large` of the legacy Makefile is not needed (no huge static arrays).

## Workflow

```
per-rank files ju03nt000.tdata.NNNN (e.g. 2048 files)
   | combine3d
   v
ju03nt000m.data  (3002x2002x1002 x 8 variables, sequential unformatted)
   |- slice_plot  -> slice plots (PostScript, convertible to PNG/PDF)
   |- polar_map   -> polar plots (FAC, convection, open-closed boundary) + raw binary maps
   |- mag3d       -> 3-D scenes (.glb / .x3d) -> viewer.html, Blender, ParaView
   |- vorticity   -> ju03nt000.vort (2 variables) -> slice_plot / bin2hdf
   '- bin2hdf     -> ju03nt000.h5 + .xmf -> ParaView / Python
```

Each tool is run as `./bin/<tool> <namelist file>`. The number of threads is set by `OMP_NUM_THREADS`.

## Physical units and precision (all tools)

The simulation works with normalized quantities; for visualization and analysis they are converted back to
physical units by multiplying the normalization constants. Write a `&units` group in the namelist file of each tool
(if absent, Jupiter values and physical units are used).

```
&units
  phys_units = .true.,        ! .false. keeps normalized values (same as the legacy tools)
  planet = 'jupiter',         ! 'jupiter' / 'saturn' / 'earth' / 'custom'
  ns = 1.0e10, bs = 4.2e-4, xs = 7.14e7, ion_mass = 1.0,   ! basic quantities (override only those given)
  ! rhos = ..., vs = ..., ts = ..., ps = ..., temps = ...  ! derived quantities may also be given
  length_unit = 'RJ',         ! distance unit of the plots (default: Jupiter RJ, Saturn RS, Earth RE)
  temp_unit = 'eV'            ! 'eV', 'keV', 'K'
/
```

**Normalization**: given the four basic quantities (number density ns [m⁻³], equatorial surface field bs [T],
planetary radius xs [m], ion mass ion_mass [amu]), the others are computed as

```
rhos = ns * ion_mass * 1.67e-27      vs = bs / sqrt(mu0 rhos)      ts = xs / vs
ps   = bs^2 / mu0 (= rhos vs^2)      temps [eV] = ps / (ns kB) / 11600
```

| planet | ns [m⁻³] | bs [T] | xs [m] | rhos [kg/m³] | vs [m/s] | ts [s] | ps [Pa] | temps [eV] | distance |
|---|---|---|---|---|---|---|---|---|---|
| jupiter (default) | 1.0e10 | 4.2e-4 | 7.14e7 | 1.67e-17 | 9.167e7 | 0.779 | 0.1403 | 8.76e7 | RJ |
| saturn | 1.0e10 | 2.0e-5 | 6.004e7 | 1.67e-17 | 4.359e6 | 13.77 | 3.185e-4 | 1.989e5 | RS |
| earth | 1.0e10 | 3.12e-5 | 6.378e6 | 1.67e-17 | 6.80e6 | 0.937 | 7.75e-4 | 4.839e5 | RE |
| custom | required | required | required | auto | auto | auto | auto | auto | `length_unit` (default R) |

The Jupiter, Saturn and Earth values are the normalization constants of the respective MHD codes; temps = ps/(ns kB)/11600.
Giving basic quantities in the namelist recomputes the derived ones (the values used and those computed are printed at run time).
The plot header (header of slice_plot and polar_map) is chosen from the planet name if left empty.

| quantity | conversion | display unit | factor (Jupiter) |
|---|---|---|---|
| rho (number density) | rho × ns | cm⁻³ | 1.0e4 |
| vx, vy, vz, absv | × vs | km/s | 9.167e4 |
| p | × ps | nPa | 1.403e8 |
| bx, by, bz, absb, absbz | × bs | nT | 4.2e5 |
| temp = p/rho | × temps (= ps/(ns kB)) | eV (keV, K) | 8.76e7 |
| pdyn = rho v² | × rhos vs² | nPa | 1.403e8 |
| beta | 2 μ0 ps p/(bs² B²) | – | |
| current density (polar_map) | × bs/(μ0 xs) | µA/m² | 4.68 |
| potential (polar_map) | × vs bs xs | kV | 2.75e9 |
| vorticity (output of vorticity) | × 1/ts | 1/s | 1.284 |
| energy flux proxy | × ps sqrt(ps/rhos) | mW/m² | |

* **Colour ranges (vmin/vmax, vecmax, iso_level, jmhd_max, ...) are given in display units.** Multiply values of old
  namelists by the factors above (e.g. B 2.5e-5 → 10.5 nT, T 1.84e-6 → 161 eV, v 6.88e-3 → 630.7 km/s, log10 rho +4)
* Simulation parameters given as input (`b0`, `omega`, `h`, ...) stay normalized
* Data files written by `combine3d` and `vorticity` (`m.data`, `.vort`) stay normalized (they are inputs of the other tools).
  The reading tools (slice_plot, bin2hdf) convert them. For vorticity files give `varunit = 'omega','omega'`
* `bin2hdf` converts to physical units by default and writes the HDF5 attributes `units`, `scale_from_normalized`
  and an XDMF `<Information Name="units">`. If your Python scripts expect normalized values, use `phys_units = .false.`
* Each tool prints the unit table and the value ranges (data range) of slices and maps at run time
* Pressure is normalized by ps = rhos vs² = bs²/μ0. A NOTE is printed if ps differs from this by more than 5%

**Precision**: the data files are single precision (about 7 digits), which limits the precision of the values.
Inside the tools, derived quantities such as p/rho, unit conversions, vorticity (differences and cross products) and
field-aligned currents (differences of B − B_dip) are computed in double precision; only results are stored in single
precision. Values in display units are within about 1e-6 to 1e6, so single precision is sufficient for storage.
ω⊥ = |ω| sqrt(1 − cos²) of the legacy vorticity lost digits where ω and B are nearly parallel; it is now computed as |ω × B|/|B|.

**Intel compilers**: the default of ifx/ifort (`-fp-model fast`) optimizes assuming that there are no NaNs, which can
break the tests for missing values (NaN). `COMPILER=intel` in the Makefile adds `-fp-model precise` (also when `FFLAGS` is given).

---

## combine3d (combination of per-rank files)

```sh
OMP_NUM_THREADS=32 ./bin/combine3d namelist/combine.nml
```

| variable | meaning |
|---|---|
| `npex,npey,npez` | number of ranks in each direction; rank = i + npex*(j + npey*k) |
| `nxx,nyy,nzz` | grid points per rank (no ghost cells) |
| `ox0,oy0,oz0` / `onx,ony,onz` | start (global grid index, 1-based) and size of the output region |
| `diri, prefix_in, suffix_in` | input `diri//prefix_in//NNN//suffix_in//RRRR` |
| `diro, prefix_out, suffix_out` | output `diro//prefix_out//NNN//suffix_out` |
| `ntap, nfiles, nstep, ndigit` | nfiles files numbered ntap, ntap+nstep, … |
| `overlap` | overlap writing with reading the next slab (default .true.) |

The cut-out by `mx=3000, nxp=1000, mxp=500` etc. of the legacy code corresponds to `ox0=1000, oy0=1000, oz0=1, onx=3002, ony=2002, onz=1002` (`namelist/combine.nml`).

**Method**

1. Only ranks overlapping the output region are read (450 of 2048 with the standard settings); edge ranks in y are read only for the rows needed
2. One rank layer in z (nzz planes) is assembled as a slab and written plane by plane
3. Rank files are read by stream access, jumping directly to the records of variable m. Record markers (4/8 bytes) and endianness are detected from the first record
4. Ranks in a slab are read in parallel with OpenMP tasks, overlapping the writing of the previous slab

| | legacy | new |
|---|---|---|
| data read (standard settings) | 1536 GB (all ranks, all variables) | about 240 GB |
| memory | 192 GB (whole 6000x4000x2000 array) | about 6 GB (two slabs) |
| parallelization | none | OpenMP (parallel rank reading + overlapped writing) |
| output | — | bit identical to the legacy code |

## vorticity

```sh
OMP_NUM_THREADS=16 ./bin/vorticity namelist/vort.nml
```

From ω = ∇×v (2nd-order central differences) computes ω_par = ω·B/|B| and ω_perp = |ω×B|/|B| and writes them to `.vort` (2 variables, same format as `m.data`).
With `out_wx, out_wy, out_wz, out_wmag` = .true. the components of ω and |ω| are also written (in the order par, perp, wx, wy, wz, wmag of those selected).

* Rolling buffer of 3 z planes. Memory: legacy about 265 GB → new about 0.35 GB
* Only the records of the 6 variables needed (v, B) are read (the legacy code read all 8 and copied them transposed)
* Computed in double precision (difference from the legacy code about 1e-6 relative). ω⊥ = |ω×B|/|B| does not lose digits where ω and B are nearly parallel

**Differences from the legacy code**

* Boundary planes (i, j, k = 1 and the last) are written as 0. The legacy code did not compute them and wrote the values of ρ and vx there
* At |B| = 0 or ω = 0 the legacy code produced NaN. The new code gives ω_par = 0, ω_perp = |ω|
* The output uses the same record markers and endianness as the input

## bin2hdf (HDF5 + XDMF for ParaView)

```sh
OMP_NUM_THREADS=2 ./bin/bin2hdf namelist/bin2hdf_mhd.nml    # 8 MHD variables
OMP_NUM_THREADS=2 ./bin/bin2hdf namelist/bin2hdf_vort.nml   # 2 vorticity variables
```

Merges the legacy HDF5 converters for MHD data and vorticity. The number and names of variables are given by `nvar`, `varname`.
HDF5 dataset names and shape (nz2 x ny2 x nx2 in C order) are the same as before.

| variable | meaning |
|---|---|
| `filetype` | `'native'` (default): native data type in the file. `'be'`: big endian as the legacy tools |
| `out_format` | `'hdf5'` (default) / `'raw'`: raw binary (.raw) and an XDMF file referring to it, without HDF5 |
| `nslab` | number of z planes written at once (default 32) |
| `compress` | gzip level (0: none). Compression is slower |
| `origin` | `'legacy'` (default): Origin values of the legacy tools / `'planet'`: planet center at the origin [planetary radii] |

* Legacy: a whole variable (24 GB) was read before writing, and all data were byte-swapped to big endian on writing. New: written as hyperslabs of `nslab` planes, without conversion. Memory 24 GB → 1.5 GB
* Reading the next slab overlaps the HDF5 write. HDF5 is always called by one thread, so no thread-safe HDF5 is needed; 2 threads are enough
* Values are identical to the legacy tools. h5py and ParaView handle byte order automatically
* `out_format='raw'` is for systems without HDF5. Readable with the ParaView XDMF reader (Xdmf2)

## slice_plot (slice plots)

```sh
OMP_NUM_THREADS=8 ./bin/slice_plot namelist/slice_legacy.nml   # same plot as the legacy slice-plot code
./scripts/ps2png.sh zxju03nt00m.ps           # PNG (one image per page)
./scripts/ps2png.sh zxju03nt00m.ps pdf       # PDF
```

One file (time step) per page. The namelist has two groups: `&slice` (files, grid, page) and `&panels` (each panel).

### Choosing slices

* `plane(p)` = `'xy'`, `'xz'`, `'yz'`
* `pos(p)` = slice position: distance from the planet center in simulation coordinates (x positive anti-sunward, z>0 north)
* Linear interpolation between grid points. `pos=0` is the same as the legacy "z=0 plane = mean of k=1,2" and "x=0 plane = mean of i=ixc, ixc+1"
* The grid position of the planet center is `i0 = nx2/2 - nxp + 0.5`, `j0 = ny2/2 + 0.5`, `k0 = 1.5` (default); give `i0, j0, k0` otherwise
* `mirror_z=.true.` (default) draws z<0 of xz and yz slices using the equatorial symmetry. Signs are `parity` (default ρ,vx,vy,p,Bz: +1, vz,Bx,By: −1). An xy slice with `pos<0` is made from z>0 with the same symmetry
* The plotted range can be zoomed by `h_lo, h_hi` (horizontal), `v_lo, v_hi` (vertical) [simulation coordinates]. The aspect ratio is automatic
* yz slices are mirrored horizontally by default (horizontal axis Y = −y as in the legacy plot). Change with `hflip, vflip` (0/1)

### Choosing quantities

Give a name in `qty(p)`.

| name | quantity |
|---|---|
| `rho vx vy vz p bx by bz` | variables in the file (names can be changed by `varname`) |
| `temp` | p/ρ (as the legacy T) |
| `absb`, `absbz`, `absv` | \|B\|, \|Bz\|, \|V\| |
| `beta` | 2μ0 p/B² (2p/B² in normalized units) |
| `pdyn` | ρV² |
| `log10:<qty>`, `abs:<qty>`, `neg:<qty>` | log10, absolute value, sign change (e.g. `log10:rho`) |

* `vmin, vmax`: colour range. Equal values (default 0, 0): automatic (minimum to maximum)
* `palette`: 1–8 (same numbers as the legacy code). 3 = colours of the legacy plot, 8 = blue-white-red (for signed quantities)
* `nsmooth`: number of 5-point smoothing passes as in the legacy code (default 2)
* `stride`: image thinning. 2 makes large PS files about 1/4
* `vec`: `'v'` (velocity) / `'b'` (magnetic field) / `'none'`. In-plane components drawn as streamline arrows as in the legacy code. Normalized by `vecmax` (0: maximum), seed points `nvec_h, nvec_v`
* `row, col`: position on the page (row=1 at the top)
* For `.vort` files use `nvar=2, varname='par','perp', parity=1,1` (`namelist/slice_vort.nml`)

### Files and titles

* `nfiles` is the number of files (1 or more). `ntap=0, nfiles=1` processes `ju03nt000m.data` only. 0 stops with an error
* A relative `psfile` is written in the run directory (the submit directory for batch jobs). The absolute path is printed on the `PostScript:` line of the log
* List files in `datafile(i)` or use the pattern `diri, prefix, ntap, nfiles, nstep, suffix`
* Titles from `title(i)`, or automatically `title_head` + `t0 + (i−1)*dt` + `t_unit`
* `ps_per_file=.true.` writes a separate PS file per time step (`<psfile>_001.ps`, …)
* `dump_planes=.true.` writes the 2-D array of each panel (before smoothing) as `.bin` (float32) with a description `.txt`, e.g. for re-plotting in Python

### Data read and memory

The legacy code read the whole 3-D data (8 variables, 192 GB) for each time step into a 24 GB array. The new code reads only the parts needed for the slices.

| slice | data read (per variable) |
|---|---|
| xy | 2 planes (48 MB) |
| xz | 2 rows of each z plane (24 MB) |
| yz | `yz_read='full'`: all z planes (24 GB) / `'strided'`: 2 words per row |

Memory is a few slices (about 0.5 GB on the standard grid). yz slices take most of the reading time; on Lustre `yz_read='strided'` may be faster, so try both on your system.

### Differences from the legacy plot

* Ticks (every `tick`), axis names, colour bars with units and the slice position are drawn. With `label_coord='jso'` (default) the ticks show X = −x, Y = −y, Z = z (same orientation as the legacy labels); `'sim'` gives simulation coordinates
* Magnetic field vectors on xy slices are (Bx, By); the legacy code drew (Bx, −By) (nearly the same at z≈0 where By≈0 by symmetry)
* The coefficient error ((1−Z)·R → (1−R)·Z) of the uy interpolation in the arrow drawing of the legacy code is fixed
* The legacy smoothing behaviour of not smoothing the middle row is removed
* Images use 4 bits per colour as before; hex strings are built in a character buffer, so writing is fast

## polar_map (mapping to the polar ionosphere)

```sh
OMP_NUM_THREADS=32 ./bin/polar_map namelist/polar.nml
```

Merges the three legacy polar-plot programs (open-closed boundary, convection and energy flux, plotting) into one. The combined `m.data` is read directly, so neither `sdata` nor intermediate files (`.dp3`, `.dp5`) are needed.

### Method

A grid of "angle from the pole × MLT" (default 0.25° × 15 min) is taken on the ionosphere (R = 1 planetary radius), and at each point:

1. **Mapping point**: mapped from the foot point along the dipole field line to radius r_m = min(`r_map`, L) (L = 1/sin²θ). `r_map` is just outside the inner boundary (default 10 RJ for an inner boundary at 7 RJ). Inside it the field is regarded as a dipole
2. **Field-aligned current (MHD)**: J = ∇×(B − B_dip) at the mapping point, mapped to the ionosphere using the conservation of J∥/B along the field line:
   `j_up = (J·B / B²)_map × B_r,ion` (upward positive). The dipole is subtracted before differencing, so discretization errors of the dipole do not enter the current
3. **Convection**: the velocity at the mapping point (rigid corotation Ω ẑ×r subtracted by default) is converted to the ionospheric velocity with the kinematics of the dipole mapping:
   `u_θ = −(dL/dt) / (2 L^1.5 cosθ0)`, `u_φ = sinθ0 · v_φ / (r sinθ)`, `dL/dt = v_r/sin²θ − 2 cosθ v_θ / sin³θ`
4. **Potential and FAC (from convection)**: ∇·E from E = −u×B_ion; the Poisson equation on the sphere ∇²Φ = −∇·E is solved
   directly by DFT in azimuth + tridiagonal solves in colatitude (no iteration, Φ=0 at the boundary).
   `j_conv = −∇·(Σ_P E)` (uniform Pedersen conductance Σ_P). The "Field-Aligned Current" of the legacy plot corresponds to this quantity
5. **Open-closed boundary**: traced from the mapping point along B (2nd-order Runge–Kutta). For half-domain data a line reaching z=0 is closed (it returns to the conjugate point by symmetry); a line leaving the domain is open.
   At the same time ⟨p·sqrt(p/ρ)⟩ (the "Energy Flux" proxy) is averaged up to the equator with the weight r²/(r_in²+r²) of the legacy code

Near the mapping points (radius `r_map` + a few cells) the data are read at full resolution; the global field for tracing (B, ρ, p) is thinned by `stride` (standard grid, stride=3: about 4.5 GB in memory, about 40 GB read). Tracing is OpenMP parallel.

### Output

* PostScript (A4 portrait, one page per time step): (1) MHD field-aligned current (2) potential and convection vectors (3) FAC from convection (4) energy flux proxy. The open-closed boundary (thick black line) is drawn on all panels
* `<psfile>_NNN.bin` / `.txt`: all maps (float32, (ncolat, nmlt, 12)). The items are listed in the `.txt` (j_mhd, j_conv, phi, u_colat, u_mlt, topology, eflux, length, r_eq, r_map, two components of the potential flow).
  Values in the physical units of `&units` (j: µA/m², phi: kV, u: km/s (speed at R = 1 planetary radius), eflux: mW/m²). The unit of each item is also written in the `.txt`

### Main settings

| variable | meaning |
|---|---|
| `hemisphere`, `view` | `'south'` (default, seen from below = 06 left, 18 right as the legacy plot) / `'north'` (seen from above) |
| `half_domain` | `.true.`: m.data (z≥0). `.false.` for full-domain data such as the legacy sdata |
| `r_map`, `r_valid` | mapping radius and minimum radius of the mapping points used |
| `subtract_corot`, `omega` | whether corotation is subtracted; omega is vo0 of constants.f90 |
| `colat_max`, `ncolat`, `nmlt` | ionospheric grid |
| `stride`, `box_lo`, `box_hi` | thinning and region of the global field for tracing |
| `jmhd_max`, `jconv_max`, `eflux_max`, `umax` | colour / arrow ranges (0: automatic), in the units of `&units` (µA/m², mW/m², km/s) |
| `sigma_p` | uniform Pedersen conductance [S] for j_conv (= −Σ_P ∇·E) (physical units; default 1 S) |
| `jfac` | extra factor for the currents (default 1; units are converted by `&units`) |
| `jpar_delta` | difference width of J = ∇×B (multiple of the grid spacing h, default 2). With 2, numerical oscillations of wavelength 4h are not amplified by rot |
| `jpar_avg` | length along the field line over which (J·B)/B² is averaged (default 4 RJ; 0: mapping point only). J∥/B is conserved along B, so averaging does not change the physical quantity |
| `nsmooth` | smoothing passes of the plots (5-point as in the legacy code, default 2). Not applied to the .bin data |

### Implementation note: row boundaries of PostScript images

PostScript image rows must start at byte boundaries; with 4 bits per colour and an odd width (polar_map uses 401 pixels) each row needs one padding digit. Without it each row shifts by half a byte, giving a diagonal saw-tooth pattern, so ps_image in psplot pads the rows. polar_map draws with 8 bits per colour (no colour steps).

### Implementation note: small-scale oscillations of the field-aligned current

If J at the mapping point is computed by central differences of width h, oscillations of about 4 cells in the MHD solution are amplified by rot and appear on the ionospheric grid as stripes or dots.
Therefore the difference width is 2h (`jpar_delta`), the value is averaged over 4 RJ along the field line (`jpar_avg`), and the plots get two passes of the legacy 5-point smoothing (`nsmooth`).
In a test with synthetic data plus an oscillation of wavelength 4h, the noise/signal ratio of j_mhd dropped from 0.28 to 0.002 (up to 17° colatitude).

### Problems found in the legacy code (polar plots)

| part | problem | effect |
|---|---|---|
| open-closed boundary | `ia(1)=nxg=10` is passed as the number of foot points, so only 12×12 instead of 202×202 points are computed; results fill only a corner of the 202×202 array (rows 191–202, columns 1–12) | the open-closed boundary is hardly drawn. Reproduced with synthetic data |
| open-closed boundary | the inner dipole is built with `b0=+1`, while the simulation uses `b0=−1`; the two cancel at r ≈ 4.3 RJ and tracing turns back there | open/closed classification fails |
| open-closed boundary | `ardr` is not initialized | results depend on the system |
| convection | the x, y components of the velocity at the mapping point (r=10) are used directly as horizontal ionospheric components | the z component is discarded and the meridional mapping is wrong; for mapping points at colatitudes above 54.7° (foot points below about 15°) outward motion maps poleward but is drawn equatorward (wrong sign) |
| convection | the electric field is mapped only by an isotropic ratio sqrt(B) | the dipole mapping scales differently (by about 2) in the meridional and azimuthal directions, so the vorticity (= FAC) is distorted |
| plotting | the velocity is multiplied by artificial windows `1/(1+R^8)` and `(1−r²)^0.25` before taking rot | the gradients of the windows create a spurious FAC ring at the outer edge (around 70°) |
| plotting | the Poisson equation is solved by ADI with fixed parameters, stopped after 100 iterations, over the whole square (including corners outside the polar cap) | the largest structure (lowest mode) reaches only 8% of the true solution (checked with the legacy routine); hence the contours in the corners and near the edge |
| plotting | the divergence-free velocity is renormalized by `max|u|/max|v_df|` | amplitude errors are hidden in the plots |
| plotting | FAC is computed as "vorticity of the flow × uniform conductance" | not the MHD current itself. The new code mainly plots the mapped MHD J∥ and keeps the convection-based value for comparison |

### Verification (`test/polar_test.py`, part of `make test`)

Compared with analytic solutions for synthetic data of a dipole + uniform field + localized twisted current + (corotation + flow with dL/dt = const).

| item | result |
|---|---|
| j_mhd (from the analytic rot B) | relative error 2.7×10⁻³ |
| mapped convection u_colat (analytic −c/(2L√(L−1))) | relative error 2.4×10⁻³ |
| open-closed boundary (13.95° from tracing the same field in python) | within one cell (0.25°) |
| mapped corotation u_mlt = Ω sinθ, potential Ω(sin²θ − sin²θmax) | relative errors 5×10⁻⁷, 3×10⁻⁵ |
| northern and southern hemispheres (mirrored) | identical |

## mag3d (3-D visualization: successor of the legacy VRML tools)

```sh
OMP_NUM_THREADS=16 ./bin/mag3d namelist/mag3d.nml      # -> ju03nt3d.glb, ...
```

Combines Prof. Ogino's VRML subroutines and the legacy programs using them (field lines, stacked slices, isosurfaces, velocity vectors) into one program. It reads the combined `m.data` directly and writes one 3-D scene per time step.

### Output formats

| format | advantages | disadvantages | in this tool |
|---|---|---|---|
| VRML 2.0 (.wrl) | legacy assets; easy to generate as text | standard and viewers no longer developed; ASCII, large and slow | not written. Legacy files can be shown with `x3d_view.html` (X_ITE) |
| **X3D (.x3d)** | ISO successor of VRML with nearly the same node structure; text, viewpoints and animation within the standard | XML, large (similar to VRML); few applications | written with `out_x3d=.true.`; images as external PNG |
| **glTF 2.0 (.glb)** | binary, small and fast (about half of X3D); supported by three.js, Blender, ParaView, the Windows 3D viewer; passed directly to WebGL | text (labels) not in the standard (stored in extras and shown by the viewer) | **default**; shown with `viewer/viewer.html` |
| ParaView (HDF5/XDMF) | full analysis such as volume rendering and interactive streamlines | software needed to share results | via the existing `bin2hdf` |

For large data the heavy work (field-line tracing, isosurfaces) is done by Fortran on the supercomputer, and only the light geometry (.glb) is viewed in a local browser. The 3-D data themselves need not be transferred.

### Contents

| element | namelist |
|---|---|
| field lines (from a polar foot-point grid, latitude × MLT) | `fl_polar`, `fl_hemi`, `fl_lat_min/max`, `fl_nlat`, `fl_nmlt`, `fl_r_start` |
| field lines (from an equatorial grid, both directions) | `fl_equator`, `eq_xmin` … `eq_ny` |
| field lines (from any points) | `fl_npt`, `fl_px/py/pz` |
| field lines (from a 3-D grid; for lines not connected to the planet) | `fl_box`, `fl_box_lo/hi`, `fl_box_n`, `fl_keep` |
| only lines near the open-closed boundary | `fl_select='ocb'` |
| field-line colours | `fl_color='topology'` (closed: green / open: blue / solar wind: magenta / detached: grey) or a quantity name |
| slices (colour images, translucent) | `pl_plane`, `pl_pos`, `pl_qty`, `pl_vmin/vmax`, `pl_alpha`, `pl_alpha_mode` |
| stacked slices (pseudo volume, up to 3 directions at once) | `stk_plane(1:3)`, `stk_n`, `stk_lo/hi`, `stk_qty`, `stk_vmin/vmax`, `stk_tra`, `stk_alpha` (maximum opacity, default 0.4) (e.g. `stk_plane = 'yz','xz','xy'`). The viewer's "Stacks" slider changes the opacity of all sets at once (the "Slices" slider also applies) |
| selectable slices (one shown at a time in the viewer, up to 3 sets) | `sel_plane`, `sel_n`, `sel_lo/hi`, `sel_qty`, `sel_vmin/vmax`, `sel_alpha`, `sel_palette` |
| isosurfaces | `iso_qty`, `iso_level`, `iso_color`, `iso_alpha`, `iso_stride`, `iso_max_tri` (default 2 million; thinned automatically above it), `iso_smooth` (smoothing passes before extraction, default 1) |
| vector arrows (velocity, magnetic field) | `arr_plane`, `arr_pos`, `arr_vec`, `arr_n` |
| frame, axes, ticks, title, planet | `draw_frame`, `draw_planet`, `draw_inner`, `label_size` |
| planet image, Saturn's rings | `planet_texture` ('auto' / file name / 'none'), `planet_lon0`, `draw_rings`, `ring_r1/r2` |
| colour bars | shown by the viewer; `colorbar_3d=.true.` also draws them in the scene |

Quantity names are the same as in slice_plot (rho, vx … bz, temp, absb, absbz, absv, beta, pdyn, log10:<qty> etc.). For half-domain data (z≥0), field lines, slices and isosurfaces are drawn down to z<0 using parity.
Give `pl_qty(n)` for every index given in `pl_plane(n)` (otherwise an error; this catches index mistakes such as `pl_plane(2)= ..., pl_qty(1)= ...`).
Values are in the physical units of `&units` (also `pl_vmin/vmax`, `stk_vmin/vmax`, `iso_level`, `fl_vmin/vmax`, `arr_vmax`).
A label "quantity [unit] colour range" is drawn above each slice (`plane_labels=.false.` removes it).
The log shows the colour range and the data range of each slice; a WARNING is printed if there is no data, the data are constant or outside the colour range.
If a slice appears completely black or transparent, look here first.

### Field lines not connected to the planet

Field lines traced in both directions from the equator, from points or from a 3-D grid (`fl_box`) are classified by where their two ends go.

| layer name | class | colour (`fl_color='topology'`) |
|---|---|---|
| `fieldlines_closed` | both ends on the planet | green |
| `fieldlines_open` | one end on the planet, one leaving the domain | blue |
| `fieldlines_solarwind` | both ends leave the domain (solar wind, IMF) | magenta |
| `fieldlines_detached` | neither reaching the planet nor leaving the domain (closed loops, lines reaching the maximum length) | grey |

`fl_keep` selects the lines kept (`'unconnected'` = solarwind + detached). Lines from the polar foot points (`fl_polar`) are always drawn.
A closed loop returning to its starting point is stopped after one turn. Example:

```
fl_box = .true., fl_box_lo = -60.0, -60.0, -30.0, fl_box_hi = 150.0, 60.0, 30.0, fl_box_n = 8, 6, 3,
fl_keep = 'unconnected',
```

Positions are in simulation coordinates, so negative x is sunward. To look only at the solar wind, restrict `fl_box_lo(1)` / `fl_box_hi(1)` to the sunward side.

### Selectable slices

As for stacked slices, `sel_n` parallel slices are made in advance, and the viewer's "Selectable slices" slider chooses the one shown
("Show all" shows all). Without a colour range the minimum to maximum of the global field is used for all planes. The slice positions are stored in the node extras
(`group`, `index`, `pos`) of the glb.

```
sel_plane(1) = 'yz', sel_n(1) = 21, sel_lo(1) = -50.0, sel_hi(1) = 250.0, sel_qty(1) = 'temp',
sel_vmin(1) = 0.0, sel_vmax(1) = 300.0,
```

### Colour bars

The colour ranges of slices, stacked slices, selectable slices, field lines coloured by a quantity and arrows are stored in `scenes[0].extras.legends` of the glb,
and viewer.html shows them at the lower right (only for visible layers; the "Colour bars" checkbox hides them). They are included in "Save image".
For Blender or X3D, `colorbar_3d = .true.` also draws colour bars in the scene (below the frame).

### Planet images

`planet_texture = 'auto'` (default) maps `textures/<planet>.jpg` (or .png) onto the sphere (north pole at the top of the image,
the image center = longitude 0 facing the Sun; rotate with `planet_lon0` [deg]). For Saturn the rings (C to A ring, 1.24–2.27 RS) are also drawn.

| file | content |
|---|---|
| `textures/earth.jpg` | Natural Earth I (public domain, naturalearthdata.com) resized to 1024×512 |
| `textures/jupiter.jpg` | **schematic** image of the bands and the Great Red Spot (generated by `textures/make_textures.py`; not a real image) |
| `textures/saturn.jpg` | **schematic** image of the bands (same) |

To use real images, replace them with equirectangular (width:height = 2:1) JPEG/PNG images of the same names, or give
`planet_texture = '/path/to/map.jpg'`. Public-domain images include the Cassini global map of Jupiter (PIA07782) from
NASA/JPL Photojournal. Images are embedded in the glb (for X3D `<outbase>_img1.jpg` is written).

### Method

* The global field is read with `stride` thinning and the region within `r_fine` of the planet at full resolution, only the variables needed (row by row, so about 1/stride² of the data is read)
* Field lines are traced with 2nd-order Runge–Kutta and thinned according to curvature and length (as the legacy VRML tools). From the foot point to `fl_r_start` dipole field lines are drawn
* Isosurfaces by marching tetrahedra. Vertices are shared per grid edge and normals come from the gradient of the scalar, giving smooth shading and small files
* Field-line tracing, slice sampling and isosurface extraction are OpenMP parallel
* Display coordinates (X sunward, Z north) are rotated to the Y-up coordinates of glTF/X3D. The legacy VRML wrote (x, z, y), which is a mirror image

### Viewing

* **glb**: open `viewer/viewer.html` in a browser and choose .glb files with "Open files" (several files = time animation), or drag & drop.
  Layer visibility, opacity of slices / stacks / isosurfaces, selectable-slice sliders, views (from the Sun, tail, dusk, north, oblique), white background, colour bars and PNG export.
  Via an http server, `viewer.html?files=a_001.glb,a_002.glb` opens files directly (`python3 -m http.server`).
  Blender (File → Import → glTF) gives high-quality images and movies. ParaView 5.10 or later can also open them
* **x3d / legacy wrl**: `viewer/x3d_view.html?file=demo.x3d` (via an http server, with the PNG files in the same place). Legacy VRML output (.wrl) can also be shown on this page (the Inline + TimeSensor mechanism used for time series also works)
* Both viewers load three.js / X_ITE from a CDN, so an internet connection is needed (for offline use see the comments at the top of each file)

### Performance (synthetic data, 2-core machine)

| data | settings | time | memory | output |
|---|---|---|---|---|
| 722×482×242 × 8 variables (2.7 GB) | stride 2, 1920 field lines, 2 slices, 1 isosurface (380k triangles), arrows | 1 thread 9.7 s / 2 threads 6.1 s | 0.44 GB | .glb 13 MB, .x3d 23 MB |

A legacy VRML output with about 1000 field lines and 3 slices was 4.4 MB of ASCII. On the standard grid (3002×2002×1002) the global field with stride 4 takes 0.38 GB per variable.

### Verification

`make test` writes .glb and .x3d from synthetic data and checks the glTF structure with the official validator (gltf-validator, if node is available; 0 errors and warnings). Browser display was checked with Chromium (three.js 0.170 / X_ITE 16.5).

---

## File formats

* `m.data` / `.vort`: Fortran sequential unformatted, `do m=1,nvar; do k=1,nz2; write(u) f(1:nx2,1:ny2,k)`. One record = one z plane (nx2*ny2 float32 values)
* The tools read them by stream access and jump directly to record (m, k), at `(irec-1)*(4*nx2*ny2 + 2*marker) + marker + 1` (irec = (m−1)*nz2 + k)
* The record-marker size and endianness are detected from the leading and trailing markers (files written with gfortran `-frecord-marker=8` or `-fconvert=big-endian` can also be read)

## Tests

```sh
make test
```

Each tool is run on small synthetic data and compared with analytic solutions or independent calculations (python3 + numpy required).
combine3d (normal / big endian / 8-byte markers × 1, 2, 4 threads) is compared with the generating formula, vorticity with a numpy calculation,
and polar_map with analytic solutions. bin2hdf (raw, HDF5) must be bit identical to the input.
Conversion to physical units is checked for bin2hdf (raw/XDMF/HDF5 attributes) and slice_plot (ratio to the normalized version).
The slice_plot output is written to `test/work/s/new.ps` (and `new.png` if ghostscript is available) for visual checks.

## Performance (small tests)

Measured on a 2-core machine, dropping the page cache each time. Values on production systems (Lustre, many cores) need to be measured separately.

| tool | test size | legacy | new (1 thread) | new (2 threads) |
|---|---|---|---|---|
| combine3d | 256 ranks, input 4.1 GB → output 545 MB | 6.6–8.1 s, 549 MB | 0.8–1.7 s, 69 MB | 1.1–1.4 s |
| vorticity | 512x256x130 | 1.24 s, 718 MB | 0.71 s, 10 MB | 0.77 s |
| bin2hdf | 512x256x130 x 8 variables | 1.0–1.2 s | 0.64–0.99 s | 0.79–1.07 s |
| slice_plot | 302x202x102, 4 panels | — (reads all data) | 0.2 s | — |

The main gains of combine3d are less data read (1/6 with the standard settings) and less memory, so the difference should be larger on production systems.
OpenMP helps when concurrent streams pay off on a parallel file system. Measure with `OMP_NUM_THREADS` = 8, 16, 32 … and use the fastest.

## Maintenance notes

* **Changing grids or paths**: edit only the namelists
* **Adding a quantity**: add a case to `mark_needs` (variables needed) and `eval_qty` (formula) in `slice_plot.f90`, and to `qty_parse` / `qty_needs` / `qty_eval` in `mhdfield.f90` for mag3d
* **Adding a colour map**: add a case to `palette_rgb` in `psplot.f90`
* **Parallelization**: do not use functions returning deferred-length strings or internal writes with run-time formats inside parallel regions (races in gfortran). `make_fname` and `zpad` return fixed-length strings; use them with `trim()`
* **Older Intel compilers**: avoid updating host-associated arrays from internal procedures (ifx 2023 failed to do this correctly); the variables needed by mag3d are therefore collected directly in the main program
* **Record markers of NEC/Fujitsu**: detected automatically. If detection fails the tools stop with "record length mismatch" (a wrong grid size gives the same error)
