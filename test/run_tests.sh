#!/bin/bash
#==============================================================================
# Regression tests: run each tool on small synthetic data and compare with analytic solutions / independent calculations
#
#   usage   :  make test   (or: make && cd test && ./run_tests.sh)
#   requires: gfortran (to generate test data), python3 + numpy
#   optional: HDF5 output is checked if h5py is available
#==============================================================================
set -u
BIN=$(cd ../bin && pwd)
FC=${FC:-gfortran}
W=$(pwd)/work
rm -rf "$W" && mkdir -p "$W"
pass=0; fail=0
ok()  { echo "  [PASS] $1"; pass=$((pass+1)); }
ng()  { echo "  [FAIL] $1"; fail=$((fail+1)); }

#------------------------------------------------------------------ combine3d
echo "== combine3d"
mkdir -p "$W/c" && cd "$W/c" && mkdir -p in outnew
cat > combine.nml <<'EOF'
&combine
 npex=2, npey=3, npez=4, nxx=10, nyy=8, nzz=6, nb=8,
 ox0=4, oy0=7, oz0=1, onx=10, ony=12, onz=12,
 diri='in/', diro='outnew/', prefix_in='ju03nt', ntap=0, nfiles=1
/
EOF
for conv in "" "-fconvert=big-endian" "-frecord-marker=8"; do
  rm -f in/*
  $FC $conv ../../gen_ranks.f90 -o gen && ./gen
  for t in 1 2 4; do
    OMP_NUM_THREADS=$t "$BIN/combine3d" > log_$t.txt 2>&1
    if python3 - <<'EOF'
import numpy as np
# gen_ranks.f90: f = m*1e6 + i + 100 j + 1e4 k + 0.25 (global grid indices); output covers (4..13, 7..18, 1..12)
nx, ny, nz = 10, 12, 12
raw = np.fromfile('outnew/ju03nt000m.data', np.uint8)
rec = 8 + 4*nx*ny
assert raw.size == 8*nz*rec, raw.size
i = np.arange(nx) + 4; j = np.arange(ny) + 7
for m in range(1, 9):
    for k in range(1, nz + 1):
        r = (m-1)*nz + k - 1
        a = np.frombuffer(raw[r*rec+4:r*rec+4+4*nx*ny].tobytes(), '<f4').reshape(ny, nx)
        ref = m*1e6 + i[None, :] + 100.0*j[:, None] + 1e4*k + 0.25
        assert np.array_equal(a, ref.astype(np.float32)), (m, k)
EOF
    then ok "input[$conv] threads=$t : matches the generated values"
    else ng "input[$conv] threads=$t"; fi
  done
done

#------------------------------------------------------------------ vorticity
echo "== vorticity"
mkdir -p "$W/v" && cd "$W/v" && mkdir -p new
$FC ../../gen_m.f90 -o genm && ./genm
printf "&vort\n nx2=40, ny2=30, nz2=20, diri='./', diro='new/', ntap=186\n/\n" > vort.nml
"$BIN/vorticity" > log.txt 2>&1
if python3 - <<'EOF'
import numpy as np
nx, ny, nz, h = 40, 30, 20, 0.15
def load(fn, nv):
    raw = np.fromfile(fn, np.uint8); rec = 8 + 4*nx*ny
    return np.stack([np.frombuffer(raw[r*rec+4:r*rec+4+4*nx*ny].tobytes(), '<f4').reshape(ny, nx)
                     for r in range(nv*nz)]).reshape(nv, nz, ny, nx).astype(np.float64)
f = load('ju02nt186m.data', 8); w = load('new/ju02nt186.vort', 2)
vx, vy, vz, bx, by, bz = f[1], f[2], f[3], f[5], f[6], f[7]
d = lambda a, ax: (np.roll(a, -1, ax) - np.roll(a, 1, ax))/(2*h)      # axes: 0=z, 1=y, 2=x
w1 = d(vz, 1) - d(vy, 0); w2 = d(vx, 0) - d(vz, 2); w3 = d(vy, 2) - d(vx, 1)
bm = np.sqrt(bx**2 + by**2 + bz**2); ok = bm > 0
par = np.where(ok, (w1*bx + w2*by + w3*bz)/np.where(ok, bm, 1), 0)
cx, cy, cz = w2*bz - w3*by, w3*bx - w1*bz, w1*by - w2*bx
perp = np.where(ok, np.sqrt(cx**2 + cy**2 + cz**2)/np.where(ok, bm, 1), np.sqrt(w1**2 + w2**2 + w3**2))
s = (slice(1, -1),)*3
for name, ref, new in [('par', par, w[0]), ('perp', perp, w[1])]:
    e = np.abs(new[s] - ref[s]).max()/np.abs(ref[s]).max()
    print('  %s: max relative difference %.1e' % (name, e)); assert e < 1e-5
EOF
then ok "vorticity agrees with an independent numpy calculation"; else ng "vorticity (numpy)"; fi

#------------------------------------------------------------------ bin2hdf (raw)
echo "== bin2hdf (raw format)"
mkdir -p "$W/h" && cd "$W/h" && mkdir -p new
cp ../v/ju02nt186m.data .
printf "&bin2hdf\n nx2=40, ny2=30, nz2=20, nvar=8, nxp=5, diri='./', diro='new/', ntap=186, nslab=7, out_format='raw'\n/\n&units phys_units=.false. /\n" > b.nml
"$BIN/bin2hdf" b.nml > log.txt 2>&1
$FC ../../cmp_raw.f90 -o cmp_raw
if ./cmp_raw ju02nt186m.data new/ju02nt186.raw 40 30 20 8; then ok "raw output identical to input"
else ng "raw output"; fi
if grep -q 'Format="Binary"' new/ju02nt186.xmf; then ok "XDMF (binary) written"; else ng "XDMF"; fi
# conversion to physical units (default): RHO x ns*1e-6 [cm^-3], VX x vs*1e-3 [km/s], P x ps*1e9 [nPa], BX x bs*1e9 [nT]
mkdir -p phys
printf "&bin2hdf\n nx2=40, ny2=30, nz2=20, nvar=8, nxp=5, diri='./', diro='phys/', ntap=186, nslab=7, out_format='raw'\n/\n" > bp.nml
"$BIN/bin2hdf" bp.nml > logp.txt 2>&1
if python3 - <<'EOF'
import numpy as np
n = np.fromfile('new/ju02nt186.raw', '<f4').reshape(8, -1).astype(np.float64)
p = np.fromfile('phys/ju02nt186.raw', '<f4').reshape(8, -1).astype(np.float64)
fac = [1e10*1e-6, 9.167e7*1e-3, 9.167e7*1e-3, 9.167e7*1e-3, 0.1403*1e9, 4.2e-4*1e9, 4.2e-4*1e9, 4.2e-4*1e9]
err = max(np.abs(p[m] - n[m]*fac[m]).max()/max(np.abs(n[m]*fac[m]).max(), 1e-30) for m in range(8))
print('  max relative error of the unit conversion: %.1e' % err)
x = open('phys/ju02nt186.xmf').read()
assert err < 1e-6 and 'Value="cm^-3"' in x and 'Value="nT"' in x
EOF
then ok "physical units in raw/XDMF (bin2hdf)"; else ng "physical units (bin2hdf)"; fi

# HDF5 (if bin2hdf was built with HDF5 and python3+h5py exists)
printf "&bin2hdf\n nx2=40, ny2=30, nz2=20, nvar=8, nxp=5, diri='./', diro='new/', ntap=186, nslab=7\n/\n&units phys_units=.false. /\n" > b5.nml
printf "&bin2hdf\n nx2=40, ny2=30, nz2=20, nvar=8, nxp=5, diri='./', diro='phys/', ntap=186, nslab=7\n/\n" > b5p.nml
if "$BIN/bin2hdf" b5.nml > log5.txt 2>&1 && python3 -c "import h5py" 2>/dev/null; then
  if python3 - <<'EOF'
import h5py, numpy as np
f = h5py.File('new/ju02nt186.h5', 'r')
raw = np.fromfile('new/ju02nt186.raw', '<f4').reshape(8, 20, 30, 40)
names = ['RHO','VX','VY','VZ','P','BX','BY','BZ']
assert all(np.array_equal(f[n][...], raw[i]) for i, n in enumerate(names))
EOF
  then ok "HDF5 datasets identical"; else ng "HDF5 datasets"; fi
  "$BIN/bin2hdf" b5p.nml > log5p.txt 2>&1
  if python3 - <<'EOF'
import h5py
f = h5py.File('phys/ju02nt186.h5', 'r')
u = f['BX'].attrs['units']; u = u.decode() if isinstance(u, bytes) else u
assert u == 'nT' and abs(f['BX'].attrs['scale_from_normalized'] - 4.2e5) < 1e-6
EOF
  then ok "HDF5 unit attributes"; else ng "HDF5 unit attributes"; fi
else
  echo "  [SKIP] HDF5 (bin2hdf built without HDF5, or no h5py)"
fi

#------------------------------------------------------------------ slice_plot
echo "== slice_plot"
mkdir -p "$W/s" && cd "$W/s"
$FC -O2 ../../gen_mag.f90 -o gen_mag && ./gen_mag
cat > slice.nml <<'EOF'
&slice
  ntap = 171, nfiles = 1, title(1) = 'test data',
  psfile = 'new.ps', nx2 = 302, ny2 = 202, nz2 = 102, nxp = 50, h = 1.5,
  dump_planes = .true.
/
&panels
  npanel = 4,
  plane = 'xy','yz','xy','yz', qty = 'absbz','absb','temp','temp',
  vmax = 2.5e-5, 2.5e-5, 1.84e-6, 1.84e-6,
  vec = 'b','b','v','v', vecmax = 2.5e-5, 2.5e-5, 6.88e-3, 6.88e-3,
  nvec_h = 40, 20, 40, 20, row = 1,1,2,2, col = 1,2,1,2, label = 'B','B','T','T',
/
&units phys_units = .false. /
EOF
if "$BIN/slice_plot" > log.txt 2>&1 && [ -s new.ps ]; then ok "new.ps written"; else ng "slice_plot"; fi
# physical units: the same plot in nT, eV, km/s, compared with the dumped planes of the normalized version
sed -e "s/psfile = 'new.ps'/psfile = 'phys.ps'/" -e "s/vmax = 2.5e-5, 2.5e-5, 1.84e-6, 1.84e-6/vmax = 10.5, 10.5, 161.2, 161.2/" \
    -e "s/vecmax = 2.5e-5, 2.5e-5, 6.88e-3, 6.88e-3/vecmax = 10.5, 10.5, 630.7, 630.7/" -e "/&units/d" slice.nml > slice_phys.nml
if "$BIN/slice_plot" slice_phys.nml > logp.txt 2>&1 && python3 - <<'EOF'
import numpy as np
fac = [4.2e5, 4.2e5, 8.76e7, 8.76e7]
for p in range(4):
    a = np.fromfile('new_f001_p%02d.bin' % (p+1), '<f4').astype(np.float64)
    b = np.fromfile('phys_f001_p%02d.bin' % (p+1), '<f4').astype(np.float64)
    e = np.abs(b - a*fac[p]).max()/np.abs(a*fac[p]).max()
    assert e < 1e-6, (p, e)
assert 'eV' in open('phys.ps', errors='ignore').read()
EOF
then ok "slice_plot physical units (B x bs, T x temps)"; else ng "slice_plot physical units"; fi
if command -v gs > /dev/null; then
  gs -q -dNOPAUSE -dBATCH -sDEVICE=png16m -sPAPERSIZE=a4 -r80 -sOutputFile=new.png new.ps && echo "  PNG: $W/s/new.png"
fi

#------------------------------------------------------------------ mag3d
echo "== mag3d (3D: glTF / X3D)"
cd "$W/s"
cat > mag3d.nml <<'EOF2'
&mag3d
  datafile(1) = 'ju02nt171m.data', title(1) = 'test', outbase = 'test3d', out_glb = .true., out_x3d = .true.,
  nx2 = 302, ny2 = 202, nz2 = 102, nxp = 50, h = 1.5, stride = 2, r_fine = 12.0,
  fl_nlat = 5, fl_nmlt = 12, fl_equator = .true., eq_nx = 4, eq_ny = 3,
  pl_plane(1) = 'xz', pl_qty(1) = 'temp', pl_plane(2) = 'xy', pl_qty(2) = 'log10:rho', pl_alpha_mode(2) = 'value',
  iso_qty(1) = 'absb', iso_level(1) = 1.0e-4, arr_plane = 'xy', arr_vec = 'v', arr_n = 10,
  stk_plane = 'yz', stk_n = 3, stk_lo = -20.0, stk_hi = 60.0, stk_qty = 'temp',
  sel_plane(1) = 'xz', sel_n(1) = 4, sel_lo(1) = -20.0, sel_hi(1) = 20.0, sel_qty(1) = 'absb',
  fl_box = .true., fl_box_lo = -40.0, -30.0, -10.0, fl_box_hi = 60.0, 30.0, 10.0, fl_box_n = 3, 3, 2, fl_keep = 'unconnected',
/
EOF2
if "$BIN/mag3d" mag3d.nml > mag3d.log 2>&1 && [ "$(head -c 4 test3d.glb)" = "glTF" ] && grep -q "</X3D>" test3d.x3d \
   && [ -s test3d_tex1.png ]; then
  ok "mag3d wrote test3d.glb / test3d.x3d ($(grep -o 'field lines: [0-9]*' mag3d.log))"
else
  ng "mag3d"; tail -5 mag3d.log
fi
if python3 - <<'EOF'
import json, struct
d = open('test3d.glb', 'rb').read()
j = json.loads(d[20:20 + struct.unpack('<I', d[12:16])[0]])
ex = j['scenes'][0]['extras']
names = [n['name'] for n in j['nodes']]
sel = [n for n in j['nodes'] if n.get('extras', {}).get('group') == 'sel1']
assert len(sel) == 4 and sorted(n['extras']['index'] for n in sel) == [0, 1, 2, 3]
assert any(L['target'] == 'plane01' for L in ex['legends']) and any(L['target'] == 'sel1_' for L in ex['legends'])
assert any(i['mimeType'] == 'image/jpeg' for i in j['images'])          # planet image
assert any(t['text'].endswith('[RJ]') for t in ex['labels'])
EOF
then ok "mag3d: selectable slices, colour bars, planet image, length unit"; else ng "mag3d extras"; fi
if command -v node > /dev/null && node -e "require('gltf-validator')" 2>/dev/null; then
  node -e "require('gltf-validator').validateBytes(new Uint8Array(require('fs').readFileSync('test3d.glb'))).then(r=>{console.log('  glTF validator: errors',r.issues.numErrors,'warnings',r.issues.numWarnings);process.exit(r.issues.numErrors?1:0)})" \
    && ok "glTF validation" || ng "glTF validation"
fi

#------------------------------------------------------------------ polar_map
echo "== polar_map (comparison with analytic solutions, needs python3+numpy)"
cd "$W/.."
if python3 -c "import numpy" 2>/dev/null; then
  if python3 polar_test.py "$BIN/polar_map"; then ok "polar_map analytic checks"; else ng "polar_map"; fi
else
  echo "  [SKIP] polar_map (no numpy)"
fi

echo
echo "PASS: $pass   FAIL: $fail"
[ $fail -eq 0 ]
