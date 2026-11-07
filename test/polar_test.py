#!/usr/bin/env python3
"""Test of polar_map (needs python3 + numpy)

Builds synthetic data with known analytic solutions (half domain z>=0, m.data format) and checks polar_map.
  B : dipole (b0=-1) + uniform field B0 z (antiparallel at the equator -> open polar cap) + a localized twist
  v : rigid corotation + c grad(L)/|grad L|^2  (dL/dt = c)
Checks
  1. field-aligned current j_mhd : value from the analytic rot B
  2. mapped convection u_colat   : -c / (2 L sqrt(L-1))   (analytic dipole mapping)
  3. open-closed boundary        : boundary from tracing the same field in python (13.95 deg)
  4. corotation only: u_mlt = Omega sin(th), potential = Omega (sin^2 th - sin^2 th_max)
  5. northern and southern (mirrored) results agree
Usage: python3 polar_test.py <polar_map executable>
"""
import subprocess, sys, os
import numpy as np

exe = os.path.abspath(sys.argv[1])
os.makedirs('work/polar', exist_ok=True)
os.chdir('work/polar')

nx2, ny2, nz2, h = 402, 402, 202, 0.3
b0, B0, alpha, w, Om = -1.0, 6.4e-5, 1e-5, 15.0, 1.38e-4
i0 = (nx2 + 1) / 2; j0 = (ny2 + 1) / 2; k0 = 1.5
x = (np.arange(1, nx2 + 1) - i0) * h
y = (np.arange(1, ny2 + 1) - j0) * h
X, Y = np.meshgrid(x, y, indexing='xy')          # (ny2, nx2) = Fortran (nx2, ny2)


def comp(m, z, c):
    r2 = X**2 + Y**2 + z * z; r = np.sqrt(r2); r5 = r2**2.5
    g = alpha * np.exp(-r2 / w**2)
    if m == 0: return np.ones_like(X)
    if m == 4: return 1e-3 / np.maximum(r, 1.0)
    if m == 5: return -3 * b0 * X * z / r5 - g * Y
    if m == 6: return -3 * b0 * Y * z / r5 + g * X
    if m == 7: return b0 * (X**2 + Y**2 - 2 * z * z) / r5 + B0
    rho2 = X**2 + Y**2
    dl = [3 * r * X / rho2 - 2 * r**3 * X / rho2**2,
          3 * r * Y / rho2 - 2 * r**3 * Y / rho2**2,
          3 * r * z / rho2]
    g2 = dl[0]**2 + dl[1]**2 + dl[2]**2
    cor = [-Om * Y, Om * X, 0 * X][m - 1]
    return cor + c * dl[m - 1] / g2


def gen(fname, c):
    with open(fname, 'wb') as fo:
        for m in range(8):
            for k in range(1, nz2 + 1):
                a = comp(m, (k - k0) * h, c).astype('<f4').tobytes()
                n = np.array([len(a)], '<i4').tobytes()
                fo.write(n + a + n)


def run(name, data, hem, corot):
    with open(name + '.nml', 'w') as f:
        f.write(f"""&polar
 datafile(1)='{data}', nx2={nx2}, ny2={ny2}, nz2={nz2}, nxp=0, h={h}, half_domain=.true.,
 hemisphere='{hem}', r_map=10.0, colat_max=20.0, ncolat=80, nmlt=96, stride=1,
 subtract_corot={corot}, psfile='{name}.ps', jpar_avg=0.0
/
&units phys_units=.false. /
""")
    subprocess.run([exe, name + '.nml'], check=True, stdout=subprocess.DEVNULL)
    return np.fromfile(name + '_001.bin', '<f4').reshape(12, 96, 80).transpose(2, 1, 0)


gen('t1m.data', 1.0)
gen('t2m.data', 0.0)
N = run('n1', 't1m.data', 'north', '.true.')
S = run('s1', 't1m.data', 'south', '.true.')
C = run('c2', 't2m.data', 'north', '.false.')

nc, nm = 80, 96
dth = np.deg2rad(20 / nc); th = (np.arange(nc) + 0.5) * dth; tm = np.deg2rad(20)
ok = True


def check(label, err, tol):
    global ok
    res = 'PASS' if err <= tol else 'FAIL'
    if err > tol: ok = False
    print(f'  [{res}] {label}: {err:.2e} (tol {tol:g})')


# 1. j_mhd
def Bf(p):
    x_, y_, z_ = p; r2 = x_*x_ + y_*y_ + z_*z_; r5 = r2**2.5; g = alpha*np.exp(-r2/w**2)
    return np.array([-3*b0*x_*z_/r5 - g*y_, -3*b0*y_*z_/r5 + g*x_, b0*(x_*x_ + y_*y_ - 2*z_*z_)/r5 + B0])


def curl(p, d=1e-3):
    dB = [(Bf(p + d*e) - Bf(p - d*e)) / (2*d) for e in np.eye(3)]
    return np.array([dB[1][2] - dB[2][1], dB[2][0] - dB[0][2], dB[0][1] - dB[1][0]])


err, ref = 0.0, 0.0
for i in range(0, nc, 4):
    for j in range(0, nm, 6):
        ell = 1 / np.sin(th[i])**2; rm = min(10.0, ell)
        if rm < 8: continue
        thm = np.arcsin(min(np.sqrt(rm/ell), 1)); ph = j*2*np.pi/nm - np.pi
        p = np.array([-rm*np.sin(thm)*np.cos(ph), -rm*np.sin(thm)*np.sin(ph), rm*np.cos(thm)])
        if p[2] < 1: continue
        Bm = Bf(p); jexp = curl(p) @ Bm / (Bm @ Bm) * (-2*b0*np.cos(th[i]))
        err = max(err, abs(N[i, j, 0] - jexp)); ref = max(ref, abs(jexp))
check('j_mhd relative error', err/ref, 1e-2)
# 2. mapped flow
L = 1/np.sin(th)**2; uexp = -1.0/(2*L*np.sqrt(L - 1))
check('u_colat relative error (colat>=2deg)', np.abs(N[8:, :, 3] - uexp[8:, None]).max()/np.abs(uexp).max(), 1e-2)
check('u_mlt (should be 0)', np.abs(N[:, :, 4]).max()/np.abs(uexp).max(), 1e-3)
# 3. open-closed boundary
ocb = [np.rad2deg(th[np.where(N[:, j, 5] == 0)[0][0]]) for j in range(nm)]
check('OCB position error [deg]', abs(np.mean(ocb) - 13.95), 0.25)
# 4. corotation
check('corotation u_mlt relative error', np.abs(C[:, :, 4] - Om*np.sin(th)[:, None]).max()/(Om*np.sin(tm)), 1e-3)
phi = Om*(np.sin(th)**2 - np.sin(tm)**2)
check('corotation potential relative error', np.abs(C[:, :, 2] - phi[:, None]).max()/np.abs(phi).max(), 1e-3)
# 5. hemispheres
f = np.isfinite(N) & np.isfinite(S)
check('north/south mirror difference (relative)', np.abs(N[f] - S[f]).max()/np.abs(N[f]).max(), 1e-5)
print('polar_map test:', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
