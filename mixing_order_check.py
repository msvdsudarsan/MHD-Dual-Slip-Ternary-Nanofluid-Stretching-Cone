import sys
sys.path.insert(0, '.')
from itertools import permutations
import numpy as np
from dataclasses import replace
from ternary_cone_master import Params, Ratios, maxwell, hamilton_crosser, solve_full, baseline

def ratios_ordered(P: Params, order) -> Ratios:
    rf, r1, r2, r3 = 997.1, 3970.0, 8933.0, 4250.0
    cf, c1, c2, c3 = 4179.0, 765.0, 385.0, 686.2
    kf, k1, k2, k3 = 0.613, 40.0, 400.0, 8.9538
    sf, s1, s2, s3 = 0.05, 3.5e7, 5.96e7, P.sig3
    bf, b1, b2, b3 = 21e-5, 0.85e-5, 1.67e-5, 0.9e-5
    phis = {1: P.p1, 2: P.p2, 3: P.p3}
    rp = {1: r1, 2: r2, 3: r3}
    cp = {1: c1, 2: c2, 3: c3}
    kp = {1: k1, 2: k2, 3: k3}
    sp = {1: s1, 2: s2, 3: s3}
    bp = {1: b1, 2: b2, 3: b3}
    n = P.n

    # density & heat capacity & viscosity: order-independent (linear/power-law) -- compute normally
    p1, p2, p3 = P.p1, P.p2, P.p3
    rnf = (1-p1)*rf + p1*r1; rhnf=(1-p2)*rnf+p2*r2; rthnf=(1-p3)*rhnf+p3*r3
    A1 = rthnf/rf
    A2 = 1.0/((1-p1)**2.5*(1-p2)**2.5*(1-p3)**2.5)
    Cf_=rf*cf; Cnf=(1-p1)*Cf_+p1*r1*c1; Chnf=(1-p2)*Cnf+p2*r2*c2; Cthnf=(1-p3)*Chnf+p3*r3*c3
    A4 = Cthnf/Cf_

    # electrical conductivity (Maxwell, order-dependent)
    sh = sf
    for i in order:
        sh = maxwell(sh, sp[i], phis[i])
    A3 = sh/sf

    # thermal conductivity (Hamilton-Crosser, order-dependent)
    kh = kf
    for i in order:
        kh = hamilton_crosser(kh, kp[i], phis[i], n)
    A5 = kh/kf

    # thermal expansion: order-independent
    Bf=rf*bf; Bnf=(1-p1)*Bf+p1*r1*b1; Bhnf=(1-p2)*Bnf+p2*r2*b2; Bthnf=(1-p3)*Bhnf+p3*r3*b3
    A7 = Bthnf/Bf

    return Ratios(A1=A1,A2=A2,A3=A3,A4=A4,A5=A5,A7=A7)


def solve_full_ordered(P, order):
    import ternary_cone_master as tcm
    R = ratios_ordered(P, order)
    # monkeypatch ratios function temporarily
    orig = tcm.ratios
    tcm.ratios = lambda PP: R
    try:
        out = tcm.solve_full(P)
    finally:
        tcm.ratios = orig
    return out, R

P0 = baseline()
labels = {1:"Al2O3",2:"Cu",3:"TiO2"}
print(f"{'Order':<20}{'A3':>12}{'A5':>12}{'Cf':>14}{'Nu':>14}{'Ns0':>12}{'Be0':>10}")
results = []
for order in permutations([1,2,3]):
    out, R = solve_full_ordered(P0, order)
    name = "-".join(labels[i] for i in order)
    print(f"{name:<20}{R.A3:12.6f}{R.A5:12.6f}{out['Cf']:14.6f}{out['Nu']:14.6f}{out['Ns0']:12.6f}{out['Be0']:10.6f}")
    results.append((name, R.A3, R.A5, out['Cf'], out['Nu'], out['Ns0'], out['Be0']))

Cfs = [r[3] for r in results]; Nus=[r[4] for r in results]
print(f"\nCf range: {min(Cfs):.6f} to {max(Cfs):.6f}  (spread {100*(max(Cfs)-min(Cfs))/min(Cfs):.3f}%)")
print(f"Nu range: {min(Nus):.6f} to {max(Nus):.6f}  (spread {100*(max(Nus)-min(Nus))/min(Nus):.3f}%)")
