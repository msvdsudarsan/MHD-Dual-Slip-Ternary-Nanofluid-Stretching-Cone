"""
ternary_cone_master.py

Reference Python implementation accompanying the paper
"MHD Dual-Slip Flow and Entropy Generation in a Cu-Al2O3-TiO2/Water
Ternary Hybrid Nanofluid over a Stretching Cone".

Cu-Al2O3-TiO2/water ternary hybrid nanofluid over a stretching cone: MHD,
velocity + thermal slip, wall suction, linearised Rosseland radiation,
uniform volumetric heat source, entropy generation.

Governing system:
    4*A2*f''' + A1*(f*f'' - 0.5*f'^2) - A3*M*f' = 0
    4*(A5+Rd)*th'' + Pr*A4*f*th' + Pr*Q*th + Br*(A2*fpp^2 + A3*M*fp^2) = 0
Boundary conditions:
    f(0) = S,  f'(0) = lam + L1*f''(0),  th(0) = 1 + L2*th'(0)
    f'(inf) = 0,  th(inf) = 0

This module is the single computational record used to generate every
numerical value, table, and figure reported in the manuscript. The
boundary-value problem is solved with scipy.integrate.solve_bvp, a
fourth-order collocation method with adaptive mesh refinement. Solution
accuracy is established through mesh/tolerance refinement, recovery of
the closed-form clear-fluid limit, and comparison of the numerically
extrapolated far-field wall values against the analytical decay rates
of Proposition 1 -- see REPRODUCIBILITY_REPORT.md and Section 3 of the
manuscript.

Run:  python ternary_cone_master.py
"""

from __future__ import annotations
import numpy as np
from scipy.integrate import solve_bvp
from dataclasses import dataclass, field, replace
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================== #
#  Parameter container
# ============================================================== #
@dataclass
class Params:
    p1: float = 0.10   # Al2O3 volume fraction
    p2: float = 0.10   # Cu volume fraction
    p3: float = 0.10   # TiO2 volume fraction
    M: float = 0.6      # magnetic parameter
    Pr: float = 6.2     # Prandtl number (water)
    Rd: float = 0.1     # radiation parameter
    Q: float = 0.05     # heat source/sink parameter
    lam: float = 0.2    # stretching-rate ratio c/a
    l1: float = 0.2     # dimensionless velocity slip (L1)
    l2: float = 0.2     # dimensionless thermal slip  (L2)
    S: float = 1.0      # wall mass transfer, S > 0 = suction
    Br: float = 1.0     # Brinkman number
    Om: float = 2.0     # dimensionless temperature-difference parameter
    n: float = 4.8      # Hamilton-Crosser shape factor
    sig3: float = 2.6e6  # sigma_TiO2 [S/m]
    einf: float = 20.0   # truncation domain
    npts: int = 320      # initial mesh points
    tol: float = 1e-8
    coupleDissipation: bool = True


@dataclass
class Ratios:
    A1: float
    A2: float
    A3: float
    A4: float
    A5: float
    A7: float


# ============================================================== #
#  Thermophysical property ratios (three-stage sequential mixture)
# ============================================================== #
def maxwell(sh, sp, phi):
    """One Maxwell electrical-conductivity stage."""
    return sh * (sp + 2 * sh - 2 * phi * (sh - sp)) / (sp + 2 * sh + phi * (sh - sp))


def hamilton_crosser(kh, kp, phi, n):
    """One Hamilton-Crosser thermal-conductivity stage (n = 3 recovers Maxwell)."""
    return kh * (kp + (n - 1) * kh - (n - 1) * phi * (kh - kp)) / (kp + (n - 1) * kh + phi * (kh - kp))


def ratios(P: Params) -> Ratios:
    # Base fluid = water. Sources: CRC Handbook of Chemistry and Physics (2020)
    # and standard nanofluid property tables (Table 1 of the manuscript).
    rf, r1, r2, r3 = 997.1, 3970.0, 8933.0, 4250.0        # kg/m^3
    cf, c1, c2, c3 = 4179.0, 765.0, 385.0, 686.2          # J/(kg K)
    kf, k1, k2, k3 = 0.613, 40.0, 400.0, 8.9538           # W/(m K)
    sf, s1, s2, s3 = 0.05, 3.5e7, 5.96e7, P.sig3          # S/m
    bf, b1, b2, b3 = 21e-5, 0.85e-5, 1.67e-5, 0.9e-5      # 1/K

    p1, p2, p3, n = P.p1, P.p2, P.p3, P.n

    # --- density (stage-wise) ---
    rnf = (1 - p1) * rf + p1 * r1
    rhnf = (1 - p2) * rnf + p2 * r2
    rthnf = (1 - p3) * rhnf + p3 * r3
    A1 = rthnf / rf

    # --- viscosity (Brinkman, stage-wise) ---
    A2 = 1.0 / ((1 - p1) ** 2.5 * (1 - p2) ** 2.5 * (1 - p3) ** 2.5)

    # --- electrical conductivity (Maxwell, stage-wise) ---
    snf = maxwell(sf, s1, p1)
    shnf = maxwell(snf, s2, p2)
    sthnf = maxwell(shnf, s3, p3)
    A3 = sthnf / sf

    # --- volumetric heat capacity (stage-wise) ---
    Cf_ = rf * cf
    Cnf = (1 - p1) * Cf_ + p1 * r1 * c1
    Chnf = (1 - p2) * Cnf + p2 * r2 * c2
    Cthnf = (1 - p3) * Chnf + p3 * r3 * c3
    A4 = Cthnf / Cf_

    # --- thermal conductivity (Hamilton-Crosser, stage-wise) ---
    knf = hamilton_crosser(kf, k1, p1, n)
    khnf = hamilton_crosser(knf, k2, p2, n)
    kthnf = hamilton_crosser(khnf, k3, p3, n)
    A5 = kthnf / kf

    # --- thermal expansion (stage-wise) ---
    Bf = rf * bf
    Bnf = (1 - p1) * Bf + p1 * r1 * b1
    Bhnf = (1 - p2) * Bnf + p2 * r2 * b2
    Bthnf = (1 - p3) * Bhnf + p3 * r3 * b3
    A7 = Bthnf / Bf

    return Ratios(A1=A1, A2=A2, A3=A3, A4=A4, A5=A5, A7=A7)


# ============================================================== #
#  ODE system and boundary conditions
# ============================================================== #
def odes(eta, y, R: Ratios, P: Params):
    f, fp, fpp, th, thp = y
    fppp = -(R.A1 * (f * fpp - 0.5 * fp ** 2) - R.A3 * P.M * fp) / (4 * R.A2)
    if P.coupleDissipation:
        diss = P.Br * (R.A2 * fpp ** 2 + R.A3 * P.M * fp ** 2)
    else:
        diss = 0.0
    thpp = -(P.Pr * R.A4 * f * thp + P.Pr * P.Q * th + diss) / (4 * (R.A5 + P.Rd))
    return np.vstack([fp, fpp, fppp, thp, thpp])


def bcs(ya, yb, P: Params):
    return np.array([
        ya[0] - P.S,
        ya[1] - P.lam - P.l1 * ya[2],
        ya[3] - 1.0 - P.l2 * ya[4],
        yb[1],
        yb[3],
    ])


# ============================================================== #
#  Solver (continuation in the truncation domain, like bvpxtend)
# ============================================================== #
def solve_full(P: Params) -> dict:
    R = ratios(P)

    out = {"ok": True, "R": R}

    e0 = min(5.0, P.einf)
    n0 = max(40, round(P.npts * e0 / P.einf))
    x = np.linspace(0, e0, n0)
    y0 = np.vstack([
        P.S + P.lam * (1 - np.exp(-x)),
        P.lam * np.exp(-x),
        -P.lam * np.exp(-x),
        np.exp(-x),
        -np.exp(-x),
    ])

    try:
        sol = solve_bvp(lambda t, y: odes(t, y, R, P),
                         lambda a, b: bcs(a, b, P),
                         x, y0, tol=P.tol, max_nodes=200000, verbose=0)
        if not sol.success:
            raise RuntimeError(sol.message)

        e = e0
        while e < P.einf - 1e-12:
            e = min(e + 2.5, P.einf)
            # extend the mesh/solution guess to the larger domain,
            # continuation strategy: extend the mesh/solution guess to the larger domain
            x_ext = np.linspace(0, e, max(len(sol.x), int(round(P.npts * e / P.einf))))
            y_ext = sol.sol(np.clip(x_ext, x_ext.min(), sol.x.max()))
            sol = solve_bvp(lambda t, y: odes(t, y, R, P),
                             lambda a, b: bcs(a, b, P),
                             x_ext, y_ext, tol=P.tol, max_nodes=200000, verbose=0)
            if not sol.success:
                raise RuntimeError(sol.message)
    except Exception as exc:  # noqa: BLE001
        out.update(ok=False, msg=str(exc), finf=np.nan, Qc=np.nan,
                    admissible=False, Cf=np.nan, Nu=np.nan, Ns0=np.nan,
                    Be0=np.nan, monotone=False, oscillatory=True)
        return out

    eta = np.linspace(0, P.einf, 4001)
    Y = sol.sol(eta)
    f, fp, fpp, th, thp = Y

    out.update(eta=eta, f=f, fp=fp, fpp=fpp, th=th, thp=thp, sol=sol)

    finf = f[-1]
    Qc = P.Pr * R.A4 ** 2 * finf ** 2 / (16 * (R.A5 + P.Rd))
    admissible = (finf > 0) and (P.Q <= Qc)
    out.update(finf=finf, Qc=Qc, admissible=admissible)

    disc = (P.Pr * R.A4 * finf) ** 2 - 16 * (R.A5 + P.Rd) * P.Pr * P.Q
    if disc >= 0:
        mth = (-P.Pr * R.A4 * finf - np.sqrt(disc)) / (8 * (R.A5 + P.Rd))
        oscillatory = False
    else:
        mth = (-P.Pr * R.A4 * finf) / (8 * (R.A5 + P.Rd))
        oscillatory = True
    mf = (-R.A1 * finf - np.sqrt((R.A1 * finf) ** 2 + 16 * R.A2 * R.A3 * P.M)) / (8 * R.A2)
    out.update(mth=mth, mf=mf, oscillatory=oscillatory)

    Cf = -R.A2 * fpp[0]
    Nu = -R.A5 * (1 + P.Rd) * thp[0]
    out.update(Cf=Cf, Nu=Nu)

    Ns = (R.A5 + P.Rd) * thp ** 2 + (P.Br / P.Om) * (R.A2 * fpp ** 2 + R.A3 * P.M * fp ** 2)
    Be = (R.A5 + P.Rd) * thp ** 2 / Ns
    out.update(Ns=Ns, Be=Be, Ns0=Ns[0], Be0=Be[0])

    monotone = bool(np.all(np.diff(th) <= 1e-8)) and (np.max(th) <= 1 + 1e-6)
    out.update(monotone=monotone)
    return out


# ============================================================== #
#  Reporting helpers
# ============================================================== #
def baseline() -> Params:
    return Params()


def report(P: Params, out: dict):
    if not out["ok"]:
        print(f"  SOLVE FAILED: {out['msg']}")
        return
    print(f"  f(einf)          = {out['finf']:12.6f}")
    print(f"  Qc (critical Q)  = {out['Qc']:12.6f}     (Q = {P.Q:.4f})")
    if out["admissible"]:
        print(f"  ADMISSIBLE: yes   theta decay rate |Re m| = {abs(out['mth']):8.5f}")
        print(f"  f' decay rate     |m_f|      = {abs(out['mf']):8.5f}")
    else:
        if out["finf"] <= 0:
            print("  ADMISSIBLE: *** NO *** (f(einf) <= 0: INCREASE S)")
        else:
            print(f"  ADMISSIBLE: *** NO *** (Q > Qc: reduce Q below {out['Qc']:.6f} or increase S)")
    if out["oscillatory"]:
        print("  WARNING: complex far-field roots -> theta oscillates.")
        print("           No truncation domain will converge. Do NOT report.")
    print(f"  theta monotone   = {int(out['monotone'])}")
    print(f"  -A2*f''(0)       = {out['Cf']:12.6f}")
    print(f"  -A5(1+Rd)th'(0)  = {out['Nu']:12.6f}")
    print(f"  Ns(0)            = {out['Ns0']:12.6f}")
    print(f"  Be(0)            = {out['Be0']:12.6f}")


# ============================================================== #
#  Table generators (reproduce every table reported in the paper)
# ============================================================== #
def admissibility_table(P0: Params):
    print("\n================ TABLE: ADMISSIBILITY (tab:admiss) ================")
    print(f"{'S':>8} {'Q':>8} {'f(einf)':>12} {'Qc':>12} {'verdict':>10}")
    Svals = [-0.5, 0, 0.5, 1.0, 1.5, 2.0]
    Qvals = [-0.1, 0, 0.05, 0.10, 0.20]
    for S in Svals:
        for Q in Qvals:
            P = replace(P0, S=S, Q=Q, einf=15.0)
            o = solve_full(P)
            if not o["ok"]:
                v, fi, qc = "fail", np.nan, np.nan
            else:
                fi, qc = o["finf"], o["Qc"]
                v = "YES" if o["admissible"] else "no"
            print(f"{S:8.2f} {Q:8.2f} {fi:12.6f} {qc:12.6f} {v:>10}")


def truncation_table(P0: Params):
    print("\n================ TABLE: TRUNCATION (tab:trunc) ================")
    print("{:>6} {:>6} {:>12} {:>16} {:>12} {:>12} {:>6}".format("p3","einf","-A2f''(0)","-A5(1+Rd)th'(0)","Ns(0)","Be(0)","adm"))
    for p3 in [0, P0.p3]:
        for e in [5, 8, 10, 12, 15, 18, 20, 25, 30]:
            P = replace(P0, p3=p3, einf=float(e), npts=max(160, round(16 * e)))
            o = solve_full(P)
            if not o["ok"]:
                continue
            print(f"{p3:6.2f} {e:6.1f} {o['Cf']:12.6f} {o['Nu']:16.6f} {o['Ns0']:12.6f} {o['Be0']:12.6f} {int(o['admissible']):6d}")


def mesh_table(P0: Params):
    print("\n================ TABLE: MESH INDEPENDENCE (tab:mesh) ================")
    print("{:>8} {:>10} {:>14} {:>16} {:>12} {:>12}".format("Npts","tol","-A2f''(0)","-A5(1+Rd)th'(0)","dCf","dNu"))
    lev = [40, 80, 160, 320, 640, 1280]
    tol = [1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-10]
    Cfp, Nup = np.nan, np.nan
    for npts, tl in zip(lev, tol):
        P = replace(P0, npts=npts, tol=tl)
        o = solve_full(P)
        if not o["ok"]:
            continue
        dC, dN = abs(o["Cf"] - Cfp), abs(o["Nu"] - Nup)
        print(f"{npts:8d} {tl:10.1e} {o['Cf']:14.8f} {o['Nu']:16.8f} {dC:12.2e} {dN:12.2e}")
        Cfp, Nup = o["Cf"], o["Nu"]


def validation_table(P0: Params):
    print("\n================ TABLE: VALIDATION (tab:valid) ================")
    print("{:>8} {:>14} {:>16}".format("lam","-f''(0)","-th'(0)"))
    for lam in [0.5, 1.0, 1.5, 2.0]:
        P = replace(P0, p1=0, p2=0, p3=0, M=0, Rd=0, Q=0, l1=0, l2=0, S=0, lam=lam)
        o = solve_full(P)
        if not o["ok"]:
            continue
        print(f"{lam:8.2f} {-o['fpp'][0]:14.6f} {-o['thp'][0]:16.6f}")


def phi3_table(P0: Params):
    print("\n================ TABLE: PHI3 SWEEP (tab:phi3) ================")
    print("{:>8} {:>12} {:>16} {:>12} {:>12} {:>10} {:>10}".format("p3","-A2f''(0)","-A5(1+Rd)th'(0)","Ns(0)","Be(0)","dCf%","dNu%"))
    p3v = [0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10]
    C0 = N0 = None
    for p3 in p3v:
        P = replace(P0, p3=p3)
        o = solve_full(P)
        if not o["ok"]:
            continue
        if C0 is None:
            C0, N0 = o["Cf"], o["Nu"]
        dC = 100 * (o["Cf"] - C0) / C0
        dN = 100 * (o["Nu"] - N0) / N0
        print(f"{p3:8.2f} {o['Cf']:12.6f} {o['Nu']:16.6f} {o['Ns0']:12.6f} {o['Be0']:12.6f} {dC:10.2f} {dN:10.2f}")


def sigma_table(P0: Params):
    print("\n============ TABLE: SIGMA_TiO2 SENSITIVITY (tab:sigsens) ============")
    print("{:>12} {:>12} {:>14} {:>16}".format("sigma3","A3","-A2f''(0)","-A5(1+Rd)th'(0)"))
    for s3 in [2.6e4, 2.6e5, 2.6e6, 2.6e7]:
        P = replace(P0, sig3=s3)
        R = ratios(P)
        o = solve_full(P)
        if not o["ok"]:
            continue
        print(f"{s3:12.2e} {R.A3:12.6f} {o['Cf']:14.6f} {o['Nu']:16.6f}")


def shape_table(P0: Params):
    print("\n================ TABLE: SHAPE FACTOR (tab:shape) ================")
    print("{:>8} {:>10} {:>14} {:>16}".format("n","A5","-A2f''(0)","-A5(1+Rd)th'(0)"))
    for n in [3.0, 3.7, 4.8, 5.7, 6.0]:
        P = replace(P0, n=n)
        R = ratios(P)
        o = solve_full(P)
        if not o["ok"]:
            continue
        print(f"{n:8.1f} {R.A5:10.6f} {o['Cf']:14.6f} {o['Nu']:16.6f}")


def sweep_table(P0: Params):
    print("\n================ TABLE: FULL SWEEP (tab:sweep) ================")
    sweeps = {
        "M": [0, 0.5, 1.0, 1.5, 2.0],
        "Rd": [0, 0.5, 1.0, 1.5, 2.0],
        "Q": [-0.2, -0.1, 0, 0.05, 0.10],
        "Pr": [4.0, 5.0, 6.2, 7.0, 8.0],
        "l1": [0, 0.2, 0.4, 0.6, 0.8],
        "l2": [0, 0.2, 0.4, 0.6, 0.8],
        "S": [0.5, 1.0, 1.5, 2.0, 2.5],
        "Br": [0.5, 1.0, 1.5, 2.0],
        "Om": [1.0, 1.5, 2.0, 2.5],
        "p1": [0, 0.02, 0.05, 0.08, 0.10],
        "p2": [0, 0.02, 0.05, 0.08, 0.10],
        "p3": [0, 0.02, 0.05, 0.08, 0.10],
    }
    for nm, vv in sweeps.items():
        print(f"\n-- {nm} --")
        print("{:>10} {:>14} {:>16} {:>12} {:>12} {:>6}".format(nm,"-A2f''(0)","-A5(1+Rd)th'(0)","Ns(0)","Be(0)","adm"))
        for v in vv:
            P = replace(P0, **{nm: v})
            o = solve_full(P)
            if not o["ok"]:
                print(f"{v:10.3f}  solve failed")
                continue
            print(f"{v:10.3f} {o['Cf']:14.6f} {o['Nu']:16.6f} {o['Ns0']:12.6f} {o['Be0']:12.6f} {int(o['admissible']):6d}")


# ============================================================== #
#  Figures
# ============================================================== #
def make_figures(P0: Params, outdir="."):
    print("\n================ FIGURES ================")
    Pb = replace(P0, p3=0)
    o0 = solve_full(Pb)
    o1 = solve_full(P0)
    if not (o0["ok"] and o1["ok"]):
        print("  skipped (solve failed)")
        return

    names = ["velocity_profile", "temperature_profile", "entropy_Ns", "bejan_number"]
    ylabels = [r"$f'(\eta)$", r"$\theta(\eta)$", r"$N_s(\eta)$", r"$Be(\eta)$"]
    Yb = [o0["fp"], o0["th"], o0["Ns"], o0["Be"]]
    Yt = [o1["fp"], o1["th"], o1["Ns"], o1["Be"]]
    xcap = [P0.einf, P0.einf, 10, 10]

    for nm, yl, yb, yt, xc in zip(names, ylabels, Yb, Yt, xcap):
        fig, ax = plt.subplots(figsize=(6.8, 5.0))
        ax.plot(o0["eta"], yb, "b--", linewidth=2.6, label=r"Binary ($\phi_3=0$)")
        ax.plot(o1["eta"], yt, "r-", linewidth=2.6, label=r"Ternary ($\phi_3=0.1$)")
        ax.set_xlabel(r"$\eta$", fontsize=16)
        ax.set_ylabel(yl, fontsize=16)
        ax.set_xlim(0, xc)
        ax.grid(True)
        ax.legend(loc="best", fontsize=12)
        fig.tight_layout()
        fig.savefig(f"{outdir}/{nm}.pdf")
        plt.close(fig)
        print(f"  saved {nm}.pdf")


def make_convergence_figure(P0: Params, outdir="."):
    print("\n================ CONVERGENCE FIGURE ================")
    elist = [5, 8, 10, 12, 15, 18, 20, 25, 30]
    Cf0, Nu0, Cf1, Nu1 = [], [], [], []
    for e in elist:
        Pb = replace(P0, p3=0, einf=float(e), npts=max(160, round(16 * e)))
        Pt = replace(P0, einf=float(e), npts=max(160, round(16 * e)))
        ob, ot = solve_full(Pb), solve_full(Pt)
        Cf0.append(ob["Cf"] if ob["ok"] else np.nan)
        Nu0.append(ob["Nu"] if ob["ok"] else np.nan)
        Cf1.append(ot["Cf"] if ot["ok"] else np.nan)
        Nu1.append(ot["Nu"] if ot["ok"] else np.nan)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    axes[0].plot(elist, Cf0, "bo--", label=r"Binary ($\phi_3=0$)")
    axes[0].plot(elist, Cf1, "rs-", label=r"Ternary ($\phi_3=0.1$)")
    axes[0].axvline(P0.einf, color="k", linestyle=":")
    axes[0].set_xlabel(r"$\eta_\infty$"); axes[0].set_ylabel(r"$C_f$")
    axes[0].grid(True); axes[0].legend(loc="best")

    axes[1].plot(elist, Nu0, "bo--", label=r"Binary ($\phi_3=0$)")
    axes[1].plot(elist, Nu1, "rs-", label=r"Ternary ($\phi_3=0.1$)")
    axes[1].axvline(P0.einf, color="k", linestyle=":")
    axes[1].set_xlabel(r"$\eta_\infty$"); axes[1].set_ylabel(r"$Nu$")
    axes[1].grid(True); axes[1].legend(loc="best")

    fig.tight_layout()
    fig.savefig(f"{outdir}/convergence.pdf")
    plt.close(fig)
    print("  saved convergence.pdf")


# ============================================================== #
#  Main driver (mirrors ternary_cone_master.m)
# ============================================================== #
def main():
    P = baseline()
    R = ratios(P)

    if P.coupleDissipation:
        print("*** Coupled viscous dissipation + Ohmic heating included in the energy equation. ***")
    else:
        print("*** Energy equation decoupled from dissipation (P.coupleDissipation = False). ***")

    print("================ EFFECTIVE PROPERTY RATIOS ================")
    print(f"  A1 (rho)      = {R.A1:12.6f}")
    print(f"  A2 (mu)       = {R.A2:12.6f}")
    print(f"  A3 (sigma)    = {R.A3:12.6f}")
    print(f"  A4 (rho cp)   = {R.A4:12.6f}")
    print(f"  A5 (k)        = {R.A5:12.6f}")
    print(f"  A7 (rho beta) = {R.A7:12.6f}")
    print(f"  A5 + Rd       = {R.A5 + P.Rd:12.6f}")

    print("\n================ BASELINE SOLVE ================")
    out = solve_full(P)
    report(P, out)

    admissibility_table(P)
    truncation_table(P)
    mesh_table(P)
    validation_table(P)
    phi3_table(P)
    sigma_table(P)
    shape_table(P)
    sweep_table(P)
    make_figures(P)
    make_convergence_figure(P)


if __name__ == "__main__":
    main()
