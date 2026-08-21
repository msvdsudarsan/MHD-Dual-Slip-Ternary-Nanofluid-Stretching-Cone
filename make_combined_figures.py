"""
Combines the velocity and temperature profiles into a single 2-panel
figure, and the local entropy generation and Bejan number profiles into
another 2-panel figure. Produces velocity_temperature.pdf and
entropy_bejan.pdf, using the same solve_full() baseline/binary
computation as make_figures() in ternary_cone_master.py.
"""
import sys
sys.path.insert(0, '.')
from dataclasses import replace
import matplotlib.pyplot as plt
from ternary_cone_master import Params, baseline, solve_full

def make_combined_figures(P0: Params, outdir="."):
    Pb = replace(P0, p3=0)
    o0 = solve_full(Pb)
    o1 = solve_full(P0)
    if not (o0["ok"] and o1["ok"]):
        print("skipped (solve failed)")
        return

    # Panel 1: velocity + temperature
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    axes[0].plot(o0["eta"], o0["fp"], "b--", linewidth=2.2, label=r"Binary ($\phi_3=0$)")
    axes[0].plot(o1["eta"], o1["fp"], "r-", linewidth=2.2, label=r"Ternary ($\phi_3=0.1$)")
    axes[0].set_xlabel(r"$\eta$", fontsize=14)
    axes[0].set_ylabel(r"$f'(\eta)$", fontsize=14)
    axes[0].set_xlim(0, P0.einf)
    axes[0].grid(True)
    axes[0].legend(loc="best", fontsize=10)
    axes[0].set_title("(a) Velocity profile", fontsize=12)

    axes[1].plot(o0["eta"], o0["th"], "b--", linewidth=2.2, label=r"Binary ($\phi_3=0$)")
    axes[1].plot(o1["eta"], o1["th"], "r-", linewidth=2.2, label=r"Ternary ($\phi_3=0.1$)")
    axes[1].set_xlabel(r"$\eta$", fontsize=14)
    axes[1].set_ylabel(r"$\theta(\eta)$", fontsize=14)
    axes[1].set_xlim(0, P0.einf)
    axes[1].grid(True)
    axes[1].legend(loc="best", fontsize=10)
    axes[1].set_title("(b) Temperature profile", fontsize=12)

    fig.tight_layout()
    fig.savefig(f"{outdir}/velocity_temperature.pdf")
    plt.close(fig)
    print("saved velocity_temperature.pdf")

    # Panel 2: entropy + Bejan number
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    axes[0].plot(o0["eta"], o0["Ns"], "b--", linewidth=2.2, label=r"Binary ($\phi_3=0$)")
    axes[0].plot(o1["eta"], o1["Ns"], "r-", linewidth=2.2, label=r"Ternary ($\phi_3=0.1$)")
    axes[0].set_xlabel(r"$\eta$", fontsize=14)
    axes[0].set_ylabel(r"$N_s(\eta)$", fontsize=14)
    axes[0].set_xlim(0, 10)
    axes[0].grid(True)
    axes[0].legend(loc="best", fontsize=10)
    axes[0].set_title("(a) Local entropy generation", fontsize=12)

    axes[1].plot(o0["eta"], o0["Be"], "b--", linewidth=2.2, label=r"Binary ($\phi_3=0$)")
    axes[1].plot(o1["eta"], o1["Be"], "r-", linewidth=2.2, label=r"Ternary ($\phi_3=0.1$)")
    axes[1].set_xlabel(r"$\eta$", fontsize=14)
    axes[1].set_ylabel(r"$Be(\eta)$", fontsize=14)
    axes[1].set_xlim(0, 10)
    axes[1].grid(True)
    axes[1].legend(loc="best", fontsize=10)
    axes[1].set_title("(b) Bejan number", fontsize=12)

    fig.tight_layout()
    fig.savefig(f"{outdir}/entropy_bejan.pdf")
    plt.close(fig)
    print("saved entropy_bejan.pdf")


if __name__ == "__main__":
    P0 = baseline()
    make_combined_figures(P0, outdir=".")
