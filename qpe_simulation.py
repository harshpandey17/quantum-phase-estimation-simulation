
import importlib
import json
import os
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

OUT = os.path.dirname(os.path.abspath(__file__))  # save outputs next to this script
SHOTS = 4096


def out(name: str) -> str:
    return os.path.join(OUT, name)


# ----------------------------------------------------------------------------
# QPE implementation
# ----------------------------------------------------------------------------
def inverse_qft(n: int) -> QuantumCircuit:
    """Inverse QFT on n qubits (including the final swaps)."""
    qc = QuantumCircuit(n, name="QFT†")
    for j in range(n // 2):
        qc.swap(j, n - j - 1)
    for j in range(n):
        for m in range(j):
            qc.cp(-np.pi / 2 ** (j - m), m, j)
        qc.h(j)
    return qc


def qpe_circuit(unitary: QuantumCircuit, state_prep: QuantumCircuit, n_count: int) -> QuantumCircuit:
    """
    Build a QPE circuit.

    unitary    : circuit implementing U on m qubits
    state_prep : circuit preparing the eigenstate |psi> on m qubits
    n_count    : number of counting qubits (precision)
    """
    m = unitary.num_qubits
    qc = QuantumCircuit(n_count + m, n_count)
    counting = list(range(n_count))
    target = list(range(n_count, n_count + m))

    qc.compose(state_prep, qubits=target, inplace=True)           # 1. eigenstate
    qc.h(counting)                                                 # 2. superposition
    for k in range(n_count):                                       # 3. controlled-U^(2^k)
        cu = unitary.power(2**k, matrix_power=True).to_gate(label=f"U^{2**k}").control(1)
        qc.append(cu, [counting[k]] + target)
    qc.append(inverse_qft(n_count).to_gate(), counting)            # 4. inverse QFT
    qc.measure(counting, range(n_count))                           # 5. measure
    return qc


def run_qpe(unitary, state_prep, n_count, shots=SHOTS):
    """Simulate QPE; return (estimated theta, counts, circuit)."""
    qc = qpe_circuit(unitary, state_prep, n_count)
    sim = AerSimulator()
    counts = sim.run(transpile(qc, sim), shots=shots).result().get_counts()
    best = max(counts, key=counts.get)
    return int(best, 2) / 2**n_count, counts, qc


def phase_unitary(theta: float) -> QuantumCircuit:
    """U = diag(1, e^{2*pi*i*theta}); |1> is the eigenstate with phase theta."""
    u = QuantumCircuit(1)
    u.p(2 * np.pi * theta, 0)
    return u


def theory_distribution(theta: float, n: int) -> np.ndarray:
    """Analytic probability of measuring integer m for phase theta."""
    N = 2**n
    p = np.zeros(N)
    for m in range(N):
        d = theta - m / N
        p[m] = 1.0 if abs(d) < 1e-12 else abs(np.sin(np.pi * N * d) / (N * np.sin(np.pi * d))) ** 2
    return p


def counts_to_probs(counts: dict, n: int) -> np.ndarray:
    p = np.zeros(2**n)
    total = sum(counts.values())
    for bits, v in counts.items():
        p[int(bits, 2)] = v / total
    return p


def draw_circuit_image(circ: QuantumCircuit, path: str):
    """Save the circuit as a PNG. Installs pylatexenc automatically if it is missing."""
    try:
        importlib.import_module("pylatexenc")
    except ImportError:
        print("Installing pylatexenc (needed for circuit images)...")
        subprocess.call([sys.executable, "-m", "pip", "install", "pylatexenc", "-q"])
        importlib.invalidate_caches()
    try:
        fig = circ.draw("mpl", fold=-1)
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"Saved circuit diagram: {path}")
    except Exception as e:
        print(f"Could not draw image circuit ({e}). Text version:")
        print(circ.draw("text", fold=120))


# ----------------------------------------------------------------------------
# Experiments, plots, structured results
# ----------------------------------------------------------------------------
def main():
    prep = QuantumCircuit(1)
    prep.x(0)  # eigenstate |1>
    n = 5
    results = {}

    # Experiment 1: exact phase (T gate, theta = 1/8)
    th1, c1, _ = run_qpe(phase_unitary(1 / 8), prep, n)
    results["exact_theta_1_8"] = dict(n_count=n, true_theta=1 / 8, estimated_theta=th1,
                                      error=abs(th1 - 1 / 8), counts=c1)

    # Experiment 2: inexact phase (theta = 1/3)
    th2, c2, _ = run_qpe(phase_unitary(1 / 3), prep, n)
    results["inexact_theta_1_3"] = dict(n_count=n, true_theta=1 / 3, estimated_theta=th2,
                                        error=abs(th2 - 1 / 3), counts=c2)

    # Experiment 3: precision scaling
    scaling = []
    for k in range(2, 9):
        t, c, _ = run_qpe(phase_unitary(1 / 3), prep, k)
        scaling.append(dict(n_count=k, estimate=t, error=abs(t - 1 / 3),
                            resolution=1 / 2**k, top_prob=max(c.values()) / SHOTS))
    results["scaling_theta_1_3"] = scaling

    # Plot 1: histograms
    N = 2**n
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax, key, title, true in [(axes[0], "exact_theta_1_8", "Exact phase θ = 1/8 (T gate)", 1 / 8),
                                 (axes[1], "inexact_theta_1_3", "Inexact phase θ = 1/3", 1 / 3)]:
        probs = counts_to_probs(results[key]["counts"], n)
        ax.bar(range(N), probs, color="#4C72B0")
        ax.axvline(true * N, color="k", ls="--", lw=1, label=f"true θ·2ⁿ = {true * N:.2f}")
        ax.set_xticks(range(0, N, 4))
        ax.set_xticklabels([f"{m}\n{m / N:.3f}" for m in range(0, N, 4)], fontsize=7)
        ax.set_xlabel("measured integer m  /  phase m/2ⁿ")
        ax.set_ylabel("Probability")
        ax.set_title(title)
        ax.legend()
    plt.tight_layout()
    plt.savefig(out("qpe_histograms.png"), dpi=150)

    # Plot 2: simulation vs theory
    emp = counts_to_probs(c2, n)
    th_p = theory_distribution(1 / 3, n)
    x = np.arange(N) / N
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(x - 0.005, emp, width=0.01, label="Simulated", color="#4C72B0")
    ax.bar(x + 0.005, th_p, width=0.01, label="Theory", color="#DD8452")
    ax.axvline(1 / 3, color="k", ls="--", lw=1, label="true θ = 1/3")
    ax.set_xlabel("estimated phase m/2ⁿ")
    ax.set_ylabel("Probability")
    ax.set_title("QPE outcome distribution (n = 5): simulation vs theory")
    ax.legend()
    plt.tight_layout()
    plt.savefig(out("qpe_theory_vs_sim.png"), dpi=150)
    results["tv_distance_sim_vs_theory"] = float(0.5 * np.abs(emp - th_p).sum())

    # Plot 3: precision scaling
    ks = [s["n_count"] for s in scaling]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogy(ks, [s["error"] for s in scaling], "o-", label="|θ̂ − θ| (most likely outcome)")
    ax.semilogy(ks, [s["resolution"] / 2 for s in scaling], "s--", label="Bound 1/2ⁿ⁺¹")
    ax.set_xlabel("Counting qubits n")
    ax.set_ylabel("Error")
    ax.set_title("Precision scaling, θ = 1/3")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out("qpe_precision_scaling.png"), dpi=150)

    # Plot 4: circuit diagram (n = 3)
    circ = qpe_circuit(phase_unitary(1 / 8), prep, 3)
    draw_circuit_image(circ, out("qpe_circuit_n3.png"))

    # Circuit statistics
    t = transpile(qpe_circuit(phase_unitary(1 / 3), prep, 5), AerSimulator(), optimization_level=1)
    results["circuit_n5"] = dict(qubits=t.num_qubits, depth=t.depth(), ops=dict(t.count_ops()))

    # Structured results
    with open(out("qpe_results.json"), "w") as f:
        json.dump(results, f, indent=2, default=float)

    # Console summary
    print(f"Exact   (θ=1/8): estimate = {th1:.5f}, error = {abs(th1 - 1 / 8):.5f}")
    print(f"Inexact (θ=1/3): estimate = {th2:.5f}, error = {abs(th2 - 1 / 3):.5f}")
    print(f"TV distance (sim vs theory): {results['tv_distance_sim_vs_theory']:.4f}")
    print("\nn_count  estimate   error     top_prob")
    for s in scaling:
        print(f"{s['n_count']:>7}  {s['estimate']:.6f}  {s['error']:.6f}  {s['top_prob']:.3f}")
    print(f"\nSaved outputs to: {OUT}")
    plt.show()  # open all plot windows (close them to finish)


if __name__ == "__main__":
    main()