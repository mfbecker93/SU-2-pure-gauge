from pathlib import Path

import numpy as np
import pandas as pd

from analysis import (
    autocorr,
    integrated_autocorr_time,
    find_optimal_W,
)


ROOT = Path(__file__).resolve().parents[3]

INPUT = ROOT / "validation" / "aderrors" / "input" / "plaquette_beta2.csv"

df = pd.read_csv(INPUT)

x = df["plaquette"].to_numpy()

N = len(x)
mean = np.mean(x)
gamma0 = np.var(x, ddof=0)

print("N      =", N)
print("mean   =", mean)
print("gamma0 =", gamma0)

print("\n--- tau scan ---")

tau_df = integrated_autocorr_time(
    df,
    value="plaquette",
    Wmax=100,
)

print(tau_df)

print("\n--- optimal W, S=4 ---")

result = find_optimal_W(
    df,
    value="plaquette",
    Wmax=100,
    S=4,
)

print(result)

tau_int = result["tau_int"]
gamma0_corrected = result["gamma_corrected"][0]

sigma = np.sqrt(
    2 * tau_int * gamma0_corrected / N
)

ac_df = autocorr(
    df,
    value="plaquette",
    max_lag=30,
)

print(ac_df)

print("\n--- final result ---")

print("tau_int =", tau_int)
print("delta_tau =", result["delta_tau"])
print("W =", result["W"])
print("tau_exp =", result["tau_exp"])
print("gamma0 =", gamma0)
print("sigma =", sigma)


x_centered = x - np.mean(x)

print("\n--- direct autocorrelation tests ---")

for t in range(0, 10):

    numerator = np.sum(
        x_centered[:N-t] * x_centered[t:]
    )

    # Definition A: N normalization
    gamma_N = numerator / N

    # Definition B: (N-t) normalization
    gamma_Nt = numerator / (N - t)

    rho_N = gamma_N / gamma_N if t == 0 else gamma_N / (np.var(x, ddof=0))
    rho_Nt = gamma_Nt / np.var(x, ddof=0)

    print(
        f"{t:2d}  "
        f"rho_N={rho_N:.12f}  "
        f"rho_Nt={rho_Nt:.12f}"
    )

print("\n--- Python gamma ---")

xc = x - np.mean(x)

for t in range(10):
    gamma = np.sum(
        xc[:N-t] * xc[t:]
    ) / (N-t)

    print(
        t,
        " gamma = ",
        f"{gamma:.15e}"
    )    

W = 28

xc = x - mean

gamma = np.array([
    np.sum(xc[:N-t] * xc[t:]) / (N - t)
    for t in range(W + 1)
])

dbias = gamma[0] + 2.0 * np.sum(gamma[1:W])

gamma_corrected = gamma + dbias / N

print("\n--- bias correction ---")
print("W     =", W)
print("dbias =", dbias)
print("dbias/N =", dbias / N)

print("\n--- corrected gamma ---")

for t in range(10):
    print(
        t,
        f"raw={gamma[t]:.15e}",
        f"corrected={gamma_corrected[t]:.15e}",
    )    


data = {
    "Parameter": [
        r"$\langle P \rangle$",
        r"$\tau_{\mathrm{int}}$",
        r"$\Delta\tau$",
        r"$W$",
        r"$\tau_{\mathrm{exp}}$",
        r"$\gamma_0$",
        r"$\sigma$"
    ],
    "Wert": [
        mean,
        tau_int,
        result["delta_tau"],
        result["W"],
        result["tau_exp"],
        gamma0,
        sigma
    ]
}

df = pd.DataFrame(data)

print(df.to_latex(index=False, escape=False, float_format="%.10f"))    