import numpy as np
import pandas as pd


def select_scan(
    df: pd.DataFrame,
    *,
    beta: float | list[float] | None = None,
    mdsteps: int | list[int] | None = None,
    trajL: float | list[float] | None = None,
    trajectories: int | list[int] | None = None,
    exclude_beta: float | list[float] | None = None,
    exclude_mdsteps: int | list[int] | None = None,
    exclude_trajL: float | list[float] | None = None,
    exclude_trajectories: int | list[int] | None = None,
) -> pd.DataFrame:

    result = df.copy()

    filters = {
        "beta": beta,
        "mdsteps": mdsteps,
        "trajL": trajL,
        "trajectories": trajectories,
    }

    exclusions = {
        "beta": exclude_beta,
        "mdsteps": exclude_mdsteps,
        "trajL": exclude_trajL,
        "trajectories": exclude_trajectories,
    }

    for column, values in filters.items():
        if values is not None:
            if not isinstance(values, list):
                values = [values]

            result = result[result[column].isin(values)]

    for column, values in exclusions.items():
        if values is not None:
            if not isinstance(values, list):
                values = [values]

            result = result[~result[column].isin(values)]

    return result.copy()



def thermalized(
    df: pd.DataFrame,
    trajectory_cut: int,
) -> pd.DataFrame:

    return df.loc[
        df["trajectory"] > trajectory_cut
    ].copy()


def measurement_stats(
    df: pd.DataFrame,
    value: str,
    groupby: str | list[str],
) -> pd.DataFrame:
    """Calculate basic statistics for a measured quantity."""

    stats = (
        df.groupby(groupby)[value]
        .agg(
            N="count",
            mean="mean",
            std="std",
        )
        .reset_index()
    )

    stats["sem"] = stats["std"] / np.sqrt(stats["N"])

    return stats


def acceptance_stats(
    df: pd.DataFrame,
    groupby: str | list[str] = "beta",
) -> pd.DataFrame:

    return (
        df.groupby(groupby)["accepted"]
        .agg(
            N="count",
            accepted="sum",
            acceptance="mean",
        )
        .reset_index()
    )

def autocorr(
    df: pd.DataFrame,
    value: str,
    max_lag: int,
) -> pd.DataFrame:

    x = df[value].dropna().to_numpy()
    N = len(x)

    if N < 2:
        raise ValueError("Not enough data for autocorrelation.")

    if max_lag > N:
        raise ValueError(
            f"max_lag={max_lag} is larger than N={N}."
        )

    mean = np.mean(x)

    autocorrs = []

    for lag in range(max_lag):
        x1 = x[:N - lag]
        x2 = x[lag:]

        gamma = np.mean(
            (x1 - mean) * (x2 - mean)
        )

        autocorrs.append(gamma)

    autocorrs = np.asarray(autocorrs)

    if autocorrs[0] == 0:
        raise ValueError(
            f"Variance of '{value}' is zero."
        )

    autocorrs /= autocorrs[0]

    return pd.DataFrame({
        "lag": np.arange(max_lag),
        "autocorr": autocorrs,
    })


def jackknife_autocorr(df, value, t, n_blocks): #autocorr with jackknife error estimation

    N = len(df)

    if n_blocks < 2:
        raise ValueError("n_blocks must be at least 2.")

    if n_blocks > N:
        raise ValueError(
            f"n_blocks={n_blocks} cannot exceed N={N}."
        )
    block_size = N // n_blocks

    # Berechne die Autokorrelation für das Original-Dataset
    rho_original = autocorr(
        df,
        value,
        t
    )
        
    jackknife_rhos = []

    for i in range(n_blocks):

        start = i * block_size
        end = (
            (i + 1) * block_size
            if i < n_blocks - 1
            else N
        )

        # Daten ohne Block i
        df_jackknife = pd.concat([
            df.iloc[:start],
            df.iloc[end:]
        ])

        # Autokorrelation für dieses Jackknife-Sample
        rho_i = autocorr(
            df_jackknife,
            value,
            t
        )

        jackknife_rhos.append(rho_i)


    rho_jackknife = jackknife_rhos[0].copy()

    # Mittelwert der Jackknife-Samples
    rho_jackknife["autocorr"] = (
        sum(
            rho["autocorr"]
            for rho in jackknife_rhos
        )
        / n_blocks
    )

    variance = np.zeros(t)

    for rho in jackknife_rhos:

        variance += (
            rho["autocorr"].values
            - rho_jackknife["autocorr"].values
        ) ** 2
    #Varianz der Jackknife-Samples
    variance *= (n_blocks - 1) / n_blocks

    delta_rho = np.sqrt(variance)


    result = rho_original.copy()
    result["delta_autocorr"] = delta_rho

    return result


def integrated_autocorr_time(
    df: pd.DataFrame,
    value: str,
    W: int,
) -> float:

    ac = autocorr(df, value, W + 1)

    return 0.5 + ac["autocorr"].iloc[1:].sum()


def integrated_autocorr_time_scan(df, value, W):

    N = len(df)

    ac_df = autocorr(df, value, W + 1)

    tau_ints = []
    delta_taus = []

    for w in range(1, W + 1):

        tau_int = 0.5 + np.sum(
            ac_df["autocorr"][1:w + 1]
        )

        # Wolff error formula
        argument = w + 0.5 - tau_int

        if argument > 0:
            delta_tau = (
                2 * tau_int
                * np.sqrt(argument / N)
            )
        else:
            delta_tau = np.nan

        tau_ints.append(tau_int)
        delta_taus.append(delta_tau)

    return pd.DataFrame({
        "W": range(1, W + 1),
        "tau_int": tau_ints,
        "delta_tau": delta_taus
    })

def tau_exp_from_tau_int(
    tau_int: np.ndarray,
    S: float = 1,
) -> np.ndarray:

    tau_int = np.asarray(tau_int, dtype=float)

    tau_exp = np.full_like(tau_int, 1e-10)

    mask = tau_int > 0.5

    tau_exp[mask] = (
        S
        / np.log(
            (2 * tau_int[mask] + 1)
            / (2 * tau_int[mask] - 1)
        )
    )

    return tau_exp


def window(
    tau_exp: np.ndarray,
    N: int,
) -> pd.DataFrame:

    W = np.arange(1, len(tau_exp))

    g = (
        np.exp(-W / tau_exp[1:])
        - tau_exp[1:] / (np.sqrt(W) * N)
    )

    return pd.DataFrame({
        "W": W,
        "g": g,
    })


def find_optimal_W(
    df: pd.DataFrame,
    value: str,
    Wmax: int,
    S: float = 1,
) -> dict:

    N = len(df)

    tau_df = integrated_autocorr_time_scan(
        df,
        value,
        Wmax,
    )

    tau_int = tau_df["tau_int"].to_numpy()
    delta_tau = tau_df["delta_tau"].to_numpy()

    tau_exp = tau_exp_from_tau_int(
        tau_int,
        S,
    )

    g = window(
        tau_exp,
        N,
    )

    negative = np.where(g < 0)[0]

    if len(negative) == 0:
        W_opt = Wmax
    else:
        W_opt = negative[0] + 1

    idx = W_opt - 1

    return {
        "W": W_opt,
        "tau_int": tau_int[idx],
        "delta_tau": delta_tau[idx],
        "tau_exp": tau_exp[idx],
    }

def autocorr_error(
    df: pd.DataFrame,
    value: str,
    groupby: str | list[str],
    Wmax: int,
    S: float = 1,
) -> pd.DataFrame:

    results = []

    for group, group_df in df.groupby(groupby):

        group_df = group_df.dropna(subset=[value])
        N = len(group_df)

        result = find_optimal_W(
            group_df,
            value,
            Wmax,
            S,
        )

        tau_int = result["tau_int"]

        gamma0 = np.var(
            group_df[value].to_numpy(),
            ddof=0,
        )

        sigma = np.sqrt(
            2 * tau_int * gamma0 / N
        )

        row = {
            "tau_int": tau_int,
            "delta_tau": result["delta_tau"],
            "W_opt": result["W"],
            "tau_exp": result["tau_exp"],
            "gamma0": gamma0,
            "sem_autocorr": sigma,
        }

        if isinstance(groupby, str):
            row[groupby] = group
        else:
            for name, value_ in zip(groupby, group):
                row[name] = value_

        results.append(row)

    return pd.DataFrame(results)










 

