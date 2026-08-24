import numpy as np
import pandas as pd


# ============================================================
# DATA SELECTION
# ============================================================

def select_scan(
    df: pd.DataFrame,
    *,
    beta: float | list[float] | None = None,
    mdsteps: int | list[int] | None = None,
    trajL: float | list[float] | None = None,
    trajectories: int | list[int] | None = None,
    lattice: str | list[str] | None = None,
    exclude_beta: float | list[float] | None = None,
    exclude_mdsteps: int | list[int] | None = None,
    exclude_trajL: float | list[float] | None = None,
    exclude_trajectories: int | list[int] | None = None,
    exclude_lattice: str | list[str] | None = None,
) -> pd.DataFrame:
    """
    Select a specific parameter scan from a DataFrame.

    Parameters can either be a single value or a list of values.

    Example
    -------
    df_beta = select_scan(
        df,
        beta=[2, 4, 6, 8, 10],
        mdsteps=10,
        trajL=1.0,
        trajectories=1000,
        lattice="8x8x8x8",
    )
    """

    result = df.copy()

    filters = {
        "beta": beta,
        "mdsteps": mdsteps,
        "trajL": trajL,
        "trajectories": trajectories,
        "lattice": lattice,
    }

    exclusions = {
        "beta": exclude_beta,
        "mdsteps": exclude_mdsteps,
        "trajL": exclude_trajL,
        "trajectories": exclude_trajectories,
        "lattice": exclude_lattice,
    }

    # --------------------------------------------------------
    # Include filters
    # --------------------------------------------------------

    for column, values in filters.items():

        if values is not None:

            if not isinstance(values, list):
                values = [values]

            result = result[
                result[column].isin(values)
            ]

    # --------------------------------------------------------
    # Exclude filters
    # --------------------------------------------------------

    for column, values in exclusions.items():

        if values is not None:

            if not isinstance(values, list):
                values = [values]

            result = result[
                ~result[column].isin(values)
            ]

    return result.copy()


# ============================================================
# THERMALIZATION
# ============================================================

def thermalized(
    df: pd.DataFrame,
    trajectory_cut: int,
) -> pd.DataFrame:
    """
    Remove the initial thermalization trajectories.

    Keeps trajectories with

        trajectory > trajectory_cut
    """

    return df.loc[
        df["trajectory"] > trajectory_cut
    ].copy()


# ============================================================
# BASIC MEASUREMENT STATISTICS
# ============================================================

def measurement_stats(
    df: pd.DataFrame,
    value: str,
    groupby: str | list[str],
) -> pd.DataFrame:
    """
    Calculate basic statistics for a measured quantity.

    Returns
    -------
    DataFrame containing

        N
        mean
        std
        sem_naive
    """

    stats = (
        df.groupby(groupby)[value]
        .agg(
            N="count",
            mean="mean",
            std="std",
        )
        .reset_index()
    )

    stats["sem_naive"] = (
        stats["std"]
        / np.sqrt(stats["N"])
    )

    return stats


# ============================================================
# ACCEPTANCE STATISTICS
# ============================================================

def acceptance_stats(
    df: pd.DataFrame,
    groupby: str | list[str] = "beta",
) -> pd.DataFrame:
    """
    Calculate acceptance statistics.

    Assumes that the column 'accepted' contains
    boolean values or 0/1 values.
    """

    return (
        df.groupby(groupby)["accepted"]
        .agg(
            N="count",
            accepted="sum",
            acceptance="mean",
        )
        .reset_index()
    )


# ============================================================
# AUTOCORRELATION
# ============================================================

def autocorr(
    df: pd.DataFrame,
    value: str,
    max_lag: int,
) -> pd.DataFrame:
    """
    Calculate the normalized autocorrelation function

        rho(t) = Gamma(t) / Gamma(0)

    for lags

        t = 0, ..., max_lag - 1.

    The DataFrame should already correspond to one
    statistically independent ensemble definition,
    e.g. fixed beta, mdsteps, trajL and lattice.
    """

    x = (
        df[value]
        .dropna()
        .to_numpy()
    )

    N = len(x)

    if N < 2:
        raise ValueError(
            f"Not enough data for autocorrelation: N={N}."
        )

    if max_lag < 1:
        raise ValueError(
            "max_lag must be at least 1."
        )

    if max_lag > N:
        raise ValueError(
            f"max_lag={max_lag} cannot exceed N={N}."
        )

    mean = np.mean(x)

    gamma = []

    for lag in range(max_lag):

        x1 = x[:N - lag]
        x2 = x[lag:]

        gamma_lag = np.mean(
            (x1 - mean)
            * (x2 - mean)
        )

        gamma.append(gamma_lag)

    gamma = np.asarray(gamma)

    if gamma[0] == 0:
        raise ValueError(
            f"Variance of '{value}' is zero."
        )

    rho = gamma / gamma[0]

    return pd.DataFrame({
        "lag": np.arange(max_lag),
        "autocorr": rho,
    })


# ============================================================
# JACKKNIFE AUTOCORRELATION
# ============================================================

def jackknife_autocorr(
    df: pd.DataFrame,
    value: str,
    max_lag: int,
    n_blocks: int,
) -> pd.DataFrame:
    """
    Calculate the autocorrelation function and its
    Jackknife uncertainty.

    The DataFrame must already contain ONE ensemble.

    Example
    -------
    beta=6, mdsteps=10, trajL=1.0, lattice=8x8x8x8

    Parameters
    ----------
    df : pandas.DataFrame
        Measurement data.
    value : str
        Observable.
    max_lag : int
        Maximum autocorrelation lag.
    n_blocks : int
        Number of Jackknife blocks.

    Returns
    -------
    DataFrame
        lag
        autocorr
        delta_autocorr
    """

    # --------------------------------------------------------
    # Prepare data
    # --------------------------------------------------------

    df = (
        df.dropna(subset=[value])
        .sort_values("trajectory")
        .reset_index(drop=True)
    )

    N = len(df)

    if N < 2:
        raise ValueError(
            f"Not enough data for Jackknife: N={N}."
        )

    if n_blocks < 2:
        raise ValueError(
            "n_blocks must be at least 2."
        )

    if n_blocks > N:
        raise ValueError(
            f"n_blocks={n_blocks} cannot exceed N={N}."
        )

    if N % n_blocks != 0:
        raise ValueError(
            f"N={N} must be divisible by "
            f"n_blocks={n_blocks}."
        )

    if max_lag > N:
        raise ValueError(
            f"max_lag={max_lag} cannot exceed N={N}."
        )

    # --------------------------------------------------------
    # Original autocorrelation
    # --------------------------------------------------------

    rho_original = autocorr(
        df,
        value,
        max_lag,
    )

    # --------------------------------------------------------
    # Jackknife blocks
    # --------------------------------------------------------

    block_size = N // n_blocks

    jackknife_rhos = []

    for i in range(n_blocks):

        start = i * block_size
        end = (i + 1) * block_size

        df_jackknife = pd.concat(
            [
                df.iloc[:start],
                df.iloc[end:],
            ]
        ).reset_index(drop=True)

        rho_i = autocorr(
            df_jackknife,
            value,
            max_lag,
        )

        jackknife_rhos.append(rho_i)

    # --------------------------------------------------------
    # Jackknife mean
    # --------------------------------------------------------

    rho_jackknife = np.mean(
        [
            rho["autocorr"].to_numpy()
            for rho in jackknife_rhos
        ],
        axis=0,
    )

    # --------------------------------------------------------
    # Jackknife variance
    # --------------------------------------------------------

    variance = (
        (n_blocks - 1)
        / n_blocks
        * np.sum(
            [
                (
                    rho["autocorr"].to_numpy()
                    - rho_jackknife
                ) ** 2
                for rho in jackknife_rhos
            ],
            axis=0,
        )
    )

    delta_rho = np.sqrt(variance)

    # --------------------------------------------------------
    # Result
    # --------------------------------------------------------

    result = rho_original.copy()

    result["delta_autocorr"] = delta_rho

    return result


# ============================================================
# JACKKNIFE AUTOCORRELATION SCAN
# ============================================================

def jackknife_autocorr_scan(
    df: pd.DataFrame,
    value: str,
    groupby: str | list[str],
    max_lag: int,
    n_blocks: int,
) -> pd.DataFrame:
    """
    Calculate autocorrelation with Jackknife errors
    for several ensembles.

    Example
    -------
    jackknife_autocorr_scan(
        df,
        value="plaquette",
        groupby="beta",
        max_lag=30,
        n_blocks=10,
    )

    For a completely unambiguous scan one can use

        groupby=[
            "beta",
            "mdsteps",
            "trajL",
            "lattice",
        ]
    """

    results = []

    for group, group_df in df.groupby(groupby):

        rho = jackknife_autocorr(
            group_df,
            value=value,
            max_lag=max_lag,
            n_blocks=n_blocks,
        )

        if isinstance(groupby, str):

            rho[groupby] = group

        else:

            for name, value_ in zip(
                groupby,
                group,
            ):
                rho[name] = value_

        results.append(rho)

    if not results:
        raise ValueError(
            "No groups found in DataFrame."
        )

    return pd.concat(
        results,
        ignore_index=True,
    )


# ============================================================
# INTEGRATED AUTOCORRELATION TIME
# ============================================================

def integrated_autocorr_time(
    df: pd.DataFrame,
    value: str,
    Wmax: int,
) -> pd.DataFrame:
    """
    Calculate the integrated autocorrelation time

        tau_int(W)
            = 1/2 + sum_{t=1}^W rho(t)

    and the Wolff statistical uncertainty.

    The DataFrame must contain ONE ensemble.
    """

    df = (
        df.dropna(subset=[value])
        .sort_values("trajectory")
        .reset_index(drop=True)
    )

    N = len(df)

    if N < 2:
        raise ValueError(
            f"Not enough data: N={N}."
        )

    if Wmax < 1:
        raise ValueError(
            "Wmax must be at least 1."
        )

    if Wmax >= N:
        raise ValueError(
            f"Wmax={Wmax} must be smaller than N={N}."
        )

    ac_df = autocorr(
        df,
        value,
        Wmax + 1,
    )

    tau_ints = []
    delta_taus = []

    for W in range(1, Wmax + 1):

        tau_int = (
            0.5
            + ac_df["autocorr"]
            .iloc[1:W + 1]
            .sum()
        )

        # ----------------------------------------------------
        # Wolff error estimate
        # ----------------------------------------------------

        argument = (
            W
            + 0.5
            - tau_int
        )

        if argument > 0:

            delta_tau = (
                2
                * tau_int
                * np.sqrt(
                    argument / N
                )
            )

        else:

            delta_tau = np.nan

        tau_ints.append(tau_int)
        delta_taus.append(delta_tau)

    return pd.DataFrame({
        "W": np.arange(
            1,
            Wmax + 1,
        ),
        "tau_int": tau_ints,
        "delta_tau": delta_taus,
    })


# ============================================================
# INTEGRATED AUTOCORRELATION SCAN
# ============================================================

def integrated_autocorr_time_scan(
    df: pd.DataFrame,
    value: str,
    groupby: str | list[str],
    Wmax: int,
) -> pd.DataFrame:
    """
    Calculate tau_int(W) for several ensembles.

    Returns columns

        beta / other grouping variables
        W
        tau_int
        delta_tau
    """

    results = []

    for group, group_df in df.groupby(groupby):

        tau_df = integrated_autocorr_time(
            group_df,
            value=value,
            Wmax=Wmax,
        )

        if isinstance(groupby, str):

            tau_df[groupby] = group

        else:

            for name, value_ in zip(
                groupby,
                group,
            ):
                tau_df[name] = value_

        results.append(tau_df)

    if not results:
        raise ValueError(
            "No groups found in DataFrame."
        )

    return pd.concat(
        results,
        ignore_index=True,
    )


# ============================================================
# EXPONENTIAL AUTOCORRELATION TIME
# ============================================================

def tau_exp_from_tau_int(
    tau_int: np.ndarray,
    S: float = 1,
) -> np.ndarray:
    """
    Convert integrated autocorrelation time to
    exponential autocorrelation time.

        tau_exp =
            S / log((2 tau_int + 1)
                    /(2 tau_int - 1))

    Returns np.nan where tau_int <= 0.5 (formula undefined).
    """

    tau_int = np.asarray(tau_int, dtype=float)

    tau_exp = np.full_like(tau_int, 1e-10)   # Ist das so sinnvoll? Muss ich ein besseres Kriiterium um ungultige Werte zu erkennen finden?

    mask = tau_int > 0.5

    tau_exp[mask] = (
        S
        / np.log(
            (2 * tau_int[mask] + 1)
            / (2 * tau_int[mask] - 1)
        )
    )

    return tau_exp


# ============================================================
# WINDOW FUNCTION
# ============================================================

def window_function(
    W: np.ndarray,
    tau_exp: np.ndarray,
    N: int,
) -> np.ndarray:
    """
    Wolff window function

        g(W) =
            exp(-W / tau_exp)
            - tau_exp / (sqrt(W) * N)

    Returns np.nan for invalid entries.
    """

    W = np.asarray(W, dtype=float)
    tau_exp = np.asarray(tau_exp, dtype=float)

    g = np.full_like(W, np.nan, dtype=float)

    valid = (
        (W > 0)                    # fix: parentheses restore correct precedence
        & np.isfinite(tau_exp)
        & (tau_exp > 0)
    )

    g[valid] = (
        np.exp(-W[valid] / tau_exp[valid])
        - tau_exp[valid] / (np.sqrt(W[valid]) * N)
    )

    return g


# ============================================================
# OPTIMAL WINDOW
# ============================================================

def find_optimal_W(
    df: pd.DataFrame,
    value: str,
    Wmax: int,
    S: float = 1,
) -> dict:
    """
    Find the optimal window W_opt via the Wolff criterion:
    the first W where g(W) < 0.

    Falls back to Wmax if no such W is found.
    """

    # fix: compute N after dropna, consistent with integrated_autocorr_time
    N = len(df[value].dropna())

    tau_df = integrated_autocorr_time(
        df,
        value,
        Wmax,
    )

    tau_int = tau_df["tau_int"].to_numpy()
    delta_tau = tau_df["delta_tau"].to_numpy()

    tau_exp = tau_exp_from_tau_int(tau_int, S)

    W_arr = np.arange(1, Wmax + 1)                  # fix: construct W array
    g = window_function(W_arr, tau_exp, N)           # fix: correct call

    negative = np.where(g < 0)[0]

    if len(negative) == 0:
        idx = Wmax - 1
    else:
        idx = int(negative[0])

    W_opt = idx + 1

    return {
        "W": W_opt,
        "tau_int": tau_int[idx],
        "delta_tau": delta_tau[idx],
        "tau_exp": tau_exp[idx],
    }


# ============================================================
# FINAL AUTOCORRELATION-CORRECTED STATISTICS
# ============================================================

def autocorr_error(
    df: pd.DataFrame,
    value: str,
    groupby: str | list[str],
    Wmax: int,
    S: float = 1,
) -> pd.DataFrame:
    """
    Calculate autocorrelation-corrected errors
    for several ensembles.

    The final error is

        sigma_auto =
            sqrt(2 * tau_int * Gamma(0) / N)

    where Gamma(0) = Var(X).
    """

    results = []

    for group, group_df in df.groupby(groupby):

        group_df = (
            group_df
            .dropna(subset=[value])
            .sort_values("trajectory")
            .reset_index(drop=True)
        )

        N = len(group_df)

        if N < 2:
            raise ValueError(
                f"Not enough data in group {group}: N={N}."
            )

        result = find_optimal_W(
            group_df,
            value,
            Wmax,
            S,
        )

        tau_int = result["tau_int"]
        x = group_df[value].to_numpy()

        gamma0 = np.var(x, ddof=0)

        sem_autocorr = np.sqrt(
            2 * tau_int * gamma0 / N
        )

        row = {
            "mean": np.mean(x),          # added: callers almost always need this
            "tau_int": tau_int,
            "delta_tau": result["delta_tau"],
            "W_opt": result["W"],
            "tau_exp": result["tau_exp"],
            "gamma0": gamma0,
            "sem_autocorr": sem_autocorr,
        }

        if isinstance(groupby, str):
            row[groupby] = group
        else:
            for name, value_ in zip(groupby, group):
                row[name] = value_

        results.append(row)

    return pd.DataFrame(results)
