import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# ============================================================
# Allgemeine Hilfsfunktion
# ============================================================

def savefig(fig, filename, dpi=300):
    """
    Speichert eine Figure und schließt sie anschließend.
    """
    fig.savefig(
        filename,
        dpi=dpi,
        bbox_inches="tight"
    )
    plt.close(fig)


# ============================================================
# PLAQUETTE
# ============================================================

def plot_plaquette_history(df, title=None, ax=None):
    """
    Plaquette als Funktion der Trajektorie für alle beta-Werte.
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained"
        )

    for beta in sorted(df["beta"].unique()):

        subset = (
            df[df["beta"] == beta]
            .sort_values("trajectory")
        )

        ax.plot(
            subset["trajectory"],
            subset["plaquette"],
            ".",
            label=fr"$\beta={beta}$"
        )

    ax.set_xlabel("trajectory")
    ax.set_ylabel("plaquette")
    ax.grid()
    ax.legend()

    if title:
        ax.set_title(title)

    return ax


def plot_plaquette_histogram(
    df,
    beta,
    bins=30,
    ax=None
):
    """
    Histogramm der Plaquette-Verteilung für einen beta-Wert.
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained"
        )

    values = df["plaquette"]

    mean = values.mean()
    std = values.std()

    ax.hist(
        values,
        bins=bins,
        density=True,
        alpha=0.7,
        rwidth=0.5
    )

    ax.axvline(
        mean,
        linestyle="--",
        label=fr"$\langle P\rangle={mean:.5f}$"
    )

    ax.set_xlabel(r"$P$")
    ax.set_ylabel("Probability density")

    ax.set_title(
        fr"Plaquette distribution $\beta={beta}$"
    )

    ax.set_xlim(
        mean - 5 * std,
        mean + 5 * std
    )

    ax.grid()
    ax.legend()

    return ax


def plot_plaquette_average(
    df,
    trajectories=None,
    ax=None
):
    """
    Mittelwert der Plaquette als Funktion von beta.

    df benötigt:
        beta
        mean
        sem
        sem_autocorr
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained"
        )

    dx = 0.01

    # Naiver Fehler
    ax.errorbar(
        df["beta"] - dx,
        df["mean"],
        yerr=df["sem_naive"],
        fmt="o",
        capsize=3,
        label="Naiver Fehler"
    )

    # Autokorrelations-korrigierter Fehler
    ax.errorbar(
        df["beta"] + dx,
        df["mean"],
        yerr=df["sem_autocorr"],
        fmt="s",
        capsize=3,
        label="Autokorrelations-korrigiert"
    )

    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\langle P\rangle$")

    if trajectories is not None:
        ax.set_title(
            rf"$\langle P\rangle$ vs $\beta$: "
            rf"{trajectories} trajectories"
        )
    else:
        ax.set_title(r"$\langle P\rangle$ vs $\beta$")

    ax.grid()
    ax.legend()

    return ax


# ============================================================
# AUTOKORRELATION1
# ============================================================

def plot_autocorr(
    df,
    ax=None,
    title=None,
    label=None,
    color=None
):
    """
    Autokorrelationsfunktion rho(t).
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained"
        )

    ax.errorbar(
        df["lag"],
        df["autocorr"],
        yerr=df["delta_autocorr"],
        fmt="o",
        markersize=4,
        color=color,
        ecolor=color,
        elinewidth=1.2,
        capsize=3,
        capthick=1.2,
        linestyle="none",
        label=label
    )

    ax.axhline(
        0,
        linestyle="--",
        linewidth=1,
        alpha=0.5
    )

    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\rho(t)$")

    ax.grid(True, alpha=0.3)

    if title:
        ax.set_title(title)

    if label:
        ax.legend()

    return ax


def plot_autocorr_scan(
    df,
    betas=None,
    max_lag=None,
    title=None,
    ax=None,
):
    """
    Plot der Autokorrelationsfunktion für mehrere beta-Werte.

    Erwartete Spalten:
        beta
        lag
        autocorr
        delta_autocorr
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained",
        )

    if betas is None:
        betas = sorted(df["beta"].unique())

    for beta in betas:

        subset = df[
            df["beta"] == beta
        ].sort_values("lag")

        if max_lag is not None:
            subset = subset[
                subset["lag"] <= max_lag
            ]

        plot_autocorr(
            subset,
            ax=ax,
            label=fr"$\beta={beta}$",
        )

    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\rho(t)$")
    ax.set_yscale("log")

    if title is not None:
        ax.set_title(title)

    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax


def plot_int_autocorr(
    df,
    ax=None,
    title=None,
    label=None,
    color=None,
):
    """
    Plot von tau_int(W) mit Wolff-Fehlerband.

    Erwartete Spalten:
        W
        tau_int
        delta_tau
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained",
        )

    df = df.sort_values("W")

    # tau_int(W)
    ax.plot(
        df["W"],
        df["tau_int"],
        ".-",
        color=color,
        label=label,
        linewidth=1.5,
        markersize=4,
    )

    # Wolff error band
    valid = (
        df["delta_tau"].notna()
        & (df["delta_tau"] >= 0)
    )

    ax.fill_between(
        df.loc[valid, "W"],
        (
            df.loc[valid, "tau_int"]
            - df.loc[valid, "delta_tau"]
        ),
        (
            df.loc[valid, "tau_int"]
            + df.loc[valid, "delta_tau"]
        ),
        color=color,
        alpha=0.2,
    )

    ax.set_xlabel(r"$W$")
    ax.set_ylabel(
        r"$\tau_{\mathrm{int}}(W)$"
    )

    ax.grid(True, alpha=0.3)

    if title is not None:
        ax.set_title(title)

    return ax


def plot_int_autocorr_scan(
    df,
    betas=None,
    Wmax=None,
    title=None,
    ax=None,
):
    """
    Plot von tau_int(W) für mehrere beta-Werte.

    Erwartete Spalten:
        beta
        W
        tau_int
        delta_tau
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained",
        )

    if betas is None:
        betas = sorted(df["beta"].unique())

    # Farben aus der Matplotlib-Farbpalette
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for i, beta in enumerate(betas):

        subset = (
            df[df["beta"] == beta]
            .sort_values("W")
        )

        if Wmax is not None:
            subset = subset[
                subset["W"] <= Wmax
            ]

        color = colors[i % len(colors)]

        plot_int_autocorr(
            subset,
            ax=ax,
            label=fr"$\beta={beta}$",
            color=color,
        )

    ax.set_xlabel(r"$W$")
    ax.set_ylabel(
        r"$\tau_{\mathrm{int}}(W)$"
    )

    ax.grid(True, alpha=0.3)
    ax.legend()

    if title is not None:
        ax.set_title(title)

    return ax
# ============================================================
# ACCEPTANCE
# ============================================================

def plot_acceptance(
    df_stats,
    trajectories=None,
    ax=None
):
    """
    Acceptance Rate als Funktion von beta.
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(7, 5),
            layout="constrained"
        )

    ax.plot(
        df_stats["beta"],
        df_stats["acceptance"],
        "o"
    )

    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel("Acceptance rate")

    if trajectories is not None:
        ax.set_title(
            f"Acceptance rate "
            f"({trajectories} trajectories)"
        )

    ax.grid()

    return ax


def plot_acceptance_vs_dH(
    df,
    trajectories=None,
    compare_naive=False,
    ax=None
):
    """
    Acceptance Rate gegen <Delta H>.
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained"
        )

    # Autokorrelations-korrigierter Fehler
    ax.errorbar(
            df["mean"],
            df["acceptance"],
            xerr=df["sem_autocorr"],
            fmt="o",
            capsize=7,
            label="autocorrelation corrected"
        )

    # Naiver Fehler
    if compare_naive:
        ax.errorbar(
        df["mean"],
        df["acceptance"],
        xerr=df["sem"],
        fmt="o",
        capsize=4,
        label="naive"
    )

    # beta beschriften
    for _, row in df.iterrows():

        ax.annotate(
            fr"$\beta={row['beta']}$",
            (
                row["mean"],
                row["acceptance"]
            ),
            xytext=(10, 10),
            textcoords="offset points"
        )

    ax.set_xlabel(r"$\langle \Delta H\rangle$")
    ax.set_ylabel("Acceptance rate")

    if trajectories is not None:
        ax.set_title(
            r"Acceptance rate vs $\langle \Delta H\rangle$ "
            f"({trajectories} trajectories)"
        )
    else:
        ax.set_title(
            r"Acceptance rate vs $\langle \Delta H\rangle$"
        )

    ax.margins(x=0.15, y=0.15)
    ax.grid()

    if compare_naive:
        ax.legend()

    return ax


# ============================================================
# DELTA H
# ============================================================

def plot_dH_history(
    df,
    title=None,
    ax=None
):
    """
    Delta H als Funktion der Trajektorie für alle beta.
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained"
        )

    for beta in sorted(df["beta"].unique()):

        subset = (
            df[df["beta"] == beta]
            .sort_values("trajectory")
        )

        ax.plot(
            subset["trajectory"],
            subset["dH"],
            ".",
            alpha=0.5,
            label=fr"$\beta={beta}$"
        )

    ax.set_xlabel("trajectory")
    ax.set_ylabel(r"$\Delta H$")

    ax.axhline(
        0,
        linestyle="--"
    )

    ax.grid()
    ax.legend()

    if title:
        ax.set_title(title)

    return ax


def plot_dH_histogram(
    df,
    beta=None,
    bins=40,
    trajectories=None,
    ax=None
):
    """
    Histogramm der Delta-H-Verteilung.

    Wenn beta angegeben wird:
        nur dieser beta-Wert.

    Wenn beta=None:
        alle beta-Werte.
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained"
        )

    if beta is not None:

        df_plot = df[df["beta"] == beta]

        ax.hist(
            df_plot["dH"],
            bins=bins,
            density=True,
            rwidth=0.9,
            alpha=0.7
        )

        title = (
            fr"$\Delta H$ distribution, "
            fr"$\beta={beta}$"
        )

    else:

        for b in sorted(df["beta"].unique()):

            ax.hist(
                df[df["beta"] == b]["dH"],
                bins=bins,
                density=True,
                alpha=0.5,
                label=fr"$\beta={b}$"
            )

        ax.legend()

        title = r"$\Delta H$ distribution"

    if trajectories is not None:
        title += f" ({trajectories} trajectories)"

    ax.axvline(
        0,
        linestyle="--"
    )

    ax.set_xlabel(r"$\Delta H$")
    ax.set_ylabel("Probability density")
    ax.set_title(title)

    ax.grid()

    return ax


def plot_exp_dH(
    df_stats,
    trajectories=None,
    ax=None
):
    """
    <exp(-Delta H)> als Funktion von beta.
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained"
        )

    ax.errorbar(
        df_stats["beta"],
        df_stats["mean"],
        yerr=df_stats["sem_autocorr"],
        fmt="o",
        capsize=4
    )

    ax.axhline(
        1,
        linestyle="--",
        label=r"$\langle e^{-\Delta H}\rangle=1$"
    )

    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\langle e^{-\Delta H}\rangle$")

    if trajectories is not None:
        ax.set_title(
            r"$\langle e^{-\Delta H}\rangle$ "
            f"({trajectories} trajectories)"
        )

    ax.grid()
    ax.legend()

    return ax


# ============================================================
# ALLGEMEINE AVERAGE-PLOT-FUNKTION
# ============================================================

def plot_average_with_error(
    df,
    x,
    x_label,
    y_label,
    title=None,
    compare_naive=False,
    ax=None
):
    """
    Allgemeiner Plot für Mittelwerte mit Fehlerbalken.
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(7, 5),
            layout="constrained"
        )
        # Autokorrelations-korrigierter Fehler
        ax.errorbar(
            df[x],
            df["mean"],
            yerr=df["sem_autocorr"],
            fmt="o",
            capsize=7,
            label="autocorrelation corrected"
        )
    


    # Naiver Fehler
    if compare_naive:
        ax.errorbar(
        df[x],
        df["mean"],
        yerr=df["sem"],
        fmt="o",
        capsize=4,
        label="naive"
    )


        ax.legend()

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if title:
        ax.set_title(title)

    ax.grid()

    return ax


# ============================================================
# INTEGRATOR
# ============================================================

def plot_integrator_comparison(
    df_mdsteps,
    df_trajL,
    observable,
    trajectories=None,
    reference=None,
    compare_autocorr=False,
    ax=None
):
    """
    Vergleich zwischen

        varying N_MD

    und

        varying L_traj

    als Funktion von Delta tau.
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained"
        )

    color_md = "tab:blue"
    color_traj = "tab:orange"

    # --------------------------------------------------------
    # N_MD scan
    # --------------------------------------------------------

    ax.errorbar(
        df_mdsteps["dtau"],
        df_mdsteps["mean"],
        yerr=df_mdsteps["sem"],
        fmt="o",
        color=color_md,
        ecolor=color_md,
        capsize=3,
        capthick=1,
        elinewidth=1,
        alpha=0.8,
        label=r"varying $N_{\mathrm{MD}}$"
    )

    if compare_autocorr:

        ax.errorbar(
            df_mdsteps["dtau"],
            df_mdsteps["mean"],
            yerr=df_mdsteps["sem_autocorr"],
            fmt="none",
            ecolor=color_md,
            capsize=8,
            capthick=2,
            elinewidth=3,
            alpha=0.35
        )

    # --------------------------------------------------------
    # L_traj scan
    # --------------------------------------------------------

    ax.errorbar(
        df_trajL["dtau"],
        df_trajL["mean"],
        yerr=df_trajL["sem"],
        fmt="s",
        color=color_traj,
        ecolor=color_traj,
        capsize=3,
        capthick=1,
        elinewidth=1,
        alpha=0.8,
        label=r"varying $L_{\mathrm{traj}}$"
    )

    if compare_autocorr:

        ax.errorbar(
            df_trajL["dtau"],
            df_trajL["mean"],
            yerr=df_trajL["sem_autocorr"],
            fmt="none",
            ecolor=color_traj,
            capsize=8,
            capthick=2,
            elinewidth=3,
            alpha=0.35
        )

    # --------------------------------------------------------
    # Referenzwert
    # --------------------------------------------------------

    if reference is not None:

        ax.axhline(
            reference,
            color="black",
            linestyle=":",
            linewidth=1
        )

    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------

    handles, labels = ax.get_legend_handles_labels()

    if compare_autocorr:

        handles.extend([
            Line2D(
                [0],
                [0],
                color="black",
                linewidth=1,
                label="naive error"
            ),
            Line2D(
                [0],
                [0],
                color="black",
                linewidth=3,
                alpha=0.35,
                label="autocorrelation corrected"
            )
        ])

    ax.legend(handles=handles)

    # --------------------------------------------------------
    # Labels
    # --------------------------------------------------------

    ax.set_xlabel(
        r"$\Delta\tau=L_{\rm traj}/N_{\rm MD}$"
    )

    ax.set_ylabel(observable)

    if trajectories is not None:

        ax.set_title(
            f"{observable} vs $\\Delta\\tau$ "
            f"({trajectories} trajectories)"
        )

    else:

        ax.set_title(
            f"{observable} vs $\\Delta\\tau$"
        )

    ax.grid()

    return ax


# ============================================================
# TOPOLOGIE
# ============================================================

def plot_topology_history(
    df,
    beta=None,
    trajectories=None,
    ax=None
):
    """
    Topologische Ladung Q als Funktion der Trajektorie.

    beta=None:
        alle beta-Werte

    beta=<Wert>:
        nur ein beta-Wert
    """

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(8, 5),
            layout="constrained"
        )

    if beta is None:

        for b in sorted(df["beta"].unique()):

            subset = (
                df[df["beta"] == b]
                .sort_values("trajectory")
            )

            ax.plot(
                subset["trajectory"],
                subset["Q"],
                ".-",
                label=fr"$\beta={b}$"
            )

        title = "Topological charge"

        ax.legend()

    else:

        subset = (
            df[df["beta"] == beta]
            .sort_values("trajectory")
        )

        ax.plot(
            subset["trajectory"],
            subset["Q"],
            ".-"
        )

        title = fr"Topological charge ($\beta={beta}$)"

    if trajectories is not None:
        title += f" ({trajectories} trajectories)"

    ax.axhline(
        0,
        color="black",
        linestyle="--",
        linewidth=1
    )

    ax.set_xlabel("Trajectory")
    ax.set_ylabel(r"$Q$")
    ax.set_title(title)

    ax.grid()

    return ax

