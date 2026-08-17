from pathlib import Path
import re
import numpy as np
import pandas as pd


# ============================================================
# File name
# ============================================================

def parse_log_name(path):
    """
    Extract simulation parameters from a log file name.
    """

    name = Path(path).name

    result = {}

    m = re.search(r"L(\d+x\d+x\d+x\d+)", name)
    if m:
        result["lattice"] = m.group(1)

    m = re.search(r"_beta([\d.]+)", name)
    if m:
        result["beta"] = float(m.group(1))

    m = re.search(r"_md(\d+)", name)
    if m:
        result["mdsteps"] = int(m.group(1))

    m = re.search(r"_trajL([\d.]+)", name)
    if m:
        result["trajL"] = float(m.group(1))

    m = re.search(r"_traj(\d+)", name)
    if m:
        result["trajectories"] = int(m.group(1))

    m = re.search(r"_ckpt(\d+)", name)
    if m:
        result["checkpoint"] = int(m.group(1))

    return result


# ============================================================
# Plaquette
# ============================================================

def parse_plaquette(logfile):
    """
    Extract plaquette measurements from one log file.
    """

    data = []

    with open(logfile) as f:

        for line in f:

            m = re.search(
                r"Plaquette:\s+\[\s*(\d+)\s*\]\s+([\d\.Ee+-]+)",
                line
            )

            if m:
                data.append({
                    "trajectory": int(m.group(1)),
                    "plaquette": float(m.group(2))
                })

    df = pd.DataFrame(data)

    # Remove duplicate MPI output
    if not df.empty:
        df = df.drop_duplicates(
            subset=["trajectory", "plaquette"]
        )

    return df


def load_plaquette_data(log_dir):
    """
    Load plaquette measurements from all log files.
    """

    log_dir = Path(log_dir)
    logs = list(log_dir.glob("*.log"))

    measurements = []

    for log in logs:

        params = parse_log_name(log)
        df = parse_plaquette(log)

        for key, value in params.items():
            df[key] = value

        df["file"] = log.name

        measurements.append(df)

    if not measurements:
        return pd.DataFrame()

    return pd.concat(
        measurements,
        ignore_index=True
    )


# ============================================================
# Delta H
# ============================================================

def parse_dh(logfile):
    """
    Extract H_after and dH from one log file.
    """

    data = []

    current_trajectory = None

    with open(logfile) as f:

        for line in f:

            # Find current trajectory from plaquette output
            m_traj = re.search(
                r"Plaquette:\s+\[\s*(\d+)\s*\]",
                line
            )

            if m_traj:
                current_trajectory = int(m_traj.group(1))

            # Find H_after and dH
            m_dh = re.search(
                r"Total H after trajectory\s*=\s*"
                r"([\d\.Ee+-]+)\s+"
                r"dH\s*=\s*([\d\.Ee+-]+)",
                line
            )

            if m_dh and current_trajectory is not None:

                data.append({
                    "trajectory": current_trajectory,
                    "H_after": float(m_dh.group(1)),
                    "dH": float(m_dh.group(2))
                })

    df = pd.DataFrame(data)

    return df


def load_dh_data(log_dir):
    """
    Load Delta-H measurements from all log files.
    """

    log_dir = Path(log_dir)
    logs = list(log_dir.glob("*.log"))

    measurements = []

    for log in logs:

        params = parse_log_name(log)
        df = parse_dh(log)

        for key, value in params.items():
            df[key] = value

        df["file"] = log.name

        measurements.append(df)

    if not measurements:
        return pd.DataFrame()

    df = pd.concat(
        measurements,
        ignore_index=True
    )

    # Remove duplicate output
    df = df.drop_duplicates(
        subset=["file", "trajectory"]
    )

    # Useful derived observable
    df["exp_minus_dH"] = np.exp(-df["dH"])

    return df


# ============================================================
# Topological charge
# ============================================================

def parse_topology(logfile):
    """
    Extract topological charge measurements.
    """

    data = []

    with open(logfile) as f:

        for line in f:

            m = re.search(
                r"Topological Charge:\s*"
                r"\[\s*(\d+)\s*\]\s*"
                r"([-0-9.eE+]+)",
                line
            )

            if m:
                data.append({
                    "trajectory": int(m.group(1)),
                    "Q": float(m.group(2))
                })

    return pd.DataFrame(data)


def load_topology_data(log_dir):
    """
    Load topological charge measurements from all log files.
    """

    log_dir = Path(log_dir)
    logs = list(log_dir.glob("*.log"))

    measurements = []

    for log in logs:

        params = parse_log_name(log)
        df = parse_topology(logfile=log)

        for key, value in params.items():
            df[key] = value

        df["file"] = log.name

        measurements.append(df)

    if not measurements:
        return pd.DataFrame()

    return pd.concat(
        measurements,
        ignore_index=True
    )


# ============================================================
# Acceptance
# ============================================================

def parse_acceptance(logfile):

    data = []

    current_trajectory = None

    with open(logfile) as f:
        for line in f:

            # tatsächliche Trajektoriennummer aus dem Log
            m_traj = re.search(
                r"Plaquette:\s+\[\s*(\d+)\s*\]",
                line
            )

            if m_traj:
                current_trajectory = int(m_traj.group(1))

            # Acceptance der aktuellen Trajektorie
            if "Metropolis_test -- ACCEPTED" in line:
                data.append({
                    "trajectory": current_trajectory,
                    "accepted": 1
                })

            elif "Metropolis_test -- REJECTED" in line:
                data.append({
                    "trajectory": current_trajectory,
                    "accepted": 0
                })

    return pd.DataFrame(data)

def load_acceptance_data(log_dir):
    """
    Load Metropolis acceptance tests from all log files.
    """

    log_dir = Path(log_dir)
    logs = list(log_dir.glob("*.log"))

    measurements = []

    for log in logs:

        params = parse_log_name(log)
        df = parse_acceptance(log)

        for key, value in params.items():
            df[key] = value

        df["file"] = log.name

        measurements.append(df)

    if not measurements:
        return pd.DataFrame()

    return pd.concat(
        measurements,
        ignore_index=True
    )


def load_all_data(log_dir):

    df_plaq = load_plaquette_data(log_dir)
    df_dh = load_dh_data(log_dir)
    df_q = load_topology_data(log_dir)
    df_acc = load_acceptance_data(log_dir)

    keys = ["file", "trajectory"]

    # Nur Messgrößen aus den zusätzlichen DataFrames verwenden
    df_dh = df_dh[
        keys + ["H_after", "dH"]
    ]

    df_q = df_q[
        keys + ["Q"]
    ]

    df_acc = df_acc[
        keys + ["accepted"]
    ]

    # Plaquette als Basis
    df = df_plaq.merge(
        df_dh,
        on=keys,
        how="left",
        validate="one_to_one"
    )

    df = df.merge(
        df_q,
        on=keys,
        how="left",
        validate="one_to_one"
    )

    df = df.merge(
        df_acc,
        on=keys,
        how="left",
        validate="one_to_one"
    )

    # Metropolis weight
    df["exp_minus_dH"] = np.exp(-df["dH"])

    return df