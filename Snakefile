configfile: "analysis/config/config.yaml"

import glob
from pathlib import Path

LOG_FILES = glob.glob(f"{config['log_dir']}/*.log")
LOG_NAMES = [Path(f).stem for f in LOG_FILES]

rule all:
    input:
        "results/combined.csv"

rule parse_log:
    input:
        log = f"{config['log_dir']}/{{logname}}.log"
    output:
        csv = "results/parsed/{logname}.csv"
    script:
        "analysis/scripts/parse_log.py"

rule merge_logs:
    input:
        expand("results/parsed/{logname}.csv", logname=LOG_NAMES)
    output:
        "results/combined.csv"
    script:
        "analysis/scripts/merge_logs.py"