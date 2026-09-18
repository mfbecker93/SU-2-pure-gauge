configfile: "analysis/config/config.yaml"

import glob
from pathlib import Path

LOG_FILES = glob.glob(f"{config['log_dir']}/*.log")
LOG_NAMES = [Path(f).stem for f in LOG_FILES]

rule all:
    input:
        "results/parsed_combined/combined.csv",
        "results/scan/scan.csv",
        "results/scan/plaquette_history.png",
        "results/scan/scan_thermalised.csv",
        "results/scan/plaquette_autocorr_jackknife.csv",
        "results/scan/plaquette_autocorr_jackknife.png",
        "results/scan/integrated_autocorr_time_scan.csv",
        "results/scan/integrated_autocorr_time_scan.png"

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
        "results/parsed_combined/combined.csv"
    script:
        "analysis/scripts/merge_logs.py"

rule select_scan:
    input:
        log = "results/parsed_combined/combined.csv"
    output:
        csv = "results/scan/scan.csv"
    script:
        "analysis/scripts/plaquette_beta_full.py"

rule plot_plaquette_history:
    input:
        csv = "results/scan/scan.csv"
    output:
        png = "results/scan/plaquette_history.png"
    script:
        "analysis/scripts/plot_plaquette_beta_full.py"  

rule thermalised:
    input:
        csv = "results/scan/scan.csv"
    output:
        csv = "results/scan/scan_thermalised.csv"
    script:
        "analysis/scripts/thermalised.py"

rule plaquette_autocorr_jackknife:
    input:
        csv = "results/scan/scan_thermalised.csv"      
    output:
        csv = "results/scan/plaquette_autocorr_jackknife.csv"
    script:
        "analysis/scripts/plaquette_autocorr_jackknife.py"  

rule plot_plaquette_autocorr_jackknife:
    input:
        csv = "results/scan/plaquette_autocorr_jackknife.csv"
    output:
        png = "results/scan/plaquette_autocorr_jackknife.png"
    script:
        "analysis/scripts/plot_plaquette_autocorr_jackknife.py"

rule integrated_autocorr_time_scan:
    input:
        csv = "results/scan/scan_thermalised.csv" 
    output:
        csv = "results/scan/integrated_autocorr_time_scan.csv"
    script:
        "analysis/scripts/integrated_autocorr_time_scan.py"

rule plot_integrated_autocorr_time_scan:
    input:
        csv = "results/scan/integrated_autocorr_time_scan.csv"
    output:
        png = "results/scan/integrated_autocorr_time_scan.png"
    script:
        "analysis/scripts/plot_integrated_autocorr_time_scan.py"