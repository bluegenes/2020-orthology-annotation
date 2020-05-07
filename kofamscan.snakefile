"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s kofamscan.snakefile --configfile test_data/test_config.yml --cores 1 --use-conda  # use -n for dry run
Function: Download and run kofamscan 
"""

import os
import sys
import pandas as pd

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

# need samples csv, with columns: "sample", "peptide_fasta"
samplelist = config["samples_csv"]
samplesDF = pd.read_csv(samplelist, dtype=str, header=0, index_col = "sample")
SAMPLES=samplesDF.index.tolist()

out_dir = config.get("out_dir", "kofamscan_out")
database_dir = config.get("db_dir", "kofamscan_data")
logs_dir = os.path.join(out_dir, "logs")

rule all:
    input: 
        expand(os.path.join(out_dir, "{sample}.txt"), sample = SAMPLES)

# download ko list
rule download_ko_list:
    input: FTP.remote("ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz", static=True, keep_local=True, immediate_close=True)
    output: 
        gz_file = os.path.join(database_dir,"ko_list.gz"),
        unzipped = os.path.join(database_dir, "ko_list"),
    log: os.path.join(logs_dir,"download_ko_list.log")
    shell: 
        """
        mv {input} {output.gz_file}
        gunzip -c {output.gz_file} > {output.unzipped} 
        """

# download hmm profiles
rule download_ko_profiles:
    input: FTP.remote("ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz", static=True, keep_local=True, immediate_close=True)
    output: 
        gz_file = os.path.join(database_dir, "profiles.tar.gz"),
        unzipped = directory(os.path.join(database_dir,"profiles"))
    log: os.path.join(logs_dir,"download_profiles.log")
    params:
        database_dir = database_dir,
    shell: 
        """
        mv {input} {output.gz_file} 2> {log}
        tar xzf {output.gz_file} --directory {params.database_dir}
        """

# for options: https://github.com/takaram/kofam_scan
rule run_kofamscan:
    input:
        fasta = lambda w: samplesDF.loc[w.sample, "peptide_fasta"],
        profile_dir = os.path.join(database_dir,"profiles"),
        ko_list = os.path.join(database_dir, "ko_list") 
    output:
        os.path.join(out_dir, "{sample}.txt")
    log:
        os.path.join(logs_dir, "kofamscan_results_{sample}.log")
    shadow: "shallow"
    threads: 28
    conda:
        "kofamscan-env.yml"
    shell:
        """
        exec_annotation --ko-list {input.ko_list} --profile {input.profile_dir} --cpu {threads} -f mapper -o {output} {input.fasta}
        """
        #{input.executable} --ko-list {input.ko_list} --profile {input.profile_dir} --cpu {threads} -f detail -o {output} {input.fasta}
