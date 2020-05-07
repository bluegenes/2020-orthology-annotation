"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s eggnog.snakefile --use-conda #--cores 26 --cluster "sbatch -t 10:00:00 -N 1 -n 26 -p bmm --mem=60gb" --jobs 5 -k --rerun-incomplete 
"""

import os
import pandas as pd

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()


# need samples csv, with columns: "sample", "peptide_fasta" 
samplelist = config["samples_csv"]
samplesDF = pd.read_csv(samplelist, dtype=str, header=0, index_col = "sample")
SAMPLES=samplesDF.index.tolist()
out_dir = config.get("out_dir", "eggnog_out")
database_dir = config.get("db_dir", "eggnog_data")
db_version=config.get("eggnog_database_version", "5.0.0")

rule anNOGtate:
    input: 
         expand(os.path.join(out_dir, "{sample}.emapper.annotations"), sample=SAMPLES),
         expand(os.path.join(out_dir, "{sample}.emapper.seed_orthologs"), sample=SAMPLES)

rule get_eggnog_dbs:
    output:
        os.path.join(database_dir, "eggnog.db"),
        os.path.join(database_dir, "eggnog_proteins.dmnd"),
    conda: "eggnog.yml"
    params:
        data_dir = database_dir
    shell:
        """
        mkdir -p {params.data_dir}
        download_eggnog_data.py -y -f --data_dir {params.data_dir}
        """

## to optimize on cluster, run eggnog mapper steps separately
# https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2
# https://github.com/eggnogdb/eggnog-mapper/issues/80

# this part is cpu intensive
rule run_eggnog_mapper_dmnd:
    input:
        nog_db=os.path.join(database_dir, "eggnog.db"),
        nog_dmnd=os.path.join(database_dir, "eggnog_proteins.dmnd"),
        pep=lambda w: samplesDF.loc[w.sample, "peptide_fasta"]
    output:
        os.path.join(out_dir, "{sample}.emapper.seed_orthologs")
    resources:
      mem_mb=10000,
      runtime=6000,
    params:
        mode="diamond",
        data_dir=directory("eggnog_data"),
        out_dir=out_dir
    log: os.path.join(out_dir, "logs", "{sample}_emapper_dmnd.log")
    benchmark: os.path.join(out_dir, "logs", "{sample}_emapper_dmnd.benchmark")
    threads: 10
    shadow: "shallow" ## this means tmpdir will be on local scratch (using --shadow-prefix) --> faster? 
    conda: "eggnog.yml"
    shell:
        """
        emapper.py --no_annot --no_file_comments -i {input.pep} --output {params.out_dir}/{wildcards.sample} -m {params.mode} --cpu {threads} --data_dir {params.data_dir} > {log} 2>&1 
        """

# this part is I/O intensive --> USE LOCAL tmpdir! (use shadow, with --shadow-prefix /$SCRATCH/ntpierce)
rule run_eggnog_mapper_annotate:
    input:
        nog_db=os.path.join(database_dir, "eggnog.db"),
        seed_orthologs=os.path.join(out_dir, "{sample}.emapper.seed_orthologs")
    output:
        os.path.join(out_dir, "{sample}.emapper.annotations"),
    params:
        mode="diamond",
        data_dir=database_dir,
        out_dir=out_dir
    log: os.path.join(out_dir, "logs", "{sample}_emapper_annot.log")
    benchmark: os.path.join(out_dir, "logs", "{sample}_emapper_annot.benchmark")
    resources:
      mem_mb=12000,
      runtime=6000,
    threads: 1
    conda: "eggnog.yml"
    shadow: "shallow" ## snakemake will make a shallow copy of dir over in the shadow-prefix dir 
    shell:
        """
        emapper.py --annotate_hits_table {input.seed_orthologs} --no_file_comments -o {params.out_dir}/{wildcards.sample} --data_dir {params.data_dir} --cpu {threads} > {log} 2>&1
        """

## unused:

# alternatively, grab a gz database version
#def download_eggnog_db_gz:
#    input: lambda wildcards: HTTP.remote(f"http://eggnogdb.embl.de/download/emapperdb-{db_version}s/eggnog.db.gz", static=True, keep_local=True, allow_redirects=True)
#    output: os.path.join(database_dir, "eggnog.db.gz")
#    shell:
#        """
#        mv {input} {output}
#        """


## if running on single node, can run all at once
#rule run_eggnog_mapper:
#    input: 
#        nog_db="eggnog_data/eggnog.db",
#        nog_dmnd="eggnog_data/eggnog_proteins.dmnd",
#        pep=get_pep
#    output:
#        os.path.join(out_dir, "{sample}.emapper.annotations"),
#        os.path.join(out_dir, "{sample}.emapper.seed_orthologs")
#    params:
#        mode="diamond",
#        data_dir=directory("eggnog_data")
#    log: os.path.join(out_dir, "logs", "{sample}_emapper.log")
#    benchmark: os.path.join(out_dir, "logs", "{sample}_emapper.benchmark")
#    threads: 10
#    conda: "eggnog.yml"
#    shell:
#        """
#        emapper.py -i {input.pep} --output {wildcards.sample} -m {params.mode} --data_dir {params.data_dir} --cpu {threads} > {log} 2>&1
#        """
