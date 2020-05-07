# 2020-orthology-annotation
run existing orthology annotation programs

# eggnog-mapperv2

    # running on a compute node:
    snakemake -s eggnog.snakefile --configfile test_data/test_config.yml --cores 1 --cluster-config cluster-config.yml --use-conda -p -k --rerun-incomplete


    # running via job scheduler:
      # 1. Set up a snakemake default or cluster profile (see http://bluegenes.github.io/Using-Snakemake_Profiles/)
      # 2. Modify "shadow-prefix" and your the `cluster-config.yml` file to reflect your cluster. Info on shadow-prefix: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#shadow-rules
      # 3. Run: snakemake -s eggnog.snakefile --configfile test_data/test_config.yml --profile slurm --jobs 3 --cluster-config cluster-config.yml --shadow-prefix /scratch/ntpierce


# kofamscan 
    # running on a compute node:
    snakemake -s kofamscan.snakefile --configfile test_data/test_config.yml --cores 1 --use-conda -p -k --rerun-incomplete


    # running via job scheduler:
      # 1. Set up a snakemake default or cluster profile (see http://bluegenes.github.io/Using-Snakemake_Profiles/)
      # 2. Modify the `cluster-config.yml` file to reflect your cluster.
      # 3. Run: snakemake -s kofamscan.snakefile --configfile test_data/test_config.yml --profile slurm --jobs 3 --cluster-config cluster-config.yml

