"""Snakemake wrapper for Kallisto quant"""

__author__ = "Sebastian Kurscheid"
__copyright__ = "Copyright 2020, Sebastian Kurscheid"
__email__ = "sebastian.kurscheid@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Creating log
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Placeholder for optional parameters
extra = snakemake.params.get("extra", "")

# Allowing for multiple FASTQ files
n = len(snakemake.input.fastq)
k = 2
pe = 'R2' in snakemake.input.fastq
assert ( 
     (n % k == 0 and pe) or (pe == False and n >= 1)
 ), "input->sample n (paired-end) elements must be divisible by 2."

#TODO: change hardcoded PE suffix to configurable version
if n == 1:
    reads = "--single {}".format(*snakemake.input.fastq)
else:
    if 'R2' in snakemake.input.fastq:
        reads = "{}".format(*snakemake.input.fastq)
    else:
        reads = "--single " + ' '.join('%-2s' for i in snakemake.input.fastq)%tuple(snakemake.input.fastq)

shell(
    "(kallisto quant "  # Tool
    "{extra} "  # Optional parameters
    "--threads={snakemake.threads} "  # Number of threads
    "--index={snakemake.input.index} "  # Input file
    "--output-dir={snakemake.output} "  # Output directory
    "{reads} )"  # Input FASTQ files
    "{log}"  # Logging
)
