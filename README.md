# de_novo-identification
Finding de novo repeats within a genome assembly and filtering out simple repeats and genes.

## Amarel
This pipeline is run on Amarel, the Rutgers University high performance computing server.
The following are the most recent specifications of the Amarel system that the pipeline was run on.
- 52 CPU-only nodes, each with 28 Xeon e5-2680v4 (Broadwell) cores + 128 GB RAM
- 20 CPU-only nodes, each with 28 Xeon e5-2680v4 (Broadwell) cores + 256 GB RAM
- 4 28-core e5-2680v4 nodes each with 2 x Nvidia Pascal P100 GPUs onboard
- 2 high-memory nodes, each with 56 e7-4830v4 (Broadwell) cores + 1.5 TB RAM
- 53 CPU-only nodes, each with 16 Intel Xeon e5-2670 (Sandy Bridge) cores + 128 GB RAM
- 5 CPU-only nodes, each with 20 Intel Xeon e5-2670 (Ivy Bridge) cores + 128 GB RAM
- 26 CPU-only nodes, each with 24 Intel Xeon e5-2670 (Haswell) cores + 128 GB RAM
- 4 CPU-only nodes, each with 16 Intel Xeon e5-2680 (Broadwell) cores + 128 GB RAM
- 3 12-core e5-2670 nodes with 8 Nvidia Tesla M2070 GPUs onboard
- 2 28-core e5-2680 nodes with 4 Quadro M6000 GPUs onboard
- 1 16-core e5-2670 node with 8 Xeon Phi 5110P accelerators onboard

For more information, please contact help@hpc.rutgers.edu or look at the Amarel user guide located at https://rutgers-oarc.github.io/amarel/.

You may also contact the Amarel director, Kristen Klepping at klepping@rutgers.edu for `sudo` permissions on your node.

## Run
The pipeline is run by using the following command.

`sbatch de_novo_repeat_shell.sh`

`de_novo_repeat_shell.sh` is a shell file that contains subcommands which runs the individual programs and tools.

## Absent files

Some files are absent from this repository due to the inability to store large files on the free version of Github.
Those files are the nucleotide sequence of the multiple fly species as well as their respective peptide sequences; all of which can be found within the FTP client of `flybase.org`.

## Program Version information
All the programs and pipelines used within this master are tht latest versions upto the last date of run. 

All of them have been kept up to date using Anaconda.

### Author Information
For more information and data, please contact me here, or at chinmay.rele@rutgers.edu.
# de_novo-identification
