# PlantGenomicsGuide
Instructions for all major step in a plant genome sequencing project

---

This document provides exemplary step-by-step commands for ONT long-read sequencing analysis, including basecalling, genome assembly, gene annotation, and final data submission. This hands-on protocol covers the entire process, from public data acquisition to final output assembly, making it accessible even to beginners. Brief explanations accompany each step to ensure that readers cannot only follow the protocol but also understand the methodology and can adapt it to their own research. 

Before the actual tools are introduced (Chapters 0 to 7), a brief introduction to general Unix-based/-like system behaviour and installation of tools is provided (Chapters I and II).

## I. General Conventions

Before starting, please note the following basic conventions:

- All commands are intended to be run on a Unix-based/-like system via a command line interface (CLI).
- Commands are generally written on a single line; however, for readability, longer commands may be split across multiple lines using a backslash `\` at the end of each line. When entered into the CLI,these commands will execute as a single line.
- `#` denotes a non-executable comment.
- Example file paths are written as `/path/to`, which should be replaced with the full path to your specific file or directory.
- Commands may include output redirection, e.g., `>> /path/to/logfile.log 2>&1` at the end to redirect both standard output (stdout) and standard error (stderr) to a log file, or `2>/path/to.err.txt` to redirect only stderr to a specified file. Some programs, such as Dorado or samtools, output results to stdout, so it is necessary to redirect these outputs to a file. Using `>` will overwrite the file, while `>>` appends to it.
- The symbol `&` at the end of the command runs it in the background. Alternatively, you may use a terminal multiplexer like `tmux` (<https://github.com/tmux/tmux/wiki>), which may be more convenient for managing long-running processes.

Links to official documentation are provided for every tool, as they often include many more features and options that cannot be covered in this document. While the examples provided here are based on our experiences, we strongly recommend consulting the official documentation and built-in help to explore options best suited for your specific project.

Helpful options to access tool documentation include:

```bash
man TOOL                   # Opens the manual page, if available
TOOL -h                    # Prints a short help message
TOOL --help                # Prints a more extensive help message
cd TOOL && nano README.md  # Many tools provide a README file in the installation directory
```
It is advisable to install and use the latest stable version of each tool to ensure optimal performance and compatibility. 

---

## II. Recommended installation of tools:

There are several common ways to install bioinformatics tools:

- **Package Manager**: Tools can be installed via system package managers like apt (Debian/Ubuntu), though these versions may be outdated.
- **Precompiled binaries**: Download pre-built binaries from official websites or GitHub repositories.
- **Python environments**: Use `conda` to create isolated environments with specific versions of tools and their dependencies.
- **Docker containers**: Containers typically include all necessary dependencies and are highly reproducible, but they tend to consume more disk space.

The preferred installation method depends on tool availability and user preference. The following sections provide basic usage instructions for each method.

### Installation with apt

On Debian/Ubuntu systems, apt can be used to install some bioinformatics tools. Basic usage includes:

```bash
sudo apt update                         # Update package lists
sudo apt upgrade                        # Upgrade installed packages
sudo apt search samtools                # Search for available packages
sudo apt install samtools seqkit        # Install one or more tools
```

Depending on your system configuration, you may need superuser rights (using sudo).


### Manual Installation

If no package or conda environment is available, you can manually install a tool. Usually, tool developers provide clear installation instructions, typically involving downloading a binary release and extracting it:

Using Precompiled Bimaries:
```bash
tar -xvzf tool.tar.xz
cd tool/
```

If binaries are not working or the latest development version is needed, you can compile from source:

```bash
git clone https://github.com/author/tool.git
cd tool/
make # Follow instructions in the README or INSTALL.md 
```

Usually, the full path to the executable binary is required to run the tool. To make the tool accessible in your shell environment, you need to add it to your `PATH` variable. For the current bash session:

```bash
export PATH="/path/to/tool:$PATH"
```

To make this change permanent across all bash sessions, add the line to your shell startup script (e.g., `~/.bashrc`). To do this, you can use a text editor like `nano`:

```bash
nano ~/.bashrc 
# Add the following line at the end of the file
export PATH="/path/to/tool:$PATH"
# Save and exit, then reload
source ~/.bashrc
```


### Python environments

As many bioinformatic tools can be installed via conda, we use this to package manager to demonstrate virtualenv. Environments which can exist along and which can have different version of the same package. Some tools are dependent on other packages (dependency) while another tool might be dependent on the same package, but with another, putatively conflicting, version. The virtual environments can handle these differences. Additionally, you do not need to export the path variable to run a command installed via conda in a venv. Conda installation guide: <https://docs.conda.io/projects/conda/en/stable/user-guide/install/linux.html>

Common commands, listing existent environments, creating and entering environments, installation and updating:

```bash
conda env list                             # List existing environments
conda create -n YOUR_ENV_NAME              # Create a new environment
conda activate YOUR_ENV_NAME               # Activate the environment
conda deactivate                           # Deactivate the current environment
conda install TOOL                         # Install a tool in the active environment
conda update -n YOUR_ENV_NAME TOOL         # Update a specific tool
conda update -n YOUR_ENV_NAME --all        # Update all tools in the environment
conda update -n base conda                 # Update conda itself
```



### Docker containers

Official Docker installation guide: <https://docs.docker.com/engine/install/>

Docker can run tools in isolated containers. An example BUSCO command:

```bash
sudo docker run --rm \
 -v /path/to:/path/to \
 ezlabgva/busco:v6.0.0_cv1 \
 busco --plot /path/to/busco_runs
```

In  this example

- `-v /path/to:/path/to` mounts your local directory inside the container, making files accessible.

- `ezlabgva/busco:v6.0.0_cv1` specifies the Docker image and version.

- Everything after the image name (`busco --plot ...`) works similarly to a local installation.

Most software projects provide ready-to-use Docker images. Using Docker can simplify installation but requires more disk space.

---

## 0. Downloading read data with prefetch and fasterq-dump

SRA-Toolkit Wiki: <https://github.com/ncbi/sra-tools/wiki>

Example datasets can be downloaded and used to practice the commands described in this manual. In this example, we use *Victoria cruziana* sequencing data, which is publicly available on NCBI’s Sequence Read Archive (SRA). The recommended tools to download the data are `prefetch` and `fasterq-dump` from the SRA Toolkit.

- `prefetch` downloads data in compressed `.sra` format.
- `fasterq-dump` converts `.sra` files to standard FASTQ format.

First, use `prefetch` to download the `.sra` files. You can specify one or multiple SRA run IDs (e.g., all IDs from BioProject PRJEB63973). Common options include:

- `-p` to show progress,

- `--max-size` to set the maximum size for downloads,

- `--output-directory` to define where downloaded files are stored.

```bash
prefetch -p \
 --max-size 250g \
 --output-directory /path/to/raw_reads/prefetch \
 ERR13955440 ERR13955441 ERR13955442 [...]
```

Next, extract FASTQ files from the downloaded .sra files using fasterq-dump.

- `e` specifies the number of parallel threads (adjust based on your available CPU cores),

- `--outdir` specifies where to write the FASTQ files.

```bash
fasterq-dump -e 30 /path/to/raw_reads/prefetch/*/*.sra --outdir /path/to/raw_reads
```

---

## 1. Basecalling

Dorado documentation: <https://github.com/nanoporetech/dorado>

Basecalling refers to converting raw signal data from ONT sequencers into nucleotide sequences. We use Oxford Nanopore’s Dorado, which provides GPU-accelerated basecalling with options for detecting modified bases.

Before basecalling, download the appropriate model for your sequencing chemistry using `dorado download`. You can list available models and options with `dorado basecaller --help`. In this example, we use the `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` model, which includes modified base detection for 5-methylcytosine (5mCG) and 5-hydroxymethylcytosine (5hmCG). The basecalling input is the folder containing your raw `.pod5` files, and the output is a BAM file. Logs are redirected to an error log file.

```bash
dorado basecaller \
 dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
 /path/to/Victoria_cruziana_run01/ \
 --modified-bases 5mCG_5hmCG \
 > /path/to/Victoria_cruziana_run01.mod.bam \
 2> /path/to/YourArchiveForBasecalling.dorado.mod.err.txt &
```

To convert the resulting BAM file into a FASTQ file (commonly used for downstream analyses), use samtools fastq.

```bash
samtools fastq -T "*" \
 /path/to/Victoria_cruziana_run01.mod.bam\
 > /path/to/Victoria_cruziana_run01.mod.fastq
```
- T "*" ensures all auxiliary tags in the BAM file are handled properly.

It is advisable to compress FASTQ files to save disk space. `gzip` compression is widely accepted by downstream bioinformatics tools. 

```bash
gzip /path/to/Victoria_cruziana_run01.mod.fastq
```
You can obtain basic quality statistics, including the read N50, using tools like seqkit (<https://bioinf.shenwei.me/seqkit/>) or a Python script (`FASTQ_stats3.py`; <https://github.com/bpucker/GenomeAssembly/blob/main/FASTQ_stats3.py>)

```bash
seqkit stats -N 50 /path/to/Victoria_cruziana_run01.mod.fastq.gz

python3 /path/to/FASTQ_stats3.py \
--in /path/to/Victoria_cruziana_run01.mod.fastq.gz
```

### 1.1 Raw read correction with HERRO

HERRO repository: <https://github.com/lbcb-sci/herro>

Raw long reads from ONT sequencing often include some errors. Correction improves accuracy before assembly. HERRO is integrated into Dorado and applies a two-step correction:

- Step 1 (CPU intensive) computes overlaps between reads,

- Step 2 (GPU intensive) corrects reads based on these overlaps.

While both steps can be conducted in one command, they can also be run separately to be able to utilize different hardware options to accelerate the running time of the tool due to proper hardware. If you have multiple runs for the same sample (e.g., multiple flowcells), you should merge them into a single file before correction. The merged file can remain compressed.

```bash
cat\
 /path/to/Victoria_cruziana_run01.mod.fastq.gz \
 /path/to/Victoria_cruziana_run02.mod.fastq.gz \
 /path/to/Victoria_cruziana_run03.mod.fastq.gz \
 > /path/to/Victoria_cruziana.fastq.gz &  
```

The first step calculates overlaps between reads and outputs a `.paf` file (specified with `--to-paf`) describing the overlaps. The output, by default, is in stdout and needs to be redirected into a file (using `>`). Use `--threads` to optimize performance based on CPU availability. Errors and logs are saved in *.doc.txt.

```bash
/path/to/dorado correct \
 /path/to/Victoria_cruziana.fastq.gz \
 --to-paf --threads 10 \
 > /path/to/Victoria_cruziana.overlaps.paf \
 2> /path/to/Victoria_cruziana.overlaps.doc.txt &
```

For the second step (GPU intensive), only uncompressed input FASTQ files are supported. To uncompress them simply use `-gunzip /path/to/Victoria_cruziana.fastq.gz`. Add `-k` if you want to keep the original compressed file. Then, run Dorado with HERRO correction. `--form-paf` specifies the overlap file, and corrected reads are output in FASTA format.

```bash
/path/to/dorado correct /path/to/Victoria_cruziana.fastq \
 --from-paf /path/to/Victoria_cruziana.overlaps.paf \
 > /path/to/Victoria_cruziana.corrected_reads.fasta \
 2> /path/to/Victoria_cruziana.corrected_reads.errors.txt &
```
---

## 2. Assembly

### 2.1 Shasta

Official Shasta documentation: <https://paoloshasta.github.io/shasta/>

Shasta is a fast and efficient assembler for ONT long reads. It uses advanced algorithms to assemble genomes quickly while requiring comparatively modest computational resources. Depending on the type of flow cell used (e.g., R9.4.1, R10.4.1) and whether the reads have undergone prior correction, different pre-configured settings (configuration files) should be used to optimize performance. Key parameters are `--input` specifying the input reads (in FASTA or FASTQ format), `--config` for suitable configuration file, `--assemblyDirectory` defines the output folder where the assembly results will be written, and `--threads` sets the number of CPU threads to use.

```bash
/path/to/binary/shasta-Linux-0.14.0 \
 --threads 10 \
 --input /path/to/Victoria_cruziana.corrected_reads.fasta \
 --config /path/to/Nanopore-r10.4.1_e8.2-400bps_sup-Herro-Jan2025_red_bp.conf \
 --assemblyDirectory /path/to/shasta/Vcruz_01
```

It is recommended to add optional arguments `--memoryBacking 2M` and `--memoryMode filesystem` to enable Shasta to back memory allocations with 2 MB huge pages and use the filesystem for temporary storage. This can speed up performance and reduce disk I/O, but may require root or superuser permissions depending on your system setup.

### 2.2 NextDenovo

Official NextDenovo documentation: <https://nextdenovo.readthedocs.io/en/latest/>

NextDenovo is a long-read assembler suitable for large genomes. It performs both read correction and assembly. NextDenovo requires a configuration file that specifies runtime parameters, making it easy to reproduce results or adjust settings for different datasets.

To run NextDenovo, you simply provide the configuration file:

```bash
/path/to/binary/nextDenovo /path/to/nextdenovo/Vcruz_01.cfg
```

The config file can be copied from the documentation. Key parameters in the config file include:
- `input_fofn`: a file listing the paths to input read file
- `input_type`: `raw` for uncorrected read or `corrected` for pre-corrected reads
- `read_type`: specifies the sequencing technology (ONT, PacBio CLR, or PacBio HiFi)
- `genome_size`: estimated genome size, which guides the assembly process, it is advised to adjust parameters regarding computation resources (`-t` and `-p` options).
- `parallel_jobs`: controls how many parallel processes NextDenovo will run
- `workdir`; directory where intermediate files and the final assembly is stored

Example configuration for *Victoria Cruziana*

```config
[General]
job_type = local
job_prefix = nextDenovo
task = all # "assemble" can be used if read correction was already performed
rewrite = yes
deltmp = yes
parallel_jobs = 6
input_type = raw # alternatively "corrected"
read_type = ont # clr, ont, hifi
input_fofn = /path/to/input.fofn # the .fofn is a file containing the input read file paths
workdir = /path/to/nextdenovo/Vcruz_01 # output directory

[correct_option]
read_cutoff = 1k
genome_size = 4g # estimated genome size
sort_options = -m 20g -t 15
minimap2_options_raw = -t 8
pa_correction = 3
correction_options = -p 15

[assemble_option]
minimap2_options_cns = -t 4
nextgraph_options = -a 1
```

Example `.fofn` file listing input FASTQ files:

```bash
/path/to/Victoria_cruziana.run1.fastq
/path/to/Victoria_cruziana.run2.fastq
/path/to/Victoria_cruziana.run3.fastq
```

### 2.3 Verkko

Official Verkko documentation: <https://github.com/marbl/verkko>

Verkko is a versatile assembler that combines multiple data types, including HERRO corrected ONT reads or PacBio HiFi reads via `--hifi`, ultra-long ONT reads via `--nano`, and Pore-C data via `--porec`.

A minimal Verkko run includes `-d` to set the output directory, `--hifi` to specify HERRO corrected ONT reads (or PacBio HiFi reads), and optional `--nano` or `--porec` for ultra-long ONT (Dorado basecalled reads) or Pore-C reads, respectively. If `--porec` is given, the correct telomere motif should be specified with `--telomere-motif`. Nowak *et al*. (2025, <https://doi.org/10.1101/2024.06.15.599162>) reported that Pore-C scaffolding produced better results when conducted with alternative tools like CPhasing. Additionally,`--local-memory` and `--local-cpus` to specify the upper limit of memory in GB (default 64) and the number of CPUs (standard all), respectively can be specified.

```bash
verkko --par-run 27 464 48 -d /path/to/verkko/Vcruz_01 --local-memory 464 --local-cpus 27\
 --hifi /path/to/Victoria_cruziana.corrected_reads.fasta
```

The `--par-run` option allows you to modify CPU, memory (GB), and runtime (hours) limits for specific pipeline stages. For example, `--par-run 27 464 48` allocates 27 CPUs, 464 GB RAM, and 48 hours for one of the internal pipeline steps. This is especially helpful if you encounter bottlenecks in particular stages of the assembly.

For more fine-grained control, consult the full list of Verkko parameters with `verkko --help`.

### 2.4 Hifiasm

Official Hifiasm documentation: <https://github.com/chhylp123/hifiasm>

Hifiasm is a modern, ultra-fast assembler primarily designed for PacBio HiFi and ONT R10 reads. For ONT reads, it is essential to use the `--ont` option. Hifiasm only accepts reads in FASTQ format. 
Key paramenters are `-t` for defining the number of CPU threads to use, `--ont` to specify ONT reads, and `-o` to specify the output file prefix.

```bash
hifiasm -t 64 --ont -o ONT.asm /path/to/Victoria_cruziana.corrected_reads.fastq.gz
```

Hifiasm automatically outputs primary and alternative assemblies as GFA files, which can be converted to FASTA using tools such as `awk` or `gfatools`.

### 2.5 CPhasing (Scaffolding with Pore-C data)

Official CPhasing documentation: <https://wangyibin.github.io/CPhasing/latest/>

CPhasing is a tool for chromosome-scale scaffolding of genome assemblies using Pore-C data. It enhances contiguity by leveraging long-range contact information. Key parameters include `-f` for draft assembly file (FASTA), `-pcd` for Pore-C reads which can be gzipped, `-t` for number of cores, `-p` for restriction enzyme site used in Pore-C library preparation, and `-n` for chromosome count as `haploid_number:ploidy` (e.g., `12:1` for 12 chromosomes, haploid).`-hcr` is an optional argument that specifies that only high confidence regions should be retained for scaffolding. 

```bash
cphasing pipeline -f /path/to/Victoria_cruziana.genome.fa\
 -pcd /path/to/Victoria_cruziana.porec.fastq.gz -t 24 -n 12:1 -hcr -p CATG
```
The output will be written to a new directory created in the current working folder. This directory contains both the improved assembly and diagnostic reports on scaffolding quality.


---


## 3. Assembly Evaluation

After running multiple assemblers, it is essential to evaluate their performance and select the most accurate contig-level assembly for downstream scaffolding. The final assembly should be comprehensively assessed to ensure completeness, accuracy, and overall quality.

### 3.1 Basic Assembly Statistics

To calculate basic statistics of the assembly such as N50, total length, number of contigs, GC content and other, the lightweight python script contig_stats.py (<https://github.com/bpucker/script_collection/blob/master/contig_stats.py>) can be used. It can optionally filter out short contigs using --min_contig_len.

```bash
python3 contig_stats.py \
 --input /path/to/Victoria_cruziana.genome.fa \
 --min_contig_len 10000 #in basepairs \
 --out /path/to/output_directory/
```

### 3.2 BUSCO (Benchmarking Universal Single-Copy Orthologs)

Documentation: <https://busco.ezlab.org/busco_userguide.html>

BUSCO evaluates genome completeness by identifying near-universal single-copy orthologs from curated datasets (OrthoDB). The first step is to select appropriate lineage dataset. To list available dayasets:

```bash
busco --list-datasets
```
The output will list current datasets which can be automatically downloaded by BUSCO. BUSCO authors suggests to use the most specific datastet to get highest-resolution analysis. For example, the full lineage of *Victoria cruziana* is (`cellular organisms/Eukaryota/Viridiplantae/Streptophyta/Streptophytina/`  
`Embryophyta/Tracheophyta/Euphyllophyta/Spermatophyta/Magnoliopsida/`  
`Nymphaeales/Nymphaeaceae/Victoria`). Based on currently available BUSCO_odb12 datasets, `eukaryota_odb12`, `viridiplantae_odb12`, and `embryophyta_odb12` datasets can be used. Thus, `embryophyta_odb12` is the current best-fitting dataset for *V. cruziana*.

```bash
- archaea_odb12 [195]
    - euryarchaeota_odb12 [282]
        - methanomicrobia_odb12 [589]
            - methanosarcinaceae_odb12 [971]
            - methanosarcina_odb12 [1620]
...
    - viridiplantae_odb12 [822]
        - embryophyta_odb12 [2026]
            - eudicotyledons_odb12 [2805]
                - brassicales_odb12 [4311]

...
```

For the actual BUSCO run, the lineage dataset can be specified with the `-l` or `--lineage_dataset` option. `--mode` is mandatory and can be `genome`, `proteins` or `transcriptome` in the case of a genome fasta, peptide fasta, or CDS FASTA, respectively. `-i` is input FASTA file and recommended, but not mandatory, is the specification of an output folder and an output folder name. We recommend the output folder to be the same per project or machine running BUSCO, with only the output folder name being changed per run. Thus, all run-specific output folders are collected in a single BUSCO output folder (The given examples are using specific run IDs). To speed up the run, more resources (number of CPUs) can be specified with `--cpu` (default `1`).

```bash
# BUSCO run for an assembly:
busco -i /path/to/Victoria_cruziana.genome.fa -m genome\
 --cpu 27 -l embryophyta_odb12 --out_path /vol/data/sam/BUSCO_runs/ -o SNM_BIR_0079

# BUSCO run for proteins derived from a structural annotation:
busco -i /path/to/Victoria_cruziana.pep.fa -m proteins --cpu 27\
 -l embryophyta_odb12 --out_path /vol/data/sam/BUSCO_runs/ -o SNM_BIR_0134
```
BUSCO reports “Complete”, “Fragmented”, and “Missing” BUSCOs to assess assembly quality.


### 3.3 LTR Assembly Index (LAI)

Documentation: <https://github.com/oushujun/LTR_retriever>
LAI measures genome assembly continuity based on intact Long Terminal Repeat Retrotransposons (LTR-RTs). This involves the following steps:


**Step 1:** Index the genome using GenomeTools (gt)

```bash
/path/to/gt suffixerator \
 -db /path/to/Victoria_cruziana.genome.fa \
 -indexname Victoria_cruziana.genome.fa \
 -tis -suf -lcp -des -ssp -sds -dna
```

**Step 2:** Identify candidate LTR elements using similarity searches

```bash
/path/to/gt ltrharvest \
 -index Victoria_cruziana.genome.fa \
 -minlenltr 100 -maxlenltr 7000 \
 -mintsd 4 -maxtsd 6 \
 -motif TGCA -motifmis 1 \
 -similar 85 -vic 10 -seed 20 \
 -seqids yes > Victoria_cruziana.genome.fa.harvest.scn
```

**Step 3:** Detects additional LTR candidates using structural features

```bash
/path/to/LTR_FINDER_parallel \
 -seq Victoria_cruziana.genome.fa \
 -threads 10 -harvest_out -size 1000000 -time 300
```

**Step 4:** Combine Results

```bash
cat Victoria_cruziana.genome.fa.harvest.scn Victoria_cruziana.genome.fa.finder.combine.scn \
 > Victoria_cruziana.genome.fa.rawLTR.scn
```

**Step 5:** LTR_retriver calculates LAI and filters high confidence LTR-RTs.

```bash
/path/to/LTR_retriever -genome Victoria_cruziana.genome.fa\
 -inharvest Victoria_cruziana.genome.fa.rawLTR.scn -threads 10 
```

The final LAI score is reported in the .out.LAI file. LAI scores between 0 to 10 indicate 'draft genome', between 10-20 'Reference genome' and more than 20 indicates 'gold standard'.

### 3.4 Mapping-based Genome Size Estimation (MGSE)

Documentation: <https://github.com/bpucker/MGSE>

MGSE estimates genome size based on k-mer coverage of reference gene sets (e.g., BUSCO genes). MGSE requires the assembly FASTA file and the corrected reads in FASTQ format. Alternatively, if you have a BAM file with the reads mapped to the assembly, that can be provided instead of both assembly and reads file. This BAM file can be generated using MINIMAP2. The reference genes can be provided as gff3 file. MGSE authors state that the BUSCO genes appears to be the best choice for reference genes. For this the `full_table.tsv` is used which can found in the BUSCO output directory.

```bash
 python3 MGSE3.py --fasta /path/to/Victoria_cruziana.genome.fa \
 --fastq /path/to/Victoria_cruziana.corrected_reads.fasta \ #here, a directory with multiple
 #fastq files can also be provided
 #alternatively, BAM file with --bam can be provided
 --out /path/to/MGSE_output/
 --ref /path/to/BUSCO_genome/run_embryophyta_odb12/full_table.tsv \
 --threads 10
```

### 3.5 Merqury

Documentation: <https://github.com/marbl/merqury>

Merqury performs a reference-free quality assessment using k-mer-based approaches to calculate consensus accuracy (QV score) and completeness. It is important to select the right k-mer size for assessing the quality. If unsure about the right k size, run `best_k.sh` script provided with MERQURY installation. It just requires the `genome_size` in bases (e.g., XXX for _Victoria_cruziana_). Next, Build k-mer database with `meryl` followed by running Merqury. 

```bash
# Calculate the right k size
sh /path/to/MERQURY/best_k.sh <genome_size>

# Count k-mers from raw reads using meryl
 meryl k=21 count /path/to/Victoria_cruziana_corrected_reads.fastq \
 output_Victoria_cruziana.meryl >>/path/to/meryl.log 2>&1 &
 
# Run Merqury to assess assembly quality
 merqury.sh output_Victoria_cruziana.meryl\
 /path/to/Victoria_cruziana.genome.fasta >> /path/to/merqury.log 2>&1 &
```

Higher QV indicates fewer errors (QV40=1 error per 10,000 bases). 

### 3.6 QUAST

Documentation: <https://quast.sourceforge.net/docs/manual.html#sec2>

QUAST provides a comprehensive summary of assembly quality, reporting metrics like N50, number of contigs, GC content, and misassemblies (if a reference genome is provided).

```bash
 /path/to/QUAST/quast.py -o /path/to/quast_results/\
 --eukaryote #default: prokaryote \
 --large #given if genome is larger than 100 Mbp \
 --circos # only give if you want to switch on the circos plotting
 --threads 10 \
 --labels "Victoria cruziana" #Assembly name to be used in reports. Quotes should be used if
 #label name include spaces.
```

If a reference genome is available, it can be added with `--reference` to detect misassemblies.

---


## 4. Structural Annotation


## 4.1 Protein-coding genes


### 4.1.1 RNA-Seq mapping with HISAT2

Documentation: <https://daehwankimlab.github.io/hisat2/>

For gene prediction pipelines like GeMoMa or BRAKER3, mapping RNA-Seq data provides essential "hints" about gene structures, such as exon-intron boundaries. HISAT2 is a splice-aware aligner that efficiently maps RNA-Seq reads to a genome assembly. Other tools like STAR (<https://github.com/alexdobin/STAR>) can be used as well. 

Before mapping, HISAT2 requires creating an index from the genome assembly FASTA file. Input genome assembly in FASTA format and output prefix (VC for Victoria_cruziana) for the generated index files are required. 

```bash
hisat2-build -p 10 \ # -p specifies the number of cores
/path/to/Victoria_cruziana.genome.fa VC
```
Indexing creates several auxiliary files (`.ht2` files) that are required for mapping.

After indexing, RNA-Seq reads (paired-end or single-end) are aligned to the genome. The input reads can be a comma separated list of file paths and may be gzip compressed. For paired end reads the arguments `-1` and `2` are used and for single-end reads `-U` is used. The output can be piped (`|`) to samtools to immediately sort the resulting mapping by genomic coordinates and save it in indexed BAM format for downstream usage. The output BAM file will be used as RNA-Seq 

```bash
hisat2 -p 10 \ 
 -x /path/to/Victoria_cruziana.genome.fa.index \
 -1 /path/to/RNA-Seq_reads_001_1.fq.gz,/path/to/RNA-Seq_reads_002_1.fq.gz \
 -2 /path/to/RNA-Seq_reads_001_2.fq.gz,/path/to/RNA-Seq_reads_002_2.fq.gz \
 | samtools sort -@ 10 \ # number of cores specified with @
 -o /path/to/Victoria_cruziana.genome_RNA-Seq_mapping.bam \
 && samtools index /path/to/Victoria_cruziana.genome_RNA-Seq_mapping.bam
```

### 4.1.2 BRAKER3

Documentation: <https://github.com/Gaius-Augustus/BRAKER>

BRAKER3 is a fully automated pipeline for annotating protein-coding genes using RNA-Seq data and/or protein sequences as hints. It combines multiple gene prediction tools like GeneMark and AUGUSTUS.

#### 4.1.2.1 Obtaining protein hints

Protein hints guide BRAKER3 to make more accurate predictions. The source for these protein hints can be the precompiled protein hints from OrthoDB 11 (<https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/> ;`Viridiplantae`). Additional peptide sequences from UniProt database can be included. Using advanced UniProt search, a broader taxon can be given via the 'Taxonomy' search field. Example, in the case of *Victoria cruziana*, its order Nymphaeales \[261007\](which translates to the search key "(taxonomy_id:261007)"). Multiple FASTA files as a comma-separated list for the `--prot_seq` parameter can be given.

#### 4.1.2.2 BRAKER3 command

Input arguments are `--genome` for genome assembly in FASTA format, `--threads` to define number of CPU threads to use, `--workingdir` to define the output directory, `--gff3` output annotations in gff3 format, `--bam` defines path to RNA-seq BAM alignment file, `--prot_seq` protein sequences for hints, and `--busco_lineage` defines the BUSCO lineage dataset for training AUGUSTUS. 

```bash
 /path/to/braker.pl --genome=/path/to/Victoria_cruziana.genome.fa --threads=10 \
 --workingdir=/path/to/braker3/Vcruz_annotation01/ --gff3 \
 --bam=/path/to/Victoria_cruziana.genome_RNA-Seq_mapping.bam \
 --prot_seq=/path/to/Viridiplantae.fa,/path/to/Uniprot_hints.fasta \
 --busco_lineage=embryophyta_odb12
```

### 4.1.3 GeMoMa

Documentation: <https://www.jstacs.de/index.php/GeMoMa> and <https://www.jstacs.de/index.php/GeMoMa-Docs>

GeMoMa uses homology-based gene prediction by aligning genes from closely related reference species and optionally includes RNA-Seq evidence.
The first step is to obtain reference genomes and their structural annotation in GFF/GFF3 file format from closely-related species. For our example *Victoria cruziana*, closely related species, *Nymphaea colorata* (GenBank ID: GCA_008831285.2) and *Nymphaea thermarum* (GenBank ID: GCA_011799765.1) are used. 

In the below example run, each reference species is given with 4 parameters. `s=own` specifies custom reference species, `i=` is an optional identifier tag for reference, the genome assembly and GFF3 file are provided with `g=` and `a=`, respectively. `t=` is used to specify the path of the target genome sequence, `outdir=` to give the output directory path, and `threads=` for specifying the number of threads to use. `r=MAPPED` and `ERE.m=` specifies that the RNA-Seq hints are mapped and the path to the mapping file. `ERE.m=` can be speciefied multiple times. Both of these arguments can be skipped if RNA-Seq data is not available. The default output is GFF file only. `pc=true` outputs FASTA coding sequences, `p=true` output FASTA protein sequences  `pgr=true` specifies that FASTA genomic regions should be returned, and `o=true` specifies that the predictions per individual reference should also be returned. This is important for further filtering steps. `Extractor.r` specifics that GeMoMa will try to repair transcript annotations that cannot be parsed . `GAF.f` specifies a filter string that filters the output. In the below example, only predictions with proper start (start==`M`) and codons (stop==`*`) are kept. Additionally, if a prediction has a certain score to length of amino acid ratio, it will be kept (score/aa>='0.75'). `AnnotationFinalizer.r=NO` option is skips renaming.

```bash
java -jar /path/to/GeMoMa-1.9.jar CLI GeMoMaPipeline \ 
 t=/path/to/Victoria_cruziana.genome.fa \
 s=own i=NyCol a=/path/to/ref_species/GCF_008831285.2/genomic.gff \
 g=/path/to/ref_species/GCF_008831285.2/GCF_008831285.2_ASM883128v2_genomic.fna \
 s=own i=NyTher a=/path/to/ref_species/GCA_011799765.1/genomic.gff \
 g=/path/to/ref_species/GCA_011799765.1/GCA_011799765.1_ASM1179976v1_genomic.fna r=MAPPED \
 ERE.m=/path/to/Victoria_cruziana.genome_RNA-Seq_mapping.bam \ 
 outdir=/path/to/gemoma/Vcruz_annotation02/ \
 pc=true pgr=true p=true o=true Extractor.r=true \
 GAF.f="start=='M' and stop=='*' and (score/aa>='0.75')" \
 AnnotationFinalizer.r=NO threads=27 >> GeMoMa.log 2>&1 \
```
The resulting output directory will have multiple files. The most important are individual predictions for each reference as a separate prediction GFF file and combined prediction ,using the filter specified, from all references. In case you do not want to merge external annotations and GeMoMa annotation is already good enough, then you can further filter it (if number of genes are too high) and proceed with the finalizing steps described in section 4.1.5


### 4.1.4 Funannotate

Documentation: <https://funannotate.readthedocs.io/en/latest/index.html>

The Funannotate pipeline can perform gene prediction, functional annotation, and comparison. Here, we will only focus on gene prediction. There are multiple steps to generate a Funannotate prediction.

1. First, optional step is to clean your assembly. If you have a haploid assembly, `funannotate clean` can remove some repetitive small, low-quality contigs:

```bash
funannotate clean -i /path/to/Victoria_genome.fa \
 -o /path/to/Victoria_cruziana_cleaned.genome.fa \
 -m 500 #min. length of contig to keep, by default 500
```
2. Second mandatory step is to softmask the repetitive elements in the assembly using `tantan`

```bash
funannotate mask -i /path/to/Victoria_cruziana_cleaned.genome.fa \
 -o /path/to/Victoria_cruziana_cleaned_softmasked.genome.fa \
 --cpus 10 # No. of cpus to use
```
3. Now, `funannotate train` uses RNA-Seq data from the species you are annotating and generates genome-guided `Trinity` assembly followed by `PASA` assembly. If RNA-Seq data is not available, this step does _de novo_ `Augustus` training. It produces the input data for funannotate predict, i.e. coord-sorted BAM alignments, trinity transcripts, and high quality PASA gff3 annotation. Single-end fastq files could also be provided with `-s` argument.

```bash
funannotate train -i /path/to/Victoria_cruziana_cleaned_softmasked.genome.fa \
 -o /path/to/funannotate_output/ #Output folder \
 -l /path/to/RNA-Seq_reads_001_1.fq.gz /path/to/RNA-Seq_reads_002_1.fq.gz \
 -r /path/to/RNA-Seq_reads_001_2.fq.gz /path/to/RNA-Seq_reads_002_2.fq.gz \
 --cpus 10 \
 --species "Victoria cruziana"
 --max_intronlen 1000000
```
4. The next step is to use `funannotate predict` to predict gene models. It automatically uses the results from the `train` step and incorporates them into the prediction, given the same output directory as `train` step. Additional, protein evidence from closely related species can be used with `--protein_evidence` argument.

```bash
funannotate predict -i /path/to/Victoria_cruziana_cleaned_softmasked.genome.fa \
 -o /path/to/funannotate_output/ #Output folder \
 -s "Victoria cruziana" --busco_db embryophyta \
 --organism other \
 --protein_evidence /path/to/ref_species/GCA_008831285.2/peptide.fasta /path/to/ref_species/GCA_011799765.1/peptide.fasta \
 --cpus 10 \
 --max_intronlen 1000000 \
```
5. Finally, `funannotate update` command can be used to add UTRs and refine gene models RNA-Seq data. 

```bash
funannotate update -i /path/to/funannotate_predict/ \
 --species "Victoria_cruziana" \
 -o path/to/funannotate_output/ \
 --cpus 10 \
```

### 4.1.5 Combining and Finalizing Annotations using GeMoMa

It is important to evaluate the quality of annotations from different sources which can vary in quality and completeness. To evaluate the completeness, BUSCO score using same lineage dataset and in protein mode should be calculated. The goal is to get the closest BUSCO score to the genome (same or higher). If any annotation recovers the genome BUSCO, it might be sufficient and do not require merging. However, if that is not the case, evaluating, filtering, and combining these annotations would be beneficial to create a final high-confidence gene set. Below is a step-by-step guide to combine annotations from different sources:

**Step 1:** Prepare RNA-Seq evidence with GeMoMa's ERE

Although all three gene annotation tools described above can incorporate `GFF3` from external sources, GeMoMa has been observed to give the best results. To compare annotations from different sources, they should be compared on the same scale. GeMoMa's gene models have some attributes that define the quality of the gene models. To add these attributes to other annotations, we need to extract RNA-Seq coverage and intron information from the RNA-Seq alignments. GeMoMa's `ERE` module is used for this where `m=` dspecifies the input RNA-Seq alignment in BAM format sorted by coordinates and `outdir=` specifies the directory where ERE will store its output.

```bash
java -jar /path/to/GeMoMa-1.9.jar CLI ERE \ 
 m=/path/to/Victoria_cruziana.genome_RNA-Seq_mapping.bam \
 outdir=/path/to/gemoma/Vcruz_RNA_evidence
```
This command creates `coverage.bedgraph` (read coverage across the genome) and `introns.gff` (detected intron positions). These files are used to add quality attributes to gene models based on expression data.

**Step 2:** Add attributes to external annotations with AnnotationEvidence

Now, gene annotations (e.g., from BRAKER3 or Funannotate) are enriched with RNA-seq based metrics. In the below command, `a=` is the external annotation GFF3 file, `g=` is the genome assembly FASTA file, `c=` defines the coverage file (BAM format and can be `UNSTRANDED`, `STRANDED`, and `NO`). In `UNSTRANDED`, only one `.bedgraph` file is needed, while for `STRANDED`, two `.bedgraph` file are required under `coverage_forward` and `coverage_reverse` attribute. `i=` points to intron evidence, and `outdir=` specifies the output directory.
Example with BRAKER3 output:
```bash
java -jar /path/to/GeMoMa-1.9.jar CLI AnnotationEvidence \ 
 a=/path/to/braker3/Vcruz_annotation01/braker.gff3 \
 g=/path/to/Victoria_cruziana.genome.fa \
 c=UNSTRANDED \
 coverage_unstranded=/path/to/gemoma/Vcruz_RNA_evidence/coverage.bedgraph \
 i=/path/to/gemoma/Vcruz_RNA_evidence/introns.gff \
 outdir=/path/to/gemoma/Vcruz_braker_anno_with_evidence/
```

This command will produce among other files, the `annotation_with_attributes.gff` file, with the added attributes. Each gene has useful attributes, such as average RNA-Seq coverage per gene (`avgCov`), fraction of transcript's introns supported by RNA-Seq (`tie`), amino acid sequence length (`aa`) among multiple others. Note that `score` is confidence score for GeMoMa predictions and therefore, external annotations do not have this. These attributes allows objective filtering in the next step based on biological evidence.

**Step 3**: Filter and merge annotations using GeMoMa GAF

Gene annotations often include low-confidence or biologically implausible predictions (e.g., unsupported isoforms, incomplete genes). GeMoMa’s GAF (GeMoMa Annotation Filter) can filter low-quality genes based on RNA-Seq support and other metrics, combine multiple annotations (e.g., GeMoMa + BRAKER3 + Funannotate), and remove redundancy and select the best gene models.
In the GeMoMa Annotation Filter (`GAF`) command. `g=` specifies input GFF3 files (can be used multiple times for multiple annotations), `f=` sets the filtering criteria, for example, `isNaN(score)` meaning score is not available, is necessary to include external annotations since they do not have the score attribute. Multiple filtering criteria can be combined with the logical use of `and`, `or` and brackets `()`. `atf=` filters alternative isoforms (optional but recommended). `outdir=` sets the output directory.
 
```bash
java -jar /path/to/GeMoMa-1.9.jar CLI GAF\ 
 g=/path/to/gemoma/Vcruz_braker_anno_with_evidence/filtered_predictions.gff\
 g=/path/to/gemoma/Vcruz_annotation02/unfiltered_predictions_from_species_0.gff\
 g=/path/to/gemoma/Vcruz_annotation02/unfiltered_predictions_from_species_1.gff\
 f="start=='M' and stop=='*' and aa>=15 and avgCov>0 and (isNaN(score) or score/aa>='3.25')"\
 atf="tie==1 or sumWeight>1" outdir=/path/to/gemoma/Vcruz_filter_01
```

In the above command, `start=='M' and stop=='*'` only keeps genes with proper start and stop codons, `aa>=15` removes predictions encoding very short peptides (likely spurious), `avgCov>0` keeps only genes with any RNA-Seq support, `is(NaN(score) or score/aa>='3.25')` keeps genes with good normalized scores for GeMoMa genes; for external tools (which don't have a GeMoMa score), the filter allows inclusion by skipping score check. `atf="tie=1 or sumWeight>6"` filters alternate transcripts to retain those with complete intron support (`tie==1`) or strong overall support (`sumWeight>1`). 
These filters can (and should) be customized based on particular project. Always check BUSCO scores (protein mode) after filtering steps to validate completeness and investigate annotation statistics (total number of genes and transcripts). It is recommended to iteratively test stricter or more relaxed filter criteria for optimal balance between gene count and completeness.

**Step 4:** Rename genes using AnnotationFinalizer

Once an annotation is merged and sufficiently filtered (reasonable number of genes and high BUSCO score), it will have inconsistent gene IDs (from different tools and species references). To have standardized gene IDs, GeMoMa's `AnnotationFinalizer` is used. `g=` is genome FASTA, `a=` is filtered annotation GFF3 file, `p=` sets the prefix for all gene IDs, `n=false` disables adding gene IDs as extra name attributes, and `outdir=` defines the output directory.

```bash
java -jar /path/to/GeMoMa-1.9.jar CLI AnnotationFinalizer \ 
 g=/path/to/Victoria_cruziana.genome.fa \
 a=/path/to/gemoma/Vcruz_filter_01/filtered_predictions.gff \
 p=Vcruz_ outdir=/path/to/gemoma/Vcruz_final_01/ n=false
```
After this step, all genes will be named uniformly, e.g., Vcruz_00001, Vcruz_00002, etc.

**Step 5:** Extract CDS and protein sequences using Extractor

Final annotations are in GFF3 format, to generate protein-coding sequences and peptide FASTA sequences, GeMoMa's `Extractor` can be used. `g=` specifies the genome FASTA file, `a=` specifies the finalized gff3 file, `p=true` and `c=true` outputs protein FASTA and CDS FASTA files, respectively. Output will be written to `outdir=`.

```bash
java -jar /path/to/GeMoMa-1.9.jar CLI Extractor \ 
 g=/path/to/Victoria_cruziana.genome.fa a=/path/to/gemoma/Vcruz_final_01/final_annotation.gff \
 p=true c=true
```

## 4.2 Non-coding genes

### 4.2.1 Identification of TE's with EDTA

Documentation:<https://github.com/oushujun/EDTA>

EDTA (Extensive de-novo TE Annotator) is a pipeline designed to comprehensively annotate transposable elements (TEs) in a genome. `--genome` defines the path to the genome assembly in FASTA format, `--cds` path to coding sequence FASTA file, this helps EDTA to make gene regions to reduce false-positive TE annotations. `--overwrite 1` allows EDTA to overwrite existing files in the directory (use with caution), `--anno 1` runs the complete annotation workflow after TE library construction, `--sensitive 1` activates sensitive mode for better detection of complex or nested TEs (slower), `--evaluates 1` performs an evaluation of annotation quality, `--threads` allocated CPU threads for parallel processing.

```bash
 /path/to/EDTA/EDTA.pl --genome /path/to/Victoria_cruziana.genome.fa \
 --cds /path/to/Victoria_cruziana.cds.fa #optional \
 --overwrite 1 \
 --anno 1 \
 --sensitive 1 \
 --evaluate 1 \
 --threads 10
```

### 4.2.2 Identification of ncRNAs with Infernal

Documentation: <http://eddylab.org/infernal>

Infernal (INFERence of RNA ALignment) identifies non-coding RNAs using covariance models (CMs), which are statistical models representing RNA secondary structures.
Rfam database has pre-available calibrated covariance models required by Infernal's cmscan which can be downloaded from <ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz> and an additional Claninfo file <https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin>.

```bash
#Install latest covariance models
 wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
 
#Uncompress Rfam Covariance models
 gunzip Rfm.cm.gz
 
#Index (compress) the Rfam.cm file for faster searching 
 /path/to/infernal/src/cmpress /path/to/Rfam.cm.gz
 
#Run cmscan to search the genome for ncRNA matches
 /path/to/infernal/src/cmscan --nohmmonly \
 --rfam --cut_ga --fmt 2 --oclan --oskip \
 --clanin Rfam.clanin -o /path/to/my_cmscan_out \
 --tblout /path/to/my_cmscan_tblout \
 /path/to/Rfam.cm /path/to/Victoria_cruziana.genome.fa
```

`--nohmmonly` uses covariance models instead of faster but less accurate HMM-only searches, `--rfam` tailors search to Rfam-specific cinventions, automatically applying recommended cutoffs, `--cut_ga` filters hits based on gathering thresholds (GA) from Rfam to reduce false positives, `--fmt 2` specifies human-readable format, `--oclan` outputs RNA family clan information, `--oskip` sking saving the full alignment for faster execution, `--clanin` supplies clan information file, `-o` is the output file for full cmscan report, `--tblout` provides tabular summary output file, `/path/to/Rfam.cm` defines input covariance models, and `/path/to/genome.fa` defines input genome to search for ncRNAs

Rfam database is a comprehensive database including sequence datasets of most non-coding RNA families. Hence, different types of ncRNAs can be identified with Infernal search using Rfam covariance model. However, it is also possible to use specialized tools for particular ncRNA's. Few of them are: 

- tRNAs:

    (a) tRNAscan-SE: <https://github.com/UCSC-LoweLab/tRNAscan-SE> 
    
- rRNAs: 

    (a) RNAmmer: <https://services.healthtech.dtu.dk/services/RNAmmer-1.2/> 
    (b) SSU-ALIGN: <http://eddylab.org/software/ssu-align/>

## 5. Functional Annotation

Functional annotation involves predicting the biological role of each protein-coding gene. This includes assigning functional terms (GO, KEGG), identifying homologs, and classifying proteins into families or pathways.

### 5.1 Orthlogy-based annotation via Reciprocal Best hits (RBH)

Reciprocal Best Hit (RBH) searches provide a reliable method to predict orthologous relationships by identifying gene pairs that are each other’s best match. To perform RBH search and transfer functional annotation, it is necessary to have a well-annotated reference. For this, for e.g., Arabidopsis thaliana can be used. The python script `identify_RBHs.py` at <https://github.com/bpucker/Nd1_PacBio/blob/master/identify_RBHs.py> can be used where `--prefix` defines the directory where RBH results will be saved, `--input1` is the peptide file of target species, `--input2` is the peptide file of reference species, and `--seqtype` is the type of sequence (`prot` for proteins and `nucl` for nucleotide sequences)

```bash
 python3 identify_RBHs.py \
 --prefix /path/to/output_directory/ \
 --input1 /path/to/Victoria_cruziana.pep.fasta \
 --input2 /path/to/Arabidopsis_thaliana.pep.fasta \
 --seqtype prot #'nucl' can be given instead
```

### 5.2 InterProScan

Documentation: <https://interproscan-docs.readthedocs.io/en/latest/>

InterProScan integrates multiple protein signature databases (Pfam, PRINTS, PROSITE, SMART, etc.) to classify proteins, predict domains, and important sites.
To run InterProScan, first obtain a copy of the tool, unpack it, setup and then finally run InterProScan. 

```bash
 mkdir interproscan #at your desired location
 cd interproscan
 wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/\
 interproscan-5.75-106.0-64-bit.tar.gz
 wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/\
 interproscan-5.75-106.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
 md5sum -c interproscan-5.75-106.0-64-bit.tar.gz.md5
# Must return *interproscan-5.75-106.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.

 tar -pxvzf interprocan-5.75-106.0-*-bit.tar.gz
 cd interprocan-5.75*

#to prepare hmm models into a format used by hmmscan
 python3 setup.py -f interprocan.properties 

#Run InterProScan
 ./interproscan.sh -i /path/to/Victoria_cruziana.pep.fasta \
 -f xml  -cpu 10 \
 -o /path/to/iprscan_results.xml
```

`-i` defines input protein FASTA file, `-o` defines the output file, `-cpu` allocates CPUs, and `-f` defines the output format and can be `xml`, `tsv`, `gff3`, and `json`. If the results are needed for `funannotate annotate` command, `xml` format is necessary.

### 5.3 Specialized functional annotation tools
There are some functional annotation tools for specific pathway or enzyme-family gene annotation such as:

- KIPEs3 (Flavonoid/Carotenoid biosynthesis): <https://github.com/bpucker/KIPEs> 

- MYB annotator (MYB transcription factors): <https://github.com/bpucker/MYB_annotator> 

- bHLH annotator (bHLH transcription factors): <https://github.com/bpucker/bHLH_annotator>

### 5.4 Funannotate Annotate

Documentation: <https://funannotate.readthedocs.io/en/latest/index.html>

Funannotate integrates multiple annotation sources, formats output for NCBI submission, and assigns functional terms. `--gff3` specifies final structural annotation (gff3 format), `--fasta` specifies genome assembly in FASTA format, `-s` specifies species name, `-o` for output directory, and `--sbt` is optional and used to define the path to the NCBI submission template. This is necessary if NCBI submission is desired, the file can be downloaded from NCBI website after filling basic details about the genome sequencing project (<https://submit.ncbi.nlm.nih.gov/genbank/template/submission/>). `-a` is an optional custom annotation TSV (gene ID, product name), `--iprscan` is the InterProScan XML file for protein domain annotation, and `--rename` is the LOCUS_TAG that can be specified when starting a WGS submission or is assigned automatically later. 

```bash
 funannotate annotate \
 --gff /path/to/Victoria_cruziana.gff3 \
 --fasta /path/to/Victoria_cruziana.genome.fasta \
 -s "Victoria cruziana" \
 -o /path/to/funannotate_annotate_output/ \
 --sbt /path/to/NCBI_template/ \
 -a /path/to/custom_annotations \
 --iprscan /path/to/iprscan_results.xml 
 --rename LOCUS_TAG \
 --busco_db embryophyta --cpus 10 
```

---

## 6. Data Submission

After completing genome assembly and annotation, the final and essential step is to share your data with the scientific community by submitting it to public databases. Depending on whether you have only a genome assembly or also annotation (gene predictions), the submission formats and requirements vary slightly for ENA and NCBI (INSDC databases).

Before submission, it is crucial to validate and clean the annotation files (typically in GFF3 format), as submission portals like ENA and NCBI have strict validation criteria. Common issues include duplicated feature locations and formatting inconsistencies. We recommend using the AGAT toolkit for this purpose:

```bash
agat_convert_sp_gxf2gxf.pl \
 -g /path/to/Victoria_cruziana.gff3 \
 -o /path/to/Victoria_cruziana_standard.gff3
```

`-g` defines input GFF3 annotation file and `o` is the output standardized GFF3 file

Then, fix any duplicated feature locations:
```bash
agat_sp_fix_features_locations_duplicated.pl \
 -f /path/to/Victoria_cruziana_standard.gff3 \
 -o /path/to/Victoria_cruziana_standard_deduplicated.gff3
```

`-f` is the input standardized GFF3 file and `-o` is the output GFF3 file with duplicated feature locations removed.

### 6.1 Submission to EMBL-EBI's ENA

Documentation: <https://ena-docs.readthedocs.io/en/latest/submit/assembly/genome.html>

**Step 1:** Convert to flat file format

For submitting genome assembly along with annotations to ENA, you must prepare a “flat file” format (*.embl file), which combines the genome sequence and annotation in a single file.
We recommend using the tool `EMBLmyGFF3` (<https://github.com/NBISweden/EMBLmyGFF3>). It takes as positional arguments the `gff` file and the `FASTA` file as well as some metadata and writes the output to a specified file (`-o`).

```bash
EMBLmyGFF3 /path/to/Victoria_cruziana_standard.gff3 /path/to/Victoria_cruziana.genome.fa\
 --topology linear --molecule_type 'genomic DNA' --transl_table 1 --species 'Victoria cruziana'\
  -o /path/to/Vcruz.embl 
```

**Step 2:** Compress the flat file

ENA requires gzipped flat files for submission. For efficient compression using multiple CPU threads, `pigz` (<https://github.com/madler/pigz>) tool can be used. The -k option is optional and keeps the original file.

```bash
pigz -k /path/to/Vcruz.embl
```

**Step 3:** Prepare the manifest file

The manifest file is a simple text file that describes your submission metadata (assembly name, organism, sample accession, etc.) and the path to the gzipped flatfile. The ENA documentation provides template examples.

**Step 4:** Validate submission using Webin-CLI

The `webin-CLI` tool can be used to validate your files before submission using the following command where `-validate` run validation without submission, `context genome` specifies that this is a genome submission, `-manifest` species the path to the manifest file, `-outputdir` specifies the directory where validation reports will be saved, and `userName` and `password` are the ENA login credentials. 

```bash
java -jar webin-cli-8.2.0.jar -validate -context genome -manifest /path/to/manifest_01.txt
 -outputdir /path/to/webin/01 -userName Webin-0000 -password ena-webin-password
```

If validation passes without errors, you can remove -validate to proceed with submission. Errors will be detailed in the .report files in the output directory.

### 6.2 Submission to NCBI

Documentation: <https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/>

For submission to NCBI, two main options are available depending on your dataset.

Option 1: Using Funannotate (Recommended if using functional annotations)

If you have already used Funannotate annotate, the command produces submission-ready .tbl and .sqn files compatible with NCBI requirements (see section 5.4). You can submit these files directly via the NCBI Genome Submission Portal.

Option 2: Using GAG + table2asn

If you prefer not to use Funannotate or only have structural annotation (GFF3), you can use the GAG toolkit to generate submission-ready files (<https://genomeannotation.github.io/GAG/>).

**Step 1:** Create .TBL file using GAG

Basic command to use GAG, which produces a .tbl file in the output folder. `--fasta` specifies path to genome assembly file, `--gff` specifies path to cleaned annotation file, `--out` specifies output directory for GAG results.

```bash
python gag.py --fasta Victoria_cruziana.genome.fasta \
 --gff Victoria_cruziana_standard.gff3 \
 --out gag_output
```

**Step 2** Convert .TBL to .SQN using table2asn

Table2asn (<https://www.ncbi.nlm.nih.gov/genbank/table2asn/>) converts the .tbl and genome FASTA into an NCBI submission-ready .sqn format.

Download table2asn:

```bash
# to get a local copy of table2asn 
wget https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz
gunzip linux64.table2asn.gz
mv linux64.table2asn table2asn
chmod +x table2asn
```

table2asn can be run using the following command where `-i` is the input directory with FASTA and TBL files, `-o` is the output directory, `-t` is the submission template (SBT file from NCBI), `-Z` generates error reports, and `-euk` specifies eukaryotic genome.

```bash
./table2asn -i /path/to/gag_output/ \
 -o /path/to/table2asn_output/ \
 -t SBT_template.txt \
 -M n \
 -Z -euk
```

In the output folder, `.stats` is the summary statistics file, `.val` is the detailed file of errors, and `.dr` is the discrepancy report which shows critical errors. Based on the errors encountered, GAG includes a number of options to remove questionable features like removing terminal N's, fixing start and stop codons, removing short introns etc. Please refer to the official documentation for the full details. 
Fix errors by adjusting GFF3/FASTA files, rerun GAG and table2asn until no major errors remain.

**Step 3:** Final submission

Once you obtain a valid .sqn file, submission is done through the NCBI Genome Submission Portal using your NCBI account (<https://submit.ncbi.nlm.nih.gov/subs/genome/>).

---
