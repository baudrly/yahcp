# yahcp (Yet Another Hi-C Pipeline) 
    
Usage: 

    ./yahcp.sh -1 reads_forward.fastq -2 reads_reverse.fastq -f genome.fa [-s size] [-o output_directory] [-e enzyme] [-q quality_min] [--duplicates] [--clean-up]
    
Generate a sparse, [GRAAL](https://github.com/koszullab/GRAAL)-compatible contact map from paired-end reads and a reference genome.
The map can also be easily visualized with [HiC-Box](https://github.com/koszullab/HiC-Box).
Information about fragments and contigs/chromosomes are stored in separate files.
The genome can either be partitioned by restriction fragments (specifying the enzyme) or fixed size chunks (specifying a number).
    
Requires [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [samtools](http://samtools.sourceforge.net/), [bedtools](https://bedtools.readthedocs.io/en/latest/) and Python (with [Biopython](https://biopython.org/) installed) to run. Optionally, [minimap2](https://lh3.github.io/minimap2/) can be used instead of bowtie2 if specified.
    
Parameters:
    
    -1 or --forward: Forward FASTQ reads
    -2 or --reverse: Reverse FASTQ reads
    -f or --fasta: Reference genome to map against in FASTA format
    -o or --output: Output directory. Defaults to the current directory.
    -e or --enzyme: Restriction enzyme if a string, or chunk size (i.e. resolution) if a number. Defaults to 5000 bp chunks.
    -q or --quality-min: Minimum mapping quality for selecting contacts. Defaults to 30.
    -d or --duplicates: If enabled, removes adapters and PCR duplicates prior to mapping. Not enabled by default.
    -s or --size: Minimum size threshold to consider contigs. Defaults to 0 (keep all contigs).
    -c or --clean-up: If enabled, removes intermediary BED files after generating the contact map. Enabled by default.
    -t or --threads: Number of threads to use for the aligner and samtools. Defaults to 1.
    -m or --minimap: Use the minimap2 aligner instead of bowtie2. Not enabled by default.
    -T or --tmp: Directory for storing intermediary BED files and temporary sort files. Defaults to the output directory.
    -h or --help: Display this help message
    
The expected files in the output directory will take the form:

* ```abs_fragments_contacts_weighted.txt```: the sparse contact map
* ```fragments_list.txt```: information about restriction fragments (or chunks)
* ```info_contigs.txt```: information about contigs or chromosomes

Please be aware that intermediary BED files can be quite bulky. You may use the ```--tmp``` option to store all such files on a different location (e.g. a scratch partition if running the pipeline on a computer cluster).
