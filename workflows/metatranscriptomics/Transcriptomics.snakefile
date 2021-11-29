################################################################################
################################################################################
################################################################################
# Snakefile for metatranscriptomics
# Raphael Eisenhofer 29/11/2021
#
################################################################################
################################################################################
################################################################################

### Setup sample wildcard:
import os
from glob import glob

SAMPLE = [os.path.basename(fn).replace("_1.fastq.gz", "")
            for fn in glob(f"2_Reads/0_Untrimmed/*_1.fastq.gz")]

print("Detected the following samples:")
print(SAMPLE)
################################################################################
### Setup the desired outputs
rule all:
    input:

################################################################################
### Filter reads with fastp
rule fastp:
    input:
        r1i = "2_Reads/0_Untrimmed/{sample}_1.fastq.gz",
        r2i = "2_Reads/0_Untrimmed/{sample}_2.fastq.gz"
    output:
        r1o = "2_Reads/1_Trimmed/{sample}_trimmed_1.fastq.gz",
        r2o = "2_Reads/1_Trimmed/{sample}_trimmed_2.fastq.gz",
        fastp_html = "2_Reads/2_fastp_reports/{sample}.html",
        fastp_json = "2_Reads/2_fastp_reports/{sample}.json"
    conda:
        "Transcriptomics_conda.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{sample}_fastp.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_fastp.log"
    message:
        "Using FASTP to trim adapters and low quality sequences for {wildcards.sample}"
    shell:
        """
        fastp \
            --in1 {input.r1i} --in2 {input.r2i} \
            --out1 {output.r1o} --out2 {output.r2o} \
            --trim_poly_g \
            --trim_poly_x \
            --n_base_limit 5 \
            --qualified_quality_phred 20 \
            --length_required 60 \
            --thread {threads} \
            --html {output.fastp_html} \
            --json {output.fastp_json} \
            --adapter_sequence CTGTCTCTTATACACATCT \
            --adapter_sequence_r2 CTGTCTCTTATACACATCT \
        &> {log}
        """
################################################################################
## Index host genomes:
rule index_ref:
    input:
        "1_References"
    output:
        bt2_index = "1_References/CattedRefs.fna.gz.rev.2.bt2l",
        catted_ref = "1_References/CattedRefs.fna.gz"
    conda:
        "1_QC.yaml"
    threads:
        40
    log:
        "3_Outputs/0_Logs/host_genome_indexing.log"
    message:
        "Concatenating and indexing host genomes with Bowtie2"
    shell:
        """
        # Concatenate input reference genomes
        cat {input}/*.gz > {input}/CattedRefs.fna.gz

        # Index catted genomes
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {output.catted_ref} {output.catted_ref} \
        &> {log}
        """
################################################################################
### Map a sample's reads to it corresponding assembly
rule assembly_mapping:
    input:
        bt2_index = "3_Outputs/2_Assemblies/{sample}_contigs.fasta.rev.2.bt2l"
    output:
        mapped_bam = "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam"
    params:
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta",
        r1 = "2_Reads/3_Host_removed/{sample}_non_host_1.fastq.gz",
        r2 = "2_Reads/3_Host_removed/{sample}_non_host_2.fastq.gz"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_mapping.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_mapping.log"
    message:
        "Mapping {wildcards.sample} to its assembly using Bowtie2"
    shell:
        """
        # Map reads to assembly using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {params.assembly} \
            -1 {params.r1} \
            -2 {params.r2} \
        | samtools sort -@ {threads} -o {output.mapped_bam}
        """
################################################################################
### Bin each sample's contigs using metaWRAP's binning module
rule metaWRAP_binning:
    input:
        "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam"
    output:
        concoct = directory("3_Outputs/4_Binning/{sample}/concoct_bins"),
        maxbin2 = directory("3_Outputs/4_Binning/{sample}/maxbin2_bins"),
        metabat2 = directory("3_Outputs/4_Binning/{sample}/metabat2_bins")
    params:
        outdir = "3_Outputs/4_Binning/{sample}",
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta",
        basename = "3_Outputs/3_Assembly_Mapping/BAMs/{sample}",
        memory = "180"
    conda:
        "2_MetaWRAP.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_binning.log"
    message:
        "Binning {wildcards.sample} contigs with MetaWRAP (concoct, maxbin2, metabat2)"
    shell:
        """
        # Create dummy fastq/assembly files to trick metaWRAP into running without mapping
        mkdir {params.outdir}/work_files

        touch {params.outdir}/work_files/assembly.fa.bwt

        echo "@" > {params.outdir}/work_files/$(basename {params.basename}_1.fastq)
        echo "@" > {params.outdir}/work_files/$(basename {params.basename}_2.fastq)

        #Symlink BAMs for metaWRAP
        ln -s `pwd`/{input} {params.outdir}/work_files/$(basename {input})

        # Run metaWRAP binning
        metawrap binning -o {params.outdir} \
            -t {threads} \
            -m {params.memory} \
            -a {params.assembly} \
            --metabat2 \
            --maxbin2 \
            --concoct \
        {params.outdir}/work_files/*_1.fastq {params.outdir}/work_files/*_2.fastq
        """
################################################################################
### Automatically refine bins using metaWRAP's refinement module
rule metaWRAP_refinement:
    input:
        concoct = "3_Outputs/4_Binning/{sample}/concoct_bins",
        maxbin2 = "3_Outputs/4_Binning/{sample}/maxbin2_bins",
        metabat2 = "3_Outputs/4_Binning/{sample}/metabat2_bins",
    output:
        stats = "3_Outputs/5_Refined_Bins/{sample}_metawrap_70_10_bins.stats",
        contigmap = "3_Outputs/5_Refined_Bins/{sample}_metawrap_70_10_bins.contigs"
    params:
        outdir = "3_Outputs/5_Refined_Bins/{sample}",
        memory = "180",
        sample = "{sample}"
    conda:
        "2_MetaWRAP.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_bin_refinement.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_bin_refinement.log"
    message:
        "Refining {wildcards.sample} bins with MetaWRAP's bin refinement module"
    shell:
        """
        # Setup checkM path
        printf "/home/projects/ku-cbd/people/rapeis/0_DBs/CHECKM" | checkm data setRoot

        metawrap bin_refinement \
            -m {params.memory} \
            -t {threads} \
            -o {params.outdir} \
            -A {input.concoct} \
            -B {input.maxbin2} \
            -C {input.metabat2} \
            -c 70 \
            -x 10

        # Rename metawrap bins to match coassembly group:
        cp {params.outdir}/metawrap_70_10_bins.stats {output.stats}
        cp {params.outdir}/metawrap_70_10_bins.contigs {output.contigmap}
        sed -i'' '2,$s/bin/bin_{params.sample}/g' {output.stats}
        sed -i'' 's/bin/bin_{params.sample}/g' {output.contigmap}
        """
################################################################################
### Calculate the number of reads that mapped to coassemblies
rule coverM_assembly:
    input:
        stats = "3_Outputs/5_Refined_Bins/{sample}_metawrap_70_10_bins.stats",
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta",
        mapped_bam = "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam"
    output:
        "3_Outputs/6_CoverM/{sample}_coverM.txt"
    params:
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{sample}_coverM.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_coverM.log"
    message:
        "Calculating assembly mapping rate for {wildcards.sample} with CoverM"
    shell:
        """
        coverm genome \
            -b {input.mapped_bam} \
            --genome-fasta-files {input.assembly} \
            -m relative_abundance \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output}
        """
