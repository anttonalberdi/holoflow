 # 30.06.20

rule get_paths:
    input:
        holopath=expand("{holopath}", holopath=config['holopath']),
        logpath=expand("{logpath}", logpath=config['logpath'])


################################################################################################################
 ############################################       COASSEMBLY     ############################################
################################################################################################################



##
# Assembly
## Coassembly is generated either with megahit or metaspades, chosen in config file. Megahit handles better big datasets.
rule assembly:
    input:
        read1="{projectpath}/MCB_00-MergedData/{group}_1.fastq",
        read2="{projectpath}/MCB_00-MergedData/{group}_2.fastq"
    output:
        "{projectpath}/MCB_01-Assembly/{group}_file_to_remove"
    params:
        coassembly=expand("{coassembly}", coassembly=config['coassembly']),
        klist_megahit=expand("{klist_megahit}", klist_megahit=config['klist_megahit']),
        klist_spades=expand("{klist_spades}", klist_spades=config['klist_spades']),
        threads=expand("{threads}", threads=config['threads']),
        assembler=expand("{assembler}", assembler=config['assembler']),
        out_dir="{projectpath}/MCB_01-Assembly/{group}_assembly",
        temp_assembly="{projectpath}/MCB_01-Assembly/{group}_assembly/temp_assembly.fa",
        memory=expand("{memory}", memory=config['memory']),
        group="{group}"

    shell:
        """
        python {rules.get_paths.input.holopath}/bin/holo-assembly.py \
        -1 {input.read1} -2 {input.read2} \
        -a {params.assembler} \
        -coa {params.coassembly} \
        -m {params.memory} \
        -t {params.threads} \
        -k_megahit {params.klist_megahit} \
        -k_spades {params.klist_spades} \
        -o {params.out_dir} \
        -empty_o {output} \
        -temp_a {params.temp_assembly} \
        -ID {params.group} \
        -log {rules.get_paths.input.logpath}
        """

##
# Assembly reformat
##Contigs shorter than specified in min contig len parameter of config will be removed and the contigs will be renamed
rule assembly_reformat:
    input:
        empt_file="{projectpath}/MCB_01-Assembly/{group}_file_to_remove"
    output:
        stats="{projectpath}/MCB_01-Assembly/{group}.stats",
        out_assembly="{projectpath}/MCB_01-Assembly/{group}.fa"
    params:
        group="{group}",
        stats_in="{projectpath}/PPR_03-MappedToReference/{group}.stats",
        min_contig_len=expand("{min_contig_len}", min_contig_len=config['min_contig_len']),
        in_assembly="{projectpath}/MCB_01-Assembly/{group}_assembly/temp_assembly.fa"
    shell:
        """
        rm {input.empt_file} && \
        python {rules.get_paths.input.holopath}/bin/holo-assembly_reformat.py \
        -ID {params.group} \
        -min_cl {params.min_contig_len} \
        -in_a {params.in_assembly} \
        -out_a {output.out_assembly} \
        -st_in {params.stats_in} \
        -st_out {output.stats} \
        -log {rules.get_paths.input.logpath}
        """


##
# Index assembly
## Coassembly is indexed with samtools and bwa
rule assembly_index:
    input:
        "{projectpath}/MCB_01-Assembly/{group}.fa"
    output:
        bt2_index="{projectpath}/MCB_01-Assembly/{group}.fa.rev.2.bt2l",
        samtools="{projectpath}/MCB_01-Assembly/{group}.fa.fai"
        # bwa_bwt="{projectpath}/MCB_01-Assembly/{group}.fa.bwt",
        # bwa_pac="{projectpath}/MCB_01-Assembly/{group}.fa.pac",
        # bwa_ann="{projectpath}/MCB_01-Assembly/{group}.fa.ann",
        # bwa_amb="{projectpath}/MCB_01-Assembly/{group}.fa.amb",
        # bwa_sa="{projectpath}/MCB_01-Assembly/{group}.fa.sa"
    threads: 40
    params:
        group="{group}"
    shell:
        """
        python {rules.get_paths.input.holopath}/bin/holo-assembly_index.py \
        -a {input} \
        -ia {output.samtools} \
        -bt2i {output.bt2_index} \
        -log {rules.get_paths.input.logpath} \
        -ID {params.group}
        """

##
# Assembly mapping
## map metagenomic reads to coassembly file to obtain differential coverage in next rule

rule assembly_mapping:
    input:
        bt2_index="{projectpath}/MCB_01-Assembly/{group}.fa.rev.2.bt2l",
        assembly="{projectpath}/MCB_01-Assembly/{group}.fa",
    output:
        directory("{projectpath}/MCB_02-AssemblyMapping/{group}")
    params:
        threads=expand("{threads}", threads=config['threads']),
        outdir="{projectpath}/MCB_02-AssemblyMapping/{group}",
        fq_path="{projectpath}/PPR_03-MappedToReference/{group}",
        group="{group}"
    shell:
        """
        python {rules.get_paths.input.holopath}/bin/holo-coassembly_mapping.py \
        -a {input.assembly} \
        -fq_path {params.fq_path} \
        -t {params.threads} \
        -obam_b {params.outdir} \
        -ID {params.group} \
        -log {rules.get_paths.input.logpath}
        """

# Unmapped reads
##  Write reads that did not map to assembly to files (for use in viral/diet analyses)

# ##
# # Prodigal ORF prediction
# ##
# #"Metagenomes - The simplest approach for metagenomes is to put all the sequences in one FASTA file and analyze them in Anonymous Mode."
# rule protein_prediction_prodigal:
#     input:
#         assembly="{projectpath}/MCB_01-Assembly/{group}.fa",
#         mapped_bams="{projectpath}/MCB_02-AssemblyMapping/{group}" # not necessary
#     output:
#         genetic_coords="{projectpath}/MCB_02-ProdigalPrediction/{group}.coords.gbk",
#         protein_translations="{projectpath}/MCB_02-ProdigalPrediction/{group}.protein_translations.faa"
#     params:
#         group="{group}"
#     shell: # Prodigal is run in "anon", Anonymous workflow
#         """
#         python {rules.get_paths.input.holopath}/bin/holo-pp_prodigal.py -i {input.assembly} -o {output.genetic_coords} -a {output.protein_translations} -ID {params.group} -log {rules.get_paths.input.logpath}
#         """

##
# Create depth table
## Create depth table from bam files
#
# rule depth_table:
#     input:
#         #genetic_coords="{projectpath}/MCB_02-ProdigalPrediction/{group}.coords.gbk",  #not actually necessary here, but used to keep order
#         mapped_bams="{projectpath}/MCB_02-AssemblyMapping/{group}"
#     output:
#         metabat_depth_file="{projectpath}/MCB_03-Binning/{group}_metabat/{group}.depth.txt",
#         maxbin_depth_file="{projectpath}/MCB_03-Binning/{group}_maxbin/{group}.depth.txt",
#         concoct_depth_file="{projectpath}/MCB_03-Binning/{group}_concoct/{group}.depth.txt"
#     params:
#         group="{group}",
#     shell:
#         """
#         python {rules.get_paths.input.holopath}/bin/holo-depth_files_coa.py \
#         -bam_p {input.mapped_bams} \
#         -cct {output.concoct_depth_file} \
#         -mtb {output.metabat_depth_file} \
#         -mxb {output.maxbin_depth_file} \
#         -ID {params.group} \
#         -log {rules.get_paths.input.logpath}
# #         """
# rule prepare_metawrap:
#     input:
#         directory("{projectpath}/MCB_02-AssemblyMapping/{group}")
#     output:
#         "{projectpath}/MCB_02-AssemblyMapping/{group}/dummy_fq_checkpoint.txt"
#     message: "Creating false fastq files to trick the metaWRAP binning module"
#     params:
#         bamdir="{projectpath}/MCB_02-AssemblyMapping/{group}",
#         dr1="${bam/.bam/_1.fastq}",
#         dr2="${bam/.bam/_2.fastq}"
#     shell:
#         """
#         for bam in {params.bamdir}/*.bam; do echo "@" > {params.dr1}; done
#         for bam in {params.bamdir}/*.bam; do echo "@" > {params.dr2}; done
#         touch {output}
#         """

rule metaWRAP_binning:
    input:
        "{projectpath}/MCB_02-AssemblyMapping/{group}"
    output:
        concoct=directory("{projectpath}/MCB_03-Binning/{group}/concoct_bins"),
        maxbin2=directory("{projectpath}/MCB_03-Binning/{group}/maxbin2_bins"),
        metabat2=directory("{projectpath}/MCB_03-Binning/{group}/metabat2_bins")
    params:
        assembly="{projectpath}/MCB_01-Assembly/{group}.fa",
        readfolder="{projectpath}/MCB_02-AssemblyMapping/{group}",
        threads=expand("{threads}", threads=config['threads']),
        memory=expand("{memory}", memory=config['memory']),
        outdir="{projectpath}/MCB_03-Binning/{group}",

    shell:
        """
        # Create dummy fastq/assembly files to trick metaWRAP into running without mapping
        mkdir {params.outdir}/work_files

        touch {params.outdir}/work_files/assembly.fa.bwt

        for bam in {input}/*.bam; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_1.fastq}}); done
        for bam in {input}/*.bam; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_2.fastq}}); done

        #Symlink BAMs for metaWRAP
        for bam in {input}/*.bam; do ln -s $bam {params.outdir}/work_files/$(basename $bam); done

        # Run metaWRAP binning
        module load metawrap-mg/1.2 && \
        metawrap binning -o {params.outdir} \
        -t {params.threads} \
        -m {params.memory} \
        -a {params.assembly} \
        --metabat2 \
        --maxbin2 \
        --concoct \
        {params.outdir}/work_files/*_1.fastq {params.outdir}/work_files/*_2.fastq
        """

##
# Binning with metabat
##

# rule binning_metabat:
#     input:
#         assembly="{projectpath}/MCB_01-Assembly/{group}.fa",
#         depth_table="{projectpath}/MCB_03-Binning/{group}_metabat/{group}.depth.txt"
#     output:
#         check_mtb="{projectpath}/MCB_03-Binning/{group}_metabat/{group}.mtb_checked_bins"
#     params:
#         base_mtb="{projectpath}/MCB_03-Binning/{group}_metabat/{group}.mtb",
#         bin_table_mtb="{projectpath}/MCB_03-Binning/{group}.bins_metabat.txt",
#         threads=expand("{threads}", threads=config['threads']),
#         group="{group}"
#     shell:
#         """
#         python {rules.get_paths.input.holopath}/bin/holo-binning_metabat.py \
#         -a {input.assembly} \
#         -d {input.depth_table} \
#         -bt {params.bin_table_mtb} \
#         -bb {params.base_mtb} \
#         -t {params.threads} \
#         -ID {params.group} \
#         -log {rules.get_paths.input.logpath}
#         """
#
#
# ##
# # Binning with maxbin
# ##
#
# rule binning_maxbin:
#     input:
#         assembly="{projectpath}/MCB_01-Assembly/{group}.fa",
#         depth_table="{projectpath}/MCB_03-Binning/{group}_maxbin/{group}.depth.txt"
#     output:
#         check_mxb="{projectpath}/MCB_03-Binning/{group}_maxbin/{group}.mxb_checked_bins"
#     params:
#         base_mxb="{projectpath}/MCB_03-Binning/{group}_maxbin/{group}.mxb",
#         bin_table_mxb="{projectpath}/MCB_03-Binning/{group}.bins_maxbin.txt",
#         threads=expand("{threads}", threads=config['threads']),
#         group="{group}"
#     shell:
#         """
#         python {rules.get_paths.input.holopath}/bin/holo-binning_maxbin.py \
#         -a {input.assembly} \
#         -d {input.depth_table} \
#         -bt {params.bin_table_mxb} \
#         -bb {params.base_mxb} \
#         -t {params.threads} \
#         -ID {params.group} \
#         -log {rules.get_paths.input.logpath}
#         """
#
#
# ##
# # Binning with Concoct
# ##
#
# rule binning_concoct:
#     input:
#         assembly="{projectpath}/MCB_01-Assembly/{group}.fa",
#         depth_table="{projectpath}/MCB_03-Binning/{group}_concoct/{group}.depth.txt"
#     output:
#         check_cct="{projectpath}/MCB_03-Binning/{group}_concoct/{group}.cct_checked_bins"
#     params:
#         base_cct="{projectpath}/MCB_03-Binning/{group}_concoct/{group}.cct",
#         bin_table_cct="{projectpath}/MCB_03-Binning/{group}.bins_concoct.txt",
#         min_cl_tobin=expand("{min_cl_tobin}", min_cl_tobin=config['min_cl_tobin']),
#         min_rl_tobin=expand("{min_rl_tobin}", min_rl_tobin=config['min_rl_tobin']),
#         threads=expand("{threads}", threads=config['threads']),
#         group="{group}"
#     shell:
#         """
#         python {rules.get_paths.input.holopath}/bin/holo-binning_concoct.py \
#         -a {input.assembly} \
#         -d {input.depth_table} \
#         -bt {params.bin_table_cct} \
#         -bb {params.base_cct} \
#         -l {params.min_cl_tobin} \
#         -r {params.min_rl_tobin} \
#         -t {params.threads} \
#         -ID {params.group} \
#         -log {rules.get_paths.input.logpath}
#         """
##
# Binning with vamb
##

# rule binning_vamb:
#     input:
#         assembly="{projectpath}/MCB_01-Assembly/{group}.fa",
#         depth_table="{projectpath}/MCB_03-Binning/{group}_metabat/{group}.depth.txt"
#     output:
#         check_vamb="{projectpath}/MCB_03-Binning/{group}_vamb/{group}.vmb_checked_bins"
#     params:
#         base_vmb="{projectpath}/MCB_03-Binning/{group}_vamb/",
#         bin_table_vmb="{projectpath}/MCB_03-Binning/{group}.bins_vamb.txt",
#         group="{group}"
#     shell:
#         """
#         python {rules.get_paths.input.holopath}/bin/holo-binning_vamb.py \
#         -a {input.assembly} \
#         -d {input.depth_table} \
#         -bt {params.bin_table_vmb} \
#         -bb {params.base_vmb} \
#         -ID {params.group} \
#         -log {rules.get_paths.input.logpath}
#         """

##
# Check binning
## If all binners created bins, then continue. If Any binner did not create bins, simply copy resulting bins from
# # another software, rename these, and continue (the result is going to be the same for DAStool, is for the sake of the pipeline)
# rule check_bins:
#     input:
#         check_mxb="{projectpath}/MCB_03-Binning/{group}_maxbin/{group}.mxb_checked_bins",
#         check_mtb="{projectpath}/MCB_03-Binning/{group}_metabat/{group}.mtb_checked_bins",
#         check_cct="{projectpath}/MCB_03-Binning/{group}_concoct/{group}.cct_checked_bins",
# #        check_vmb="{projectpath}/MCB_03-Binning/{group}_vamb/{group}.vmb_checked_bins"
#     output:
#         "{projectpath}/MCB_03-Binning/{group}_checked_bins.txt",
#     params:
#         binning_dir="{projectpath}/MCB_03-Binning",
#         group="{group}"
#     shell:
#         """
#         python {rules.get_paths.input.holopath}/bin/holo-check_bins.py \
#         --check_cct {input.check_cct} \
#         -check_mtb {input.check_mtb} \
#         -check_mxb {input.check_mxb} \
#         -binning_dir {params.binning_dir} \
#         -ID {params.group} \
#         -log {rules.get_paths.input.logpath}
#         """
# #        --check_vmb {input.check_vmb} \

##
# Bin refinement with DASTool using binning: metabat, maxbin, concoct and vamb
##

# rule das_tool:
#     input:
#         checked_bins="{projectpath}/MCB_03-Binning/{group}_checked_bins.txt",
#         assembly="{projectpath}/MCB_01-Assembly/{group}.fa"#,
#         #pproteins="{projectpath}/MCB_02-ProdigalPrediction/{group}.protein_translations.faa"
#     output:
#         directory("{projectpath}/MCB_04-BinMerging/{group}_files")
#     params:
#         threads=expand("{threads}", threads=config['threads']),
#         search_eng=expand("{search_eng}", search_eng=config['search_eng']),
#         bin_table_mxb="{projectpath}/MCB_03-Binning/{group}.bins_maxbin.txt",
#         bin_table_mtb="{projectpath}/MCB_03-Binning/{group}.bins_metabat.txt",
#         bin_table_cct="{projectpath}/MCB_03-Binning/{group}.bins_concoct.txt",
# #        bin_table_vmb="{projectpath}/MCB_03-Binning/{group}.bins_vamb.txt",
#         dastool_db=expand("{dastool_db}", dastool_db=config['dastool_db']),
#         dastool_dir="{projectpath}/MCB_04-BinMerging/{group}",
#         group="{group}"
#     shell:
#         """
#         python {rules.get_paths.input.holopath}/bin/holo-binning_dastool.py \
#         -cb {input.checked_bins} \
#         -a {input.assembly}  \
#         --bt_cct {params.bin_table_cct} \
#         -bt_mtb {params.bin_table_mtb} \
#         -bt_mxb {params.bin_table_mxb} \
#         -o {params.dastool_dir} \
#         -se {params.search_eng} \
#         -t {params.threads} \
#         -db {params.dastool_db} \
#         -ID {params.group} \
#         -log {rules.get_paths.input.logpath}
#         """
        #python {rules.get_paths.input.holopath}/bin/holo-binning_dastool.py -cb {input.checked_bins} -a {input.assembly} --bt_cct {params.bin_table_cct} -bt_mtb {params.bin_table_mtb} -bt_mxb {params.bin_table_mxb} -p {input.pproteins} -o {params.dastool_dir} -se {params.search_eng} -t {params.threads} -db {params.dastool_db} -ID {params.group} -log {rules.get_paths.input.logpath}
#        --bt_vmb {params.bin_table_vmb} \

rule metawrap_refinement:
    input:
        concoct="{projectpath}/MCB_03-Binning/{group}/concoct_bins",
        maxbin2="{projectpath}/MCB_03-Binning/{group}/maxbin2_bins",
        metabat2="{projectpath}/MCB_03-Binning/{group}/metabat2_bins",
    output:
        stats="{projectpath}/MCB_04-BinMerging/{group}_files/metawrap_70_10_bins.stats",
        workdir="{projectpath}/MCB_04-BinMerging/{group}_files",
    params:
        threads=expand("{threads}", threads=config['threads']),
        memory=expand("{memory}", memory=config['memory']),
        group="{group}"
    shell:
        """
        module load metawrap-mg/1.2 && \
        metawrap bin_refinement \
            -m {params.memory} \
            -t {params.threads} \
            -o {output.workdir} \
            -t {params.threads} \
            -A {input.concoct} \
            -B {input.maxbin2} \
            -C {input.metabat2} \
            -c 70 \
            -x 10
        # Rename metawrap bins to match coassembly group:
        sed -i'' '2,$s/bin/bin_{params.group}/g' {output.stats}
        """


rule coverm:
    input:
        "{projectpath}/MCB_04-BinMerging/{group}_files/metawrap_70_10_bins.stats"
    output:
        coverm="{projectpath}/MCB_05-CoverM/{group}_files/{group}_coverM.txt",
        covermdir=directory("{projectpath}/MCB_05-CoverM/{group}_files")
    params:
        all_mw="{projectpath}/MCB_04-BinMerging/All_files",
        groups="{projectpath}/MCB_04-BinMerging/{group}",
        assembly="{projectpath}/MCB_01-Assembly/{group}.fa",
        mapped_bams="{projectpath}/MCB_02-AssemblyMapping/{group}",
        threads=expand("{threads}", threads=config['threads']),
        group="{group}"
    shell:
        """
        echo {input}
        module load tools coverm/0.6.1 && \
        coverm genome \
            -b {params.mapped_bams}/*.bam \
            --genome-fasta-files {params.assembly} \
            -m relative_abundance \
            -t {params.threads} \
            --min-covered-fraction 0 \
            > {output.coverm}

        #Merge group metaWRAP refinement results
        mkdir {params.all_mw}

        #setup headers for combined metawrap file:
        echo -e bin' \t 'completeness' \t 'contamination' \t 'GC' \t 'lineage' \t 'N50' \t 'size' \t 'binner > header.txt

        #Cat the bin info from each group together
        for group in {params.groups}; \
            do grep -v 'contamination' "$group"/metawrap_70_10_bins.stats >> bins.txt; done

        #Merge header with bins:
        cat header.txt bins.txt > {params.all_mw}/All_metawrap_70_10_bins.stats

        #Copy bins from each group to a new folder in the 'All_files' directory
        mkdir {params.all_mw}/All_metawrap_70_10_bins

        for group in {params.groups}; \
            do for bin in "$group"/metawrap_70_10_bins/*.fa; \
                do cp $bin {params.all_mw}/All_metawrap_70_10_bins/$(basename ${{bin/bin./"${{group/_files/}}"_bin.}}); \
                done; \
                    done

        #Clean up
        rm header.txt
        rm bins.txt
        """


#This rule merges metawrap .stats files from multiple groups for use in dereplication
#It also renames and copies the bins from each group to the 'All_files' folder
# onsuccess:
#     input:
#         "{projectpath}/MCB_04-BinMerging/{group}_files/metawrap_70_10_bins.stats"
#     output:
#         "{projectpath}/MCB_04-BinMerging/All_files/metawrap_70_10_bins.stats"
#     params:
#         wd="{projectpath}/MCB_04-BinMerging"
#     shell:
#         """
#         mkdir "{params.wd}/All_files/
#
#         #setup headers for combined metawrap file:
#         echo -e bin' \t 'completeness' \t 'contamination' \t 'GC' \t 'lineage' \t 'N50' \t 'size' \t 'binner > header.txt
#
#         #Cat the bin info from each group together
#         for group in {params.wd}/*_files; \
#             do grep -v 'contamination' "$group"/metawrap_70_10_bins.stats >> bins.txt; done
#
#         #Merge header with bins:
#         cat header.txt bins.txt > {output}
#
#         #Copy bins from each group to a new folder in the 'All_files' directory
#         mkdir {params.wd}/All_files/metawrap_70_10_bins
#
#         for group in {params.wd}/*_files; \
#             do for bin in "$group"/metawrap_70_10_bins/*.fa; \
#                 do echo cp $bin {input}/metawrap_70_10_bins/$(basename ${{bin/bin./"${{group/_files/}}"_bin.}}); \
#                 done; \
#                     done
#
#         #Clean up
#         rm header.txt
#         rm bins.txt
#         """


onsuccess:
    print("Job success!")
    shell("""mail -s "workflow completed" raph.eisenhofer@gmail.com < {log}""")

onerror:
    print("An error occurred")
    shell("""mail -s "an error occurred" raph.eisenhofer@gmail.com < {log}""")



# rule singleM:
#     input:
#         "{projectpath}/MCB_04-BinMerging/All_files/metawrap_50_10_bins.stats"
#     output:
#         "{projectpath}/MCB_05-SingleM/{group}/<><><><>")
#     params:
#         threads=expand("{threads}", threads=config['threads']),
#         group="{group}"
#     conda:
#         "{holopath}/workflows/metagenomics/coassembly_binning/conda.yaml"
#     shell:
#         """
#
#         """


##
# RefineM bin refinement
##
#>refinem filter_bins <bin_dir> <outlier_output_dir>/outliers.tsv <filtered_output_dir>
# rule bin_refinement:
#     input:
#         assembly="{projectpath}/MCB_01-Assembly/{group}.fa",
#         assembly_map="{projectpath}/MCB_02-AssemblyMapping/{group}.mapped.bam",
#         check_dastool="{projectpath}/MCB_04-BinMerging/{group}_DASTool_bins"
#     output:
#         directory("{projectpath}/MCB_05-BinRefinement/{group}")
#     params:
#         dastool_bin_dir="{projectpath}/MCB_04-BinMerging/{group}_DASTool_bins",
#         threads=expand("{threads}", threads=config['threads']),
#         group="{group}"
#     shell:
#         """
#         python {rules.get_paths.input.holopath}/bin/holo-bin_refinement.py -a {input.assembly} -bam {input.assembly_map} -dastool_bd {params.dastool_bin_dir} -out_dir {output} -ID {params.group} -t {params.threads} -log {rules.get_paths.input.logpath}
#         """
