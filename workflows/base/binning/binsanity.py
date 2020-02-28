# Profile command
Binsanity-profile -i fasta_file -s {sam,bam}_file -c output_file

# Only Binsanity
Binsanity -f [/path/to/fasta] -l [FastaFile] -c [CoverageFile] -o [OutputDirectory]
    # Refinement
        # --kmer KMER           Indicate a number for the kmer calculation, the [Default: 4]
        # --refine-preference INPUTREFINEDPREF
                        #Specify a preference for refinement. [Default: -25]

# Binsanity-lc - k.mers clustering, more feasible for memory errors
# Binsanity - only coverage for clustering 


# There is Binsanity -- refine (Binsanity-CheckM for complt-Binsanity refine-Clusters)

rule binning_maxbin:
    input:
        assembly_idx="{projectpath}/05-Assembly/{sample}/{sample}.reformat.assembly.fa",
        assemblybam="{projectpath}/06-Assembly_mapping/{sample}.mapped.bam"
    output:
        dir_mxb=directory("{projectpath}/07-Binning/{sample}.maxbin"),
        #depth_file_mxb="{projectpath}/07-Binning/{sample}.maxbin/{sample}.depth.txt",
        bin_table_mxb="{projectpath}/07-Binning/{sample}.bins_maxbin.txt"
    params:
        threads=expand("{threads}", threads=config['threads'])
    run:
        shell("module unload gcc && module load tools perl/5.20.2 metabat/2.12.1 && jgi_summarize_bam_contig_depths --outputDepth {output.dir_mxb}/{sample}.depth.txt --noIntraDepthVariance {input.assemblybam}")
        shell("module unload gcc && module load tools perl/5.20.2 maxbin/2.2.7 fraggenescan/1.31 && run_MaxBin.pl -contig {input.assembly_idx} -abund {output.dir_mxb}/{sample}.depth.txt -out {output.dir_mxb} -thread {params.threads}")

        #Generate bin table
        bintable=open(output.bin_table_mxb,"a+")
        binlist = glob.glob(output.dir_mxb)
        for bin in binlist:
            binname = os.path.splitext(os.path.basename(bin))[0]+''
            with open(bin, 'r') as binfile:
               for line in binfile:
                    if line.startswith('>'):
                        contig = line.strip()
                        contig = contig.replace(">", "")
                        bintable.write("{0}\t{1}\r\n".format(contig,binname))
        bintable.close()
