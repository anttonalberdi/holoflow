
Input same as metagenomics_IB but with GROUPS!
    Either all together or per groups

        Check assemblers, megahit probably can take all paths , , ,
            Input can be string
        Metaspades probably needs all files together, snakemake .py file for assembly contemplate that and merge

Input for the first rule will be a DIRECTORY PATH (the one in common for all samples in input, otherwise create one and move them there)
