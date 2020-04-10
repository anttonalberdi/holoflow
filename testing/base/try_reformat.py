#!/usr/bin/env python
input = "/home/projects/ku-cbd/people/nurher/chick_holoflow_test/assembler_test/05-Assembly/CA16_13F1b/CA16_13F1b.assembly.fa"
output= "/home/projects/ku-cbd/people/nurher/chick_holoflow_test/assembler_test/05-Assembly/CA16_13F1b/CA16_13F1b.reformat.assembly.fa"

with open(str(input)) as f_input, open(str(output), 'w') as f_output:
    seq = ''
    contig_n = 0
    contig_len_dict = {}

    for line in f_input:
        if line.startswith('>'):

            if seq:
                if len(seq) > 1000:
                    contig_n += 1
                    contig_id = ("C_"+str(contig_n))
                    seq += ('\n')

                    f_output.write(contig_id + '\n' + seq)
                    contig_len_dict[contig_id] = len(seq)

                    seq = ''

                else:
                    seq = ''


        else:
            seq += line.strip()


    if seq:
        if len(seq) > 1000:
            contig_n += 1
            contig_id = ("C_"+str(contig_n))
            seq += ('\n')
            f_output.write(contig_id + '\n' + seq)
            contig_len_dict[contig_id] = len(seq)

        else:
            pass




print(contig_len_dict)
