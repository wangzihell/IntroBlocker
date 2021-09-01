#!/usr/bin/env python

import sys
from copy import deepcopy
import warnings
from optparse import OptionParser

def main_sample_sel (header, data, max_black_rate, num_main_sample):
    '''
    select main sample that minimize black block.
    m_b_r: maxmum black block ratio.
    num_sm: main sample number.
    '''
    main_sam_list = []
    black_block = (len(header)-3)*(len(data))
    max_black = black_block * max_black_rate
    
    if num_main_sample is not None:
        num_main_sample = int(num_main_sample)

    while(not black_block <= max_black):
        max_new_sam_id = 0
        for sam_id in range(3, len(header)):
            temp_black_block = 0
            if sam_id in main_sam_list:
                continue
            for row in data:
                for sam_id_tmp in range(3, len(row)):
                    add_switch = True
                    for x in main_sam_list:
                        if row[sam_id_tmp] == row[x]:
                            add_switch = False
                    if row[sam_id_tmp] == row[sam_id]:
                        add_switch = False
                    if add_switch:
                        temp_black_block += 1
            if temp_black_block < black_block:
                black_block = temp_black_block
                max_new_sam_id = sam_id
        if max_new_sam_id != 0:
            main_sam_list.append(max_new_sam_id)
            if num_main_sample is not None and len(main_sam_list) >= num_main_sample:
                break

    return main_sam_list


def reorder_main_sam (data, main_sample_list):
    # reorder main samples to minimize color switch time
    loop_switch = True
    N1_l = deepcopy(main_sample_list)
    N1 = count_N(deepcopy(data), N1_l)
    while loop_switch:
        loop_switch = False
        for i in range(0, len(main_sample_list)-1):
            N2_l = deepcopy(N1_l)
            N2_l[i], N2_l[i+1] = N2_l[i+1], N2_l[i]
            N2 = count_N(deepcopy(data), N2_l)

            if N2 < N1:
                N1 = N2
                N1_l = deepcopy(N2_l)
                loop_switch = True
    N1 = count_N(data, N1_l)

    return N1, N1_l


def reorder_main_sam_semi (data, main_sample_matrix):
    # reorder main samples **in each group** to minimize color switch time
    loop_switch = True
    N1_l = deepcopy(main_sample_matrix)
    N1 = count_N(deepcopy(data), [i[0] for i in N1_l])
    while loop_switch:
        loop_switch = False
        for i in range(0, len(main_sample_matrix)-1):

            # do not swith inter-group
            if N1_l[i][1] != N1_l[i+1][1]:
                continue

            N2_l = deepcopy(N1_l)
            N2_l[i], N2_l[i+1] = N2_l[i+1], N2_l[i]
            N2 = count_N(deepcopy(data), [i[0] for i in N2_l])

            if N2 < N1:
                N1 = N2
                N1_l = deepcopy(N2_l)
                loop_switch = True
    N1 = count_N(data, [i[0] for i in N1_l])

    return N1, [i[0] for i in N1_l]


def count_N (data, main_sample_list):
    jump_count = 0
    # assigne each sample into a group according to sample order of main_sample_list.
    for row in data:
        for main_sam_id in range(0, len(main_sample_list)):
            if row[main_sample_list[main_sam_id]] > 0:  ## this sample has not been assigned to a group yet. ## change >=0 to >0, cause now using 0 to represent CNV
                for main_sam_tmp in range(main_sam_id+1, len(main_sample_list)):
                    if row[main_sample_list[main_sam_tmp]] > 0:  ## this sample has not been assigned to a group yet. ## change >=0 to >0, cause now using 0 to represent CNV
                        if row[main_sample_list[main_sam_tmp]] == row[main_sample_list[main_sam_id]]:
                            row[main_sample_list[main_sam_tmp]] = -(main_sam_id+1)
                row[main_sample_list[main_sam_id]] = -(main_sam_id+1)
    # count jump number
    for i in range(3, len(data[0])):
        for j in range(0, len(data)-1):
            if data[j][i] != data[j+1][i]:
                jump_count += 1

    return jump_count


def change_order (data_matrix, sample_list):
    # sort sample according to their contribution
    data_matrix_flat = [item for sublist in data_matrix for item in sublist]
    exp_percent = list()

    for i in range(0, len(sample_list)):
        exp_percent.append((data_matrix_flat.count(-i-1)+data_matrix_flat.count(i+1))/(len(data_matrix)-1))

    tmp_index = sorted(range(len(exp_percent)), key=lambda k: exp_percent[k])[::-1]
    sample_list = [sample_list[i] for i in tmp_index]

    tmp_index = [-(i+1) for i in tmp_index] #desc 
    new_index = [-i for i in range(1, len(sample_list)+1)]

    rep_dict = dict(zip(tmp_index, new_index))
    for i in range(0, len(data_matrix)):
        for j in range(3, len(data_matrix[i])):
            if data_matrix[i][j] < 0:  # non-uniq-block only
                data_matrix[i][j] = rep_dict[data_matrix[i][j]]
            elif data_matrix[i][j] > 0:  # for block has been fliped
                data_matrix[i][j] = -rep_dict[-data_matrix[i][j]]

    return data_matrix, sample_list


def mark_singleton_block (data):
    # mark uniq block as positive number
    for i, x in enumerate(data):
        tmp_index = [j for j in set(x) if x.count(j)==1 and j!=0]  # 0 represent CNV, do not change it
        for j, y in enumerate(x):
            if j > 2 and y in tmp_index:
                data[i][j] = -data[i][j]

    return data


def data_read (infile, chr_col=1, s_col=4):
    try:
        if infile :
            IN = open(infile, 'r')
        else :
            IN = sys.stdin
    except IOError:
        sys.stderr.write("\n[Error]:\n\t File cannot be open: %s\n" % infile )
        exit(-1)
    #
    notcomm = False
    sample_header = []
    data_matrix = []
    for line in IN :
        if line[0:5] == "CHROM":
            notcomm = True
            sample_header = line.strip().split("\t")
            continue
        #
        if notcomm:
            tokens = line.strip().split("\t")
            tokens[3:] = [int(tokens[i]) for i in range(3,len(tokens))]
            data_matrix.append(tokens)
    #
    return sample_header, data_matrix


def report_stat (data_matrix, final_list, contri_sample_list, prior_sample_list, sample_header, outfile):
    # main file sample (.ms): sample name, sample label (the negative number), explain percent (including self), self explain, explain others, sample three number including positive number.
    # _t suffix indicate including singleton block (positive number)
    # name_list, self_percent, exp_percent, exp_other_percent, self_percent_t, exp_percent_t, exp_other_percent_t, sample_label = ([] for i in range(8))

    black_list = [0]*(len(sample_header)-3)  ## contains black block number of each sample.
    for i in data_matrix:
        for j in range(3, len(sample_header)):
            if i[j] > 0:
                black_list[j-3]+=1

    name_list, self_percent_t, exp_percent_t, exp_other_percent_t, sample_label = ([] for i in range(5))
    data_matrix_flat = [item for sublist in data_matrix for item in sublist]

    for i in range(0, len(final_list)):
        group_n = final_list[i]
        sample_mat = [data[group_n] for data in data_matrix]
        name_list.append(sample_header[group_n])
        sample_label.append(str(-i-1))
        #
        self_percent_t.append("{:.2f}".format((sample_mat.count(-i-1)+sample_mat.count(i+1))/(len(data_matrix)-1)))
        exp_percent_t.append("{:.2f}".format((data_matrix_flat.count(-i-1)+data_matrix_flat.count(i+1))/(len(data_matrix)-1)))
        exp_other_percent_t.append("{:.2f}".format((data_matrix_flat.count(-i-1)+data_matrix_flat.count(i+1))/(len(data_matrix)-1)-(sample_mat.count(-i-1)+sample_mat.count(i+1))/(len(data_matrix)-1)))        
    
    # final order file (.fi)
    final_sample_file = open(outfile+".fi", "w")
    final_sample_file.write("\n".join([sample_header[final_list[i]]+"\t"+str(-i-1) for i in range(0, len(final_list))])+"\n")

    # all sample file (.as): sample name, 
    all_sample_file = open(outfile+".as", "w")
    all_sample_file.write("\n".join([sample_header[final_list[i]]+"\t"+"{:.2f}".format(1-black_list[final_list[i]-3]/len(data_matrix)) for i in range(0, len(final_list))])+"\n")

    # priority order file (.pri)
    pri_sample_file = open(outfile+".pri", "w")
    pri_sample_file.write("\n".join([sample_header[i]+"\t"+str(-final_list.index(i)-1) for i in prior_sample_list])+"\n")

    # main sample file (.ms, the contribution file): 
    main_sample_file = open(outfile+".ms", 'w')
    main_sample_file.write("\n".join([sample_header[i]+"\t"+str(-final_list.index(i)-1)+"\t"+exp_percent_t[final_list.index(i)]+"\t"+self_percent_t[final_list.index(i)]+"\t"+exp_other_percent_t[final_list.index(i)] for i in contri_sample_list])+"\n")
    # main_sample_file.write("\n".join([name_list[i]+"\t"+sample_label[i]+"\t"+exp_percent_t[i]+"\t"+self_percent_t[i]+"\t"+exp_other_percent_t[i] for i in range(0, len(final_list))])+"\n")


# ===========================================
def main():
    usage = "Input format:\n" \
            "   chr\tstart\tend\tsample1\tsample2\tsample3\n" \
            "   chr1A\t1\t1000000\t0\t0\t1\n" \
            "run mode:\n" \
            "un-supervised: all samples are treated equally.\n" \
            "semi-supervised: samples within each group are reordered. Main sample file should be a two-column file. Group of each sample is given in the second column, sorted in priority.\n" \
            "supervised: sample order are fixed by given main sample list file, samples not given in main sample file are ordered randomly.\n" \
            "color mode:\n" \
            "contribution: default mode, sample color order are sorted by their contribution." \
            "priority:  sample color order are sorted by their priority."
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("-c", dest="chr_col",
                  help="column id for chromosome [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("-s", dest="s_col",
                  help="column id for first sample [default: %default]", metavar="INT",
                      default = 4)
    parser.add_option("-o", dest="outfile",
                  help="Output file, use STDOUT if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("-b", dest="max_black_ratio",
                  help="max quantity of black blocks ; ", default=0.05)
    parser.add_option("-n", dest="num_main_sample",
                  help="upper limit of number of main sample ; ", default=None)
    parser.add_option("-m", dest="main_sample",
                  help="main sample file", default=None)
    parser.add_option("-M", dest="run_mode",
                  help="main sample file", default="un-supervised")
    parser.add_option("-C", dest="color_mode",
                  help="color order mode", default="contribution")
    #
    (options, args) = parser.parse_args()
    # 
    def unique(sequence):
        seen = set()
        return [x for x in sequence if not (x in seen or seen.add(x))]
    #
    sample_header, data_matrix = data_read(options.infile, options.chr_col, options.s_col)
    #
    if options.run_mode == "un-supervised":
        if options.main_sample is not None:
            warnings.warn("in un-supervised mode, main_sample list is not used.")
        #
        main_sample_list = main_sample_sel(sample_header, data_matrix, options.max_black_ratio, options.num_main_sample)
        [main_sample_list.append(i) for i in range(3, len(sample_header)) if i not in main_sample_list]
        jump_count, reordered_main_sample_list = reorder_main_sam(data_matrix, main_sample_list)
        #
    elif options.run_mode == "semi-supervised":
        # semi-supervised method should go main sample selection first, for that reorder will find local optima too fast.
        if options.main_sample is None:
            raise ValueError('In semi-supervised mode, a two column main sample file must be given.')

        # first do large-scale selection
        main_sample_list = main_sample_sel(sample_header, data_matrix, options.max_black_ratio, options.num_main_sample)
        [main_sample_list.append(i) for i in range(3, len(sample_header)) if i not in main_sample_list]

        IN = open(options.main_sample, 'r')
        main_sample_matrix_tmp = []
        for line in IN :
            tokens = line.strip().split("\t") 
            main_sample_matrix_tmp.append([sample_header.index(tokens[0]), tokens[1]])
        
        group_name = unique([i for _,i in main_sample_matrix_tmp])
        main_sample_matrix_tmp = sorted(main_sample_matrix_tmp, key=lambda x: main_sample_list.index(x[0]))
        main_sample_matrix = sorted(main_sample_matrix_tmp, key=lambda x: group_name.index(x[1]))
        
        jump_count, reordered_main_sample_list = reorder_main_sam_semi(data_matrix, main_sample_matrix)
    #
    elif options.run_mode == "supervised":
        IN = open(options.main_sample, 'r')
        main_sample_list = []
        for line in IN :
            tokens = line.strip().split("\t") 
            main_sample_list.append(sample_header.index(tokens[0]))
        [main_sample_list.append(i) for i in range(3, len(sample_header)) if i not in main_sample_list]
        reordered_main_sample_list = main_sample_list
        jump_count = count_N(data_matrix, reordered_main_sample_list)
        #
    else:
        raise ValueError('Run mode must be one of "un-supervised", "semi-supervised", "supervised"')

    # mask cells with uniq component to positive value
    data_matrix = mark_singleton_block(data_matrix)

    if options.color_mode == "contribution":
        data_matrix, contribution_reordered_list = change_order(data_matrix, reordered_main_sample_list)
        final_list = contribution_reordered_list
        #
    elif options.color_mode == "priority":
        _, contribution_reordered_list = change_order(deepcopy(data_matrix), reordered_main_sample_list)
        final_list = reordered_main_sample_list
        #
    else:
        raise ValueError('Color mode must be one of "contribution" and "priority".')

    # change sample order according to their contribution
    
    # report (.ms, .as, .pri, .fi file)
    report_stat(data_matrix, final_list, contribution_reordered_list, reordered_main_sample_list, sample_header, options.outfile)

    # write .grp file
    grp_file = open(options.outfile+".grp", 'w')
    grp_file.write("\t".join(sample_header)+"\n")
    for line in data_matrix:
        grp_file.write("\t".join([str(i) for i in line])+"\n")

# ===========================================
if __name__ == "__main__":
    main()