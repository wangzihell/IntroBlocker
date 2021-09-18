#!/usr/bin/env python
import sys
import copy

if len(sys.argv) < 2:
	chrom = "chr1A"
else:
	chrom = sys.argv[1]

if len(sys.argv) < 3:
	grp_data_folder = "/data2/rawdata2/tetraintro/210109/data"
else:
	grp_data_folder = sys.argv[2]

if len(sys.argv) < 4:
	dist_data_folder = "/data2/rawdata2/tetraintro/210109/repr_dist_AABB"
else:
	dist_data_folder = sys.argv[3]
	
if len(sys.argv) < 5:
	new_grp_output_folder = "/data3/user3/wangwx/projs/tetraintro/210109/"
else:
	new_grp_output_folder = sys.argv[4]
	
if len(sys.argv) < 6:
	vote_dis_threshold = 593 # 593: 40alpha=beta    828: 2alpha=beta
else:
	vote_dis_threshold = int(sys.argv[5])
	
if len(sys.argv) < 7:
	vote_range = "3"
else:
	vote_range = sys.argv[6]
	
grp_file = grp_data_folder + "/" + chrom + ".grp"
ms_file = grp_data_folder + "/" + chrom + ".ms"
dist_order_file = dist_data_folder + "/name.txt"

f = open(ms_file, 'r')
ms_mat = f.readlines()
f.close()
ms_dic = {}
for i in range(len(ms_mat)):
	ms_line = ms_mat[i].strip().split("\t")
	grp_label = str(abs(int(ms_line[1])))
	if grp_label in ms_dic:
		print("Grp label" + grp_label + "appears twice in " + chrom + ".ms" + " file.")
		exit()
	ms_dic[grp_label] = ms_line[0]

f = open(dist_order_file, 'r')
do_mat = f.readlines()
f.close()
dist_index_dic = {}
for i in range(len(do_mat)):
	do_line = do_mat[i].strip()
	dist_index_dic[do_line] = i

f = open(grp_file, 'r')
sample_list = (f.readline()).strip().split("\t")
grp_mat = f.readlines()
f.close()

new_grp_mat = copy.deepcopy(grp_mat)

for i in range(len(grp_mat)):
	grp_mat[i] = (grp_mat[i]).strip().split("\t")

for i in range(len(grp_mat)):
	grp_line = grp_mat[i]
	new_grp_line = (new_grp_mat[i]).strip().split("\t")
	dist_file = dist_data_folder + "/" + chrom + "." + str(i+1) + ".dist.txt"
	dist_file_io = open(dist_file, "r")
	dist_mat = dist_file_io.readlines()
	dist_file_io.close()
	for j in range(3, len(grp_line)):
		dist_line = (dist_mat[dist_index_dic[sample_list[j]]]).strip().split("\t")
		try:
			curr_grp = abs(int(grp_line[j]))
		except:
			print("Error in reading grp label: line " + str(i+1) + "  col " + str(j))
			exit()
		if curr_grp != 0:
			# if not CNV
			if int(dist_line[dist_index_dic[ms_dic[str(curr_grp)]]]) > vote_dis_threshold:
				# if alpha > beta
				l_vote_bin_index = i - int(vote_range)
				if l_vote_bin_index < 0:
					l_vote_bin_index = 0
				r_vote_bin_index = i + int(vote_range)
				if r_vote_bin_index > len(grp_mat)-1:
					r_vote_bin_index = len(grp_mat)-1
				# count votes
				vote_grp_list = []
				vote_count_list = []
				for k in range(l_vote_bin_index, r_vote_bin_index+1):
					vote_grp_line = grp_mat[k]
					curr_vote_grp = abs(int(vote_grp_line[j]))
					if curr_vote_grp == 0:
						continue
					try:
						vote_grp_index = vote_grp_list.index(curr_vote_grp)
					except:
						vote_grp_list.append(curr_vote_grp)
						vote_count_list.append(0)
						vote_grp_index = vote_grp_list.index(curr_vote_grp)
					vote_count_list[vote_grp_index] += 1
				if len(vote_grp_list) > 0:
					max_vote = max(vote_count_list)
					vote_max_grp_list = []
					for k in range(len(vote_count_list)):
						if vote_count_list[k] == max_vote:
							vote_max_grp_list.append(vote_grp_list[k])
					if len(vote_max_grp_list) > 1:
						distance_list = []
						for k in range(len(vote_max_grp_list)):
							distance_list.append(int(dist_line[dist_index_dic[ms_dic[str(vote_max_grp_list[k])]]]))
						min_dist = min(distance_list)
						dist_min_grp_list = []
						for k in range(len(distance_list)):
							if distance_list[k] == min_dist:
								dist_min_grp_list.append(vote_max_grp_list[k])
						if len(dist_min_grp_list) > 1:
							try:
								_ = dist_min_grp_list.index(curr_grp)
								print("Warning: block on sample " + sample_list[j] + " at bin No." + str(i) + " of " + chrom  + " have " + str(len(dist_min_grp_list)) + "grps that vote same counts and have same distances. Grp didnot changed because its old grp was contained in that " + str(len(dist_min_grp_list)) + " grps.")
							except:
								new_grp_line[j] = str(-dist_min_grp_list[0])
								print("Warning: block on sample " + sample_list[j] + " at bin No." + str(i) + " of " + chrom  + " have " + str(len(dist_min_grp_list)) + "grps that vote same counts and have same distances. Grp was changed by the first grp in the " + str(len(dist_min_grp_list)) + " grps.")
					else:
						if curr_grp == vote_max_grp_list[0]:
							print("Message: block on sample " + sample_list[j] + " at bin No." + str(i) + " of " + chrom  + " have not been changed because its old grp is the winer grp in voting.")
						else:
							new_grp_line[j] = str(-vote_max_grp_list[0])
							print("Message: block on sample " + sample_list[j] + " at bin No." + str(i) + " of " + chrom  + " have been changed successfully by voteing.")
				else:
					print("Warning: block on sample " + sample_list[j] + " at bin No." + str(i) + " of " + chrom  + " was faild in voting for no valid votes there.")
			else:
				# alpha < beta
				pass
		else:
			# CNV
			pass
	new_grp_line = "\t".join(new_grp_line) + "\n"
	new_grp_mat[i] = new_grp_line
	
# finished
f = open(new_grp_output_folder + "/" + chrom + ".smoothed.grp", 'w')
f.write("\t".join(sample_list) + "\n")
f.writelines(new_grp_mat)
f.close()
