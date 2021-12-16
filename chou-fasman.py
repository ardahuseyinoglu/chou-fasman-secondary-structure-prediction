# Arda Huseyinoglu

import sys
from itertools import groupby
#import seaborn as sn
#import pandas as pd
#import matplotlib.pyplot as plt



# handle command line args
measure_flag = 0
if len(sys.argv) == 2:
	input_file_name = sys.argv[1]
if len(sys.argv) == 3:
	input_file_name = sys.argv[1]
	gt_ss_file_name = sys.argv[2]
	measure_flag = 1


with open(input_file_name) as file:
    lines = file.readlines()
    lines = [line.rstrip() for line in lines]

protein_header = lines[0]
seq = lines[1]


pt = {}
pt['A'] = {'Pa':142,   'Pb': 83,    'Pt': 66,   'f0': 0.06,    'f1': 0.076,   'f2': 0.035,   'f3': 0.058}
pt['R'] = {'Pa':98,    'Pb': 93,    'Pt': 95,   'f0': 0.070,   'f1': 0.106,   'f2': 0.099,   'f3': 0.085}
pt['N'] = {'Pa':101,   'Pb': 54,    'Pt': 146,  'f0': 0.147,   'f1': 0.110,   'f2': 0.179,   'f3': 0.081}
pt['D'] = {'Pa':67,    'Pb': 89,    'Pt': 156,  'f0': 0.161,   'f1': 0.083,   'f2': 0.191,   'f3': 0.091}
pt['C'] = {'Pa':70,    'Pb': 119,   'Pt': 119,  'f0': 0.149,   'f1': 0.050,   'f2': 0.117,   'f3': 0.128}
pt['E'] = {'Pa':151,   'Pb': 37,    'Pt': 74,   'f0': 0.056,   'f1': 0.060,   'f2': 0.077,   'f3': 0.064}
pt['Q'] = {'Pa':111,   'Pb': 110,   'Pt': 98,   'f0': 0.074,   'f1': 0.098,   'f2': 0.037,   'f3': 0.098}
pt['G'] = {'Pa':57,    'Pb': 75,    'Pt': 156,  'f0': 0.102,   'f1': 0.085,   'f2': 0.190,   'f3': 0.152}
pt['H'] = {'Pa':100,   'Pb': 87,    'Pt': 95,   'f0': 0.140,   'f1': 0.047,   'f2': 0.093,   'f3': 0.054}
pt['I'] = {'Pa':108,   'Pb': 160,   'Pt': 47,   'f0': 0.043,   'f1': 0.034,   'f2': 0.013,   'f3': 0.056}
pt['L'] = {'Pa':121,   'Pb': 130,   'Pt': 59,   'f0': 0.061,   'f1': 0.025,   'f2': 0.036,   'f3': 0.070}
pt['K'] = {'Pa':114,   'Pb': 74,    'Pt': 101,  'f0': 0.055,   'f1': 0.115,   'f2': 0.072,   'f3': 0.095}
pt['M'] = {'Pa':145,   'Pb': 105,   'Pt': 60,   'f0': 0.068,   'f1': 0.082,   'f2': 0.014,   'f3': 0.055}
pt['F'] = {'Pa':113,   'Pb': 138,   'Pt': 60,   'f0': 0.059,   'f1': 0.041,   'f2': 0.065,   'f3': 0.065}
pt['P'] = {'Pa':57,    'Pb': 55,    'Pt': 152,  'f0': 0.102,   'f1': 0.301,   'f2': 0.034,   'f3': 0.068}
pt['S'] = {'Pa':77,    'Pb': 75,    'Pt': 143,  'f0': 0.120,   'f1': 0.139,   'f2': 0.125,   'f3': 0.106}
pt['T'] = {'Pa':83,    'Pb': 119,   'Pt': 96,   'f0': 0.086,   'f1': 0.108,   'f2': 0.065,   'f3': 0.079}
pt['W'] = {'Pa':108,   'Pb': 137,   'Pt': 96,   'f0': 0.077,   'f1': 0.013,   'f2': 0.064,   'f3': 0.167}
pt['Y'] = {'Pa':69,    'Pb': 147,   'Pt': 114,  'f0': 0.082,   'f1': 0.065,   'f2': 0.114,   'f3': 0.125}
pt['V'] = {'Pa':106,   'Pb': 170,   'Pt': 50,   'f0': 0.062,   'f1': 0.048,   'f2': 0.028,   'f3': 0.053}



def find_regions_of(element):
    
    arr = []

    for i in range(0, len(seq) - 6 + 1):

        counter = 0
        start = i
        end = i+5
        for j in range(i, i+6):
            if pt[seq[j]][element] > 100:
                counter += 1

        if counter >= 4:

            right_pos = i+6
            left_pos = i-1

            while (right_pos != len(seq)):

                counter_right = 0
                for k in range(right_pos-3, right_pos+1):
                    if pt[seq[k]][element] < 100:
                        counter_right += 1

                if counter_right == 4:
                    end = right_pos - 4
                    break
                else:
                    right_pos += 1


            while (left_pos != -1):

                counter_left = 0
                for k in range(left_pos, left_pos+4):
                    if pt[seq[k]][element] < 100:
                        counter_left += 1

                if counter_left == 4:
                    start = left_pos + 4
                    break
                else:
                    left_pos -= 1

            arr.append([start, end])
    
    return arr

	
def converter(arr):
    for i in range(len(seq)):
        if arr[i] == 'a':
            arr[i] = 'H'
        elif arr[i] == 'b':
            arr[i] = 'E'
        elif arr[i]  == 't':
            arr[i]  = 'T'
    return arr

	
	
# alpha and beta assignments with extensions
alpha_regions = find_regions_of('Pa')
beta_regions = find_regions_of('Pb')
alpha = ['_'] * len(seq)
beta = ['_'] * len(seq)

for a_region in alpha_regions:
    for pos in range(a_region[0], a_region[1] + 1):
        alpha[pos] = 'a'
        
for a_region in beta_regions:
    for pos in range(a_region[0], a_region[1] + 1):
        beta[pos] = 'b'

alpha_H = alpha.copy()
beta_E = beta.copy()



# find overlapping regions
overlapping_regions = []
flag = 0
start = 0

for pos in range(len(seq)):
    is_overlap = (alpha[pos] == 'a') and (beta[pos] == 'b')
    if is_overlap and (flag == 0):
        start = pos
        flag = 1
        region_size = 1
    elif is_overlap and (flag == 1):
        region_size += 1
    elif (not is_overlap) and (flag == 1):
        flag = 0
        end = start + region_size - 1
        overlapping_regions.append([start,end])
    elif (not is_overlap) and (flag == 0):
        pass


		
# for overlapping regions, remove assignments from alpha and beta according to sum of alpha and sum of beta
for a_region in overlapping_regions:
    count_a = 0
    count_b = 0
    for pos in range(a_region[0], a_region[1] + 1):
        count_a += pt[seq[pos]]['Pa']
        count_b += pt[seq[pos]]['Pb']
    if count_a > count_b:
        for pos in range(a_region[0], a_region[1] + 1):
            beta[pos] = '_'
    elif count_b > count_a:
        for pos in range(a_region[0], a_region[1] + 1):
            alpha[pos] = '_'



# finally, merge alpha and beta predictions into a single array			
ab_prediction = []

for pos in range(len(seq)):
    ab_prediction.append(alpha[pos])

for pos in range(len(seq)):
    if beta[pos] == 'b':
        ab_prediction[pos] = 'b'



# check if size of alpha or beta regions is smaller then 5. If so, remove assignments for this region
grouped_ab = [[k, sum(1 for i in g)] for k,g in groupby(ab_prediction)]
for i in range(1, len(grouped_ab)):
    grouped_ab[i][1] = grouped_ab[i][1] + grouped_ab[i-1][1]

for i in range(len(grouped_ab)):
    number_of_occ = -1
    if i == 0:
        number_of_occ = grouped_ab[i][1]
    else:
        number_of_occ = grouped_ab[i][1] - grouped_ab[i-1][1]
        
    if (grouped_ab[i][0] != '_') and (number_of_occ < 5):
        for k in range(grouped_ab[i-1][1], grouped_ab[i][1]):
            ab_prediction[k] = '_'


			
# turn prediction		
arr = []
turns_T = ['_'] * len(seq) 
for i in range(0, len(seq) - 4 + 1):
    
    total_a = 0
    total_b = 0
    total_t = 0
    bend = 1
    
    for j in range(i, i+4):
        total_a += pt[seq[j]]['Pa'] 
        total_b += pt[seq[j]]['Pb'] 
        total_t += pt[seq[j]]['Pt']
        col = 'f' + str(j-i)
        bend *= pt[seq[j]][col]
    
    cond1 = (bend > 0.000075)
    cond2 = ((total_t / 4) > 100)
    cond3 = (total_t > total_a)
    cond4 = (total_t > total_b)
    
    if cond1 and cond2 and cond3 and cond4:
        for k in range(i,i+4):
            ab_prediction[k] = 't'
            turns_T[k] = 'T'
			



"""
# if removal needed after turn prediction:

grouped_ab = [[k, sum(1 for i in g)] for k,g in groupby(ab_prediction)]
for i in range(1, len(grouped_ab)):
    grouped_ab[i][1] = grouped_ab[i][1] + grouped_ab[i-1][1]

for i in range(len(grouped_ab)):
    number_of_occ = -1
    if i == 0:
        number_of_occ = grouped_ab[i][1]
    else:
        number_of_occ = grouped_ab[i][1] - grouped_ab[i-1][1]
        
    if ((grouped_ab[i][0] == 'a') or (grouped_ab[i][0] == 'b')) and (number_of_occ < 5):
        for k in range(grouped_ab[i-1][1], grouped_ab[i][1]):
            ab_prediction[k] = '_'
"""




with open("output.txt", "w") as text_file:
	text_seq = ''.join(seq)
	text_alpha = ''.join(converter(alpha_H))
	text_beta = ''.join(converter(beta_E))
	text_turn = ''.join(turns_T)
	text_pred = ''.join(converter(ab_prediction))
	
	text_file.write(protein_header)
	text_file.write('\n')
	text_file.write(f'Protein Sequence\t\t\t\t:{text_seq}')
	text_file.write('\n')
	text_file.write(f'Alpha-Helix Assignment\t\t\t\t:{text_alpha}')
	text_file.write('\n')
	text_file.write(f'Beta-Strand Assignment\t\t\t\t:{text_beta}')
	text_file.write('\n')
	text_file.write(f'Turn Assignment\t\t\t\t\t:{text_turn}')
	text_file.write('\n')
	text_file.write(f'Secondary Structure Prediction\t\t\t:{text_pred}')
	

print("\n*output.txt created\n")

if measure_flag == 1:
	# read ground-truth data from the file
	with open(gt_ss_file_name) as file:
		lines = file.readlines()
		lines = [line.rstrip() for line in lines]

	gt_regions = []
	for line in lines:
		splitted_line = line.split()
		
		element = splitted_line[0][0]
		if element == 'S':
			element = 'E'
		
		gt_regions.append([element, [int(splitted_line[1]) - 1, int(splitted_line[2]) - 1]])

	gt_ss_sequence = ['_'] * len(seq)
	for a_region in gt_regions:
		for pos in range(a_region[1][0], a_region[1][1] + 1):
			gt_ss_sequence[pos] = a_region[0]


	for i in range(len(seq)):
		if ab_prediction[i] == 'a':
			ab_prediction[i] = 'H'
		elif ab_prediction[i] == 'b':
			ab_prediction[i] = 'E'
		elif ab_prediction[i]  == 't':
			ab_prediction[i]  = 'T'

	pred_ss_sequence = ab_prediction.copy()



	# performance calculations
	def conf_matrix_row(gt_element):
		count_pred_as_H = 0
		count_pred_as_E = 0
		count_pred_as_T = 0

		for i in range(len(seq)):
			if gt_ss_sequence[i] == gt_element:
				if pred_ss_sequence[i] == 'H':
					count_pred_as_H += 1
				if pred_ss_sequence[i] == 'E':
					count_pred_as_E += 1
				if pred_ss_sequence[i] == 'T':
					count_pred_as_T += 1
					
		return count_pred_as_H, count_pred_as_E, count_pred_as_T


	hh, he, ht = conf_matrix_row('H')
	eh, ee, et = conf_matrix_row('E')
	th, te, tt = conf_matrix_row('T')

	tp_h = hh
	tp_e = ee
	tp_t = tt

	tn_h = ee+et+te+tt
	tn_e = hh+ht+th+tt
	tn_t = hh+he+eh+ee

	fp_h = eh+th
	fp_e = he+te
	fp_t = ht+et

	fn_h = he+ht
	fn_e = eh+et
	fn_t = th+te

	prec_h = tp_h/(tp_h+fp_h)
	prec_e = tp_e/(tp_e+fp_e)
	prec_t = tp_t/(tp_t+fp_t)

	recall_h = tp_h/(tp_h+fn_h)
	recall_e = tp_e/(tp_e+fn_e)
	recall_t = tp_t/(tp_t+fn_t)

	f1_h = (2*prec_h*recall_h)/(prec_h+recall_h)
	f1_e = (2*prec_e*recall_e)/(prec_e+recall_e)
	f1_t = (2*prec_t*recall_t)/(prec_t+recall_t)

	acc_h = (tp_h+tn_h)/(tp_h+tn_h+fp_h+fn_h)
	acc_e = (tp_e+tn_e)/(tp_e+tn_e+fp_e+fn_e)
	acc_t = (tp_t+tn_t)/(tp_t+tn_t+fp_t+fn_t)


	print('\t\t\t    Predicted')
	print('\t\t\tH\tE\tT')
	print(f'\t\tH\t{hh}\t{he}\t{ht}')
	print(f'Ground Truth\tE\t{eh}\t{ee}\t{et}')
	print(f'\t\tT\t{th}\t{te}\t{tt}')

	print("\n")
	print("\tPrec\t\tRecall\t\tF1\t\tAcc")
	print(f"H\t{round(prec_h,4)}\t\t{round(recall_h,4)}\t\t{round(f1_h,4)}\t\t{round(acc_h,4)}")
	print(f"E\t{round(prec_e,4)}\t\t{round(recall_e,4)}\t\t{round(f1_e,4)}\t\t{round(acc_e,4)}")
	print(f"T\t{round(prec_t,4)}\t\t{round(recall_t,4)}\t\t{round(f1_t,4)}\t\t{round(acc_t,4)}")

	print(f'\nOverall Accuracy: {(hh+ee+tt)/(hh+he+ht+eh+ee+et+th+te+tt)}')


"""
# confusion matrix
df_cm = pd.DataFrame([[hh,he,ht], [eh, ee, et], [th, te, tt]], index = ['H', 'E', 'T'], columns = ['H', 'E', 'T'])
ax = sn.heatmap(df_cm, annot=True)
plt.title("Confusion Matrix\n", fontsize =15)
plt.xlabel('Predicted', fontsize = 12)
plt.ylabel('Ground Truth', fontsize = 12)
plt.show()
"""





