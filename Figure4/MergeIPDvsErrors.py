# Merge motifs_and_errors and Windows_Collected_F (errors and IPD)
import pandas as pd
import sys
import numpy as np

errors = open(sys.argv[1], 'rt') # Using collapsed file
ipds = open(sys.argv[2], 'rt') # Using mf file
divergence = open(sys.argv[3], 'rt') # Using collapsed
diversity = open(sys.argv[4], 'rt') # Using collapsed


IPD_rows = []
IPD_colnames = ['chrom', 'start', 'end', 'motif', 'length', 'meanIPD']
errors_rows = []
errors_colnames = ['chrom', 'start','end', 'SNP', 'INS', 'DEL']

for line in ipds:
    line = line.strip()
    array = line.split('\t')
    #chrom, start, end, motif, length = array[0], array[1], array[2], array[3], array[4]
    motif = 'GQuadPlus'
    chrom, start, end, length = array[0][3:], array[1], array[2], array[3]
    start = str(int(float(start)))
    end = str(int(float(start)))
#    print(chrom, start, end, length)
    meanIPD = np.mean([float(i) for i in array[4:] if i != 'NA'])
    IPD_rows.append([chrom,start,end,motif,length,meanIPD])

df_IPD = pd.DataFrame(data=IPD_rows, columns=IPD_colnames)

#print('IPD')
#print(df_IPD.head())
#print(df_IPD.shape)

for line in errors:
    line = line.strip()
    array = line.split('\t')
#    array = line.split(' ')
#    print(array)
#    chrom, start, end, errors = array[0][3:], array[1], array[2], array[3]
    chrom, start, end, errors = array[0], array[1], array[2], array[3]
    total_nouse, SNP, INS, DEL = errors.split(' ')
    #chrom, start, end, total_nouse, SNP, INS, DEL = array[0], array[1], array[2], array[3], array[4], array[5], array[6]
    start = str(int(float(start)))
    end = str(int(float(start)))
    errors_rows.append([chrom, start, end, SNP, INS, DEL])

df_ERRORS = pd.DataFrame(data=errors_rows, columns=errors_colnames)

#print('ERRORS')
#print(df_ERRORS.head())

df_IPDvsERRORS = pd.merge(df_IPD, df_ERRORS, on=['chrom','start', 'end'])

#print('IPD|ERRORS')
#print(df_IPDvsERRORS.head())

divergence_rows = []
divergence_rows_colnames = ['chrom', 'start','end', 'divergence_SNP', 'divergence_INS', 'divergence_DEL']
                                                                                                                                            61,1          36%
for line in divergence:
    line = line.strip()
    array = line.split('\t')
    chrom, start, end, total_nouse, SNP, INS, DEL = array[0], array[1], array[2], array[3], array[4], array[5], array[6]
    start = str(int(float(start)))
    end = str(int(float(start)))
    divergence_rows.append([chrom, start, end, SNP, INS, DEL])

df_divergence = pd.DataFrame(data=divergence_rows, columns=divergence_rows_colnames)

df_IPDvsERRORSvsDiv = pd.merge(df_IPDvsERRORS, df_divergence, on=['chrom','start', 'end'])

#print('IPD|ERRORS|DIVERGENCE')
#print(df_IPDvsERRORSvsDiv.head())

diversity_rows = []
diversity_colnames = ['chrom', 'start','end', 'diversity_SNP', 'diversity_INS', 'diversity_DEL']

for line in diversity:
    line = line.strip()
    array = line.split('\t')
    chrom, start, end, total_nouse, SNP, INS, DEL = array[0], array[1], array[2], array[3], array[4], array[5], array[6]
    start = str(int(float(start)))
    end = str(int(float(start)))
    diversity_rows.append([chrom, start, end, SNP, INS, DEL])

df_diversity = pd.DataFrame(data=diversity_rows, columns=diversity_colnames)

#print('DIVERISTY')
#print(df_diversity.head())

df_IPDvsERRORSvsDiv2 = pd.merge(df_IPDvsERRORSvsDiv, df_diversity, on=['chrom','start', 'end'])

#print('IPD|ERRORS|DIVERGENCE|DIVERSITY')
#print(df_IPDvsERRORSvsDiv2.head())

#for index, row in df_IPDvsERRORSvsDiv.iterrows():
#    print(row['chrom']+'\t'+row['start']+'\t'+row['end']+'\t'+row['motif']+'\t'+row['length']
#         +'\t'+str(row['meanIPD'])+'\t'+row['SNP']+'\t'+row['INS']+'\t'+row['DEL']
#         +'\t'+row['div_SNP']+'\t'+row['div_INS']+'\t'+row['div_DEL'])

print('CHROM'+'\t'+'START'+'\t'+'END'+'\t'+'MOTIF'+'\t'+'MEANIPD'+'\t'+'PACBIOSNP'+'\t'+'DIVERGENCESNP'+'\t'+'DIVERSITYSNP')
for index, row in df_IPDvsERRORSvsDiv2.iterrows():
    print(row['chrom']+'\t'+row['start']+'\t'+row['end']+'\t'+row['motif']+'\t'+str(row['meanIPD'])+
          '\t'+row['SNP']+'\t'+row['divergence_SNP']+'\t'+row['diversity_SNP'])

