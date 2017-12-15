import os
path = 'DECODE/'
filelist = []

for filename in os.listdir(path):
    if filename[0:3] == 'Pro' and filename[-3:] == 'vcf':
        filelist.append(path+filename)
        
def extract_calls(DecodeFile):
    variants = []
    for line in DecodeFile:
        if line[0] != "#" and 'PASS' in line:
            array = line.split('\t')
            chrom, pos, id_nouse, ref, alt = array[0:5]
            if len(ref) == len(alt):
                VT = 'SNP'
            if len(ref) < len(alt):
                VT = 'INS'
            if len(ref) > len(alt):
                VT = 'DEL'
            variants.append([chrom,pos,VT,ref,alt])
    return variants

DecodeVariants = []
output = open(path+'DecodeDNM.gff', 'w')

for DecodeFile in filelist:
    DecodeVariants = DecodeVariants + extract_calls(open(DecodeFile,'rt'))
    
print(len(DecodeVariants))
    
for variant in DecodeVariants:
    chrom, pos, VT, ref, alt = variant
    output.write(chrom[3:]+'\t1kG\t'+VT+'\t'+pos+'\t'+pos+'\t.\t.\t.\n')

print('DONE')