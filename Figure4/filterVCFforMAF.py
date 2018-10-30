import re
from operator import itemgetter
import sys

def parse_vcf_byline(vcf_line):
    vcf_line = vcf_line.strip()
    array = line.split('\t')
    if len(array) > 9 and line[0] != '#':
        chrom, position, ref, alts = array[0], array[1], array[3], array[4]
        info = array[7]
        VT = re.search('[V][T][=][A-Z,\,]+', info).group()[3:]
        AF = re.search('[\;][A][F][=][0-9,\.]+', info).group()[4:]
        return [chrom,position,VT,AF,ref,alts]

def filters(AF, Freq):
    if type(AF) is float:
        if Freq == 'High' and AF >= 0.05:
            return True
        if Freq == 'Low' and 0.01 <= AF < 0.05:
            return True
    else:
        return 'Error!'

def apply_filters(parsed_vcf_byline, Freq):
    #print(parsed_vcf_byline)
    output = []
    chrom, position, VT, AF, ref, alts = parsed_vcf_byline

    if ',' not in alts:
        alt = alts
        flank=5
        start = int(position) - flank - 1
        end = int(position) + len(ref) + flank -1
        alt_AF = float(AF)
        ref_AF = 1 - alt_AF

        if len(ref) != len(alt):
            VT = 'INDEL'
        else:
            VT = 'SNP'


        if alt_AF > ref_AF:
            if filters(ref_AF, Freq) is not None:
                cat_freqs = str(ref_AF) + '|' + str(alt_AF)
                if VT == 'SNP':
                    output.append([chrom,position,position,cat_freqs,ref,alt,VT])
                else:
                    output.append([chrom,start,end,cat_freqs,ref,alt,VT])
        else:
            if filters(alt_AF, Freq) is not None:
                cat_freqs = str(ref_AF) + '|' + str(alt_AF)
                if VT == 'SNP':
                    output.append([chrom,position,position,cat_freqs,ref,alt,VT])
                else:
                    output.append([chrom,start,end,cat_freqs,ref,alt,VT])

        return output

infile = open(sys.argv[1], 'r')

lowfreqsnp = []
lowfreqindel = []

highfreqsnp = []
highfreqindel = []

for line in infile:
    parsed_vcf_byline = parse_vcf_byline(line)
    if parsed_vcf_byline:
        filtered = apply_filters(parsed_vcf_byline,'Low')
        if filtered:
            for call in filtered:
                #print(call)
                if call[6] == 'SNP':
                    lowfreqsnp.append(call)
                if call[6] == 'INDEL':
                    lowfreqindel.append(call)

        filtered = apply_filters(parsed_vcf_byline,'High')
        if filtered:
            for call in filtered:
                #print(call)
                if call[6] == 'SNP':
                    highfreqsnp.append(call)
                if call[6] == 'INDEL':
                    highfreqindel.append(call)

out1 = open('lowfreqsnp'+sys.argv[2], 'w')
out2 = open('lowfreqindel'+sys.argv[2], 'w')

for snp in lowfreqsnp:
    chrom, start, end, freq, ref, alt, VT = snp
    out1.write('chr'+str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(freq)+'\t'+str(ref)+'\t'+str(alt)+'\n')

for indel in lowfreqindel:
    chrom, start, end, freq, ref, alt, VT = indel
    out2.write('chr'+str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(freq)+'\t'+str(ref)+'\t'+str(alt)+'\n')


out3 = open('highfreqsnp'+sys.argv[2], 'w')
out4 = open('highfreqindel'+sys.argv[2], 'w')

for snp in highfreqsnp:
    chrom, start, end, freq, ref, alt, VT = snp
    out3.write('chr'+str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(freq)+'\t'+str(ref)+'\t'+str(alt)+'\n')

for indel in highfreqindel:
    chrom, start, end, freq, ref, alt, VT = indel
    out4.write('chr'+str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(freq)+'\t'+str(ref)+'\t'+str(alt)+'\n')