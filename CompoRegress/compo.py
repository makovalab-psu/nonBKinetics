import sys

bed = open(sys.argv[1], 'r')
seq = open(sys.argv[2], 'r')

i = 0
compos = {}

def occurrences(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count

for line in seq:
        if line[0] == '>':

                seqline = seq.next().strip()
                seqline = seqline.replace('a', 'A')
                seqline = seqline.replace('t', 'T')
                seqline = seqline.replace('g', 'G')
                seqline = seqline.replace('c', 'C')

                A = seqline.count('A')
                T = seqline.count('T')
                G = seqline.count('G')
                C = seqline.count('C')

                total = A + T + G + C

                compo = 'A' + str(A) + 'T' + str(T) + 'G' + str(G) + 'C' + str(C)

                AA = occurrences(seqline, 'AA')
                AT = occurrences(seqline, 'AT')
                AG = occurrences(seqline, 'AG')
                AC = occurrences(seqline, 'AC')

                TA = occurrences(seqline, 'TA')
                TT = occurrences(seqline, 'TT')
                TG = occurrences(seqline, 'TG')
                TC = occurrences(seqline, 'TC')

                GA = occurrences(seqline, 'GA')
                GT = occurrences(seqline, 'GT')
                GG = occurrences(seqline, 'GG')
                GC = occurrences(seqline, 'GC')

                CA = occurrences(seqline, 'CA')
                CT = occurrences(seqline, 'CT')
                CG = occurrences(seqline, 'CG')
                CC = occurrences(seqline, 'CC')


                #compo = 'AA'+str(AA)+'AT'+str(AT)+'AG'+str(AG)+'AC'+str(AC)+'TA'+str(TA)+'TT'+str(TT)+'TG'+str(TG)+'TC'+str(TC)+'GA'+str(GA)+'GT'+str(GT)+'GG'+str(GG)+'GC'+str(GC)+'CA'+str(CA)+'CT'+str(CT)+'CG'+str(CG)+'CC'+str(CC)

                compos[i] = compo
                i += 1

for i, line in enumerate(bed):
        if i in compos:
                print(str(compos[i]) +'\t'+ str(line.strip()))