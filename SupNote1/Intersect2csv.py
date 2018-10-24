#Intersect2csv.py
import pandas as pd

intersectfile = open('GQuadPlus2kflanks_chr21Ends.intersect', 'rt')

colnames = ['chrom', 'motif', 'window_start', 'window_end', 'score', 'motif_start', 'motif_end', 'event_coord', 'ref', 'alt','event_coord_inwindow', 'event_relative_pos']
intersects = []

for line in intersectfile:
    line = line.strip()
    array = line.split('\t')
    # in next line, not using allele frequency since using a merged-uniq gff file
    chrom, motif, window_start, window_end, score, motif_start, motif_end, event_coord, ref, alt = array[0], array[2], int(array[3]), int(array[4]), array[5], int(array[7]), int(array[8]), int(array[12]), array[16], array[17]
    if event_coord < motif_start:
        event_coord_inwindow = str(motif_start - event_coord)
        event_relative_pos = 'L'
    if event_coord > motif_end:
        event_coord_inwindow = str(event_coord - motif_end)
        event_relative_pos = 'R'
    if motif_start <= event_coord <= motif_end:
        event_coord_inwindow = str(event_coord - motif_start)
        event_relative_pos = 'C'
    intersect = [chrom, motif, window_start, window_end, score, motif_start, motif_end, event_coord, ref, alt, event_coord_inwindow, event_relative_pos]
    intersects.append(intersect)

df = pd.DataFrame(data=intersects, columns=colnames)
pd.DataFrame.to_csv(df, 'GQuadPlus2kflanks_chr21Ends.tab', sep = '\t')