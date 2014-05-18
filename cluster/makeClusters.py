from pandas.io.parsers import *
import numpy as np
import mysql.connector
from pandas.io.sql import *

df = read_csv(open('selected.csv'))
df.index = df['seq']
cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')

NUM_SGRNA = len(df.iloc[:,0])
proportions = df.iloc[:, 1:] + 1  # skip first column which contains 20-mer sequence name
proportions = proportions.iloc[:,:].apply(lambda x:x/sum(x),1)  # Normalize to proportions

sgRNAs_with_gene = read_frame('select * from sgRNAs where gene_name!=" "', cnx, index_col='seq')
df = df.merge(sgRNAs_with_gene, 'left', left_index=True, right_index=True)
column_names = {'control':('1_aerobic_t30_1_control', '1_aerobic_t60_1_control', '1_aerobic_t180_1_control'),
              'chlor':('1_aerobic_t30_2_chlor', '1_aerobic_t60_2_chlor', '1_aerobic_t180_2_chlor'),
              'nor':('1_aerobic_t30_6_nor', '1_aerobic_t60_6_nor', '1_aerobic_t180_6_Nor')  }
conditions = column_names.keys()[1:]  # assume first entry at index 0 is the control
scores = np.zeros((NUM_SGRNA, len(conditions)))
is_missing = np.zeros(NUM_SGRNA)
for idx, condition in enumerate(conditions):
    wts = np.empty((NUM_SGRNA, 3))
    for j in xrange(3):
        wts[:,j] = df.loc[:, column_names['control'][j]] + df.loc[:, column_names[condition][j]]

    wt_sums = np.sum(wts, axis=1)
    is_missing = np.logical_or(wt_sums==0, is_missing)  # element-wise logical OR
    for j in xrange(3):
        scores[:,idx] += (wts[:,j]/wt_sums) * (proportions.loc[:, column_names[condition][j]] / proportions.loc[:, column_names['control'][j]])

print "Done computing scores"

def makeOutputFile (outFile, filtered_scores):
    outFile.write('GENE')
    for condition in conditions:
        outFile.write('\t'+condition)
    outFile.write("\n")

    for j in xrange(filtered_scores.shape[0]):
        outFile.write(df.iloc[j,0])
        for score in filtered_scores[j, :]:
            outFile.write('\t')
            outFile.write(str(score))
        outFile.write('\n')

sgRNAs_with_assigned_gene = ~df.loc[:,'gene_name'].notnull().values

makeOutputFile(open('filtered.txt', 'w'), scores[~is_missing, :])
makeOutputFile(open('sgRNAs_for_genes_only.txt', 'w'), scores[(~is_missing) & sgRNAs_with_assigned_gene, :])
print np.sum(is_missing),'sgRNA\'s skipped'