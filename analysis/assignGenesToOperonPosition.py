import mysql.connector
from Bio import SeqIO
import pickle

refFasta = SeqIO.read("../sequence.fasta", 'fasta')

cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='t')
cursor = cnx.cursor(buffered=True)
cursor.execute('SELECT operon_name, firstgeneposleft,lastgeneposright,operon_strand from OPERON')
ret = []
genePos={}  # Store found values to here
for (name, left, right,strand) in cursor:
    geneCursor=cnx.cursor(buffered=True)
    if strand=='forward':
        orderBy='asc'
    else:
        orderBy='desc'
    query="SELECT gene_name from GENE where gene_posleft>=%s and gene_posright <=%s and gene_strand='%s' ORDER BY gene_posleft %s"%(left,right,strand,orderBy)
    geneCursor.execute(query)
    curPos=1
    #print query
    for geneName in geneCursor:
        geneName=str(geneName[0])
        if geneName in genePos and genePos[geneName]!=0:
            if genePos[geneName] != curPos:
                print geneName,' has conflict! setting to 0'
                genePos[geneName]=0
        else:
            genePos[geneName]=curPos
        curPos+=1

System.exit(1)
for (k,v) in genePos.iteritems():
    #query="UPDATE GENE SET operon_pos=%s where gene_name='%s'"%(v,str(k))
    #print query
    cursor.execute("UPDATE GENE SET operon_pos=%s where gene_name=%s",(v,k))

'''READ ME BELOW'''
# Need to manually set remaining nulls to 1