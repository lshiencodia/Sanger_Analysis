#!/opt/anaconda/bin//python3

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
import os
import glob
from Bio.SeqRecord import SeqRecord
import sys

# In[34]:

# alignment = AlignIO.read("APR2-8-psiR3_result.aln", "clustal")
# help(AlignIO)
#read the reference sequences
refseq = SeqIO.read(sys.argv[1], "fasta")
#print(refseq.name)
refstart=refseq.seq[0:10]
refend=refseq.seq[-10:]
#print(refstart)
#print(refend)

# In[2]:

list_dna=[]
list_translate=[]
for FileList in glob.glob("*psiR3_result.aln"):
#for FileList in glob.glob("*3a-psiR3_result.aln"):
#for FileList in glob.glob("ClpS2_V1_48-9_psiR3_result.aln"):
    alignment = AlignIO.read(FileList, "clustal")
    startt=0
    endtt=0
   
    #need to create a new instance, somehow the sequences is not destroyed
    colonyseq = Seq("GATCGATC")
    targetseq = Seq("AGTCAGTC")

    for r in alignment:
        if r.id==refseq.name:
            startt=r.seq.find(refstart)
            endt=r.seq.find(refend)+len(refend)
            targetseq=r.seq[startt:endt]
            #print("startt: "+str(startt)+" end:"+str(endt))
        
    for r in alignment:
        if r.id!=refseq.name and startt!=0 and endt!=0:
            colonyid=r.id
            colonyseq=r.seq[startt:endt]

            
    #print(str(len(colonyseq))+" from "+FileList)

    #only deal with same length sequencing reads (no indel)
    if (colonyseq!="GATCGATC"):
      if (len(str(targetseq).replace("-",""))==len(str(colonyseq).replace("-",""))):
          N_number=colonyseq.count('N')
          seqi=SeqRecord(seq=colonyseq,id=colonyid,description="; "+str(N_number)+" Ns from ab1")
          if N_number<10:
            colonyseq_list=list(colonyseq)
            colonyseq_list_copy=[]
            for i in range(len(colonyseq)):
              if(colonyseq_list[i]!="N"):
                  colonyseq_list_copy.append(colonyseq[i])
              else:
                  colonyseq_list_copy.append(targetseq[i])
            colonyseq_fix=Seq("".join(colonyseq_list_copy))
      
#            print(colonyid)
            if "-" not in colonyseq_fix:
              seq1=SeqRecord(seq=colonyseq_fix.translate(),id="trans_"+colonyid,description="; "+str(N_number)+" Ns from ab1")
              #if seq1 not in list_translate:
              list_translate.append(seq1)
              print(FileList+" good ")
      
          else:
            print(FileList+" failed_#N10 "+str(N_number))
      else:
            seqi=SeqRecord(seq=colonyseq,id=colonyid,description="; "+str(len(str(colonyseq).replace("-","")))+" indels")
            print(FileList+" indels "+str(len(str(colonyseq).replace("-",""))-len(str(targetseq).replace("-",""))))
            #print(FileList+" indels "+str(len(str(colonyseq).replace("-","")))+" target "+ str(len(str(targetseq).replace("-",""))))

      list_dna.append(seqi)

# In[3]:
#print(list_translate)

SeqIO.write(list_dna, "all_dna.fa", "fasta")
SeqIO.write(list_translate, "translate_all.fa", "fasta")
# SeqIO.write(list_translate, "translate_"+FileList.replace(".aln",".fa"), "fasta")
#print("target")
#print(targetseq)
#print(str(targetseq).replace("-",""))
seqtarget=SeqRecord(seq=Seq(str(targetseq).replace("-","")).translate(),id="trans_target")
SeqIO.write(seqtarget, "target.fa", "fasta")


# In[6]:

alignment = AlignIO.read("translate_all.fa", "fasta")

#find everything that has mutations


#select subset of the sequences
#designed_pos=[28,30,31,32,33,34,35,36,39,58,61,62,95]
designed_pos=[28,39,58,61,95,31,34,30,62,32,33,35,36]
designed_ID=[i-2 for i in designed_pos]


# In[10]:

edited=alignment[:,designed_ID[0]:(designed_ID[0]+1)]
for i in range(1,len(designed_ID)):
    edited+=alignment[:,designed_ID[i]:(designed_ID[i]+1)]
# print(edited)
AlignIO.write(edited,"designed_pos.fa","fasta")

resfile=open("designed_pos.position","w")
wt_seq_pos=[]
for i in designed_ID:
	wt_seq_pos.append(str(seqtarget[i])+str(i+2))
resfile.write(",".join(wt_seq_pos))
resfile.close()


# In[ ]:



