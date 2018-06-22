#convert to csv files
cat target.fa translate_all.fa > translate_all_target.fa
seqkit fx2tab translate_all_target.fa > translate_all_target.csv
seqkit fx2tab designed_pos.fa > designed_pos.csv
awk '{print $1}' designed_pos.csv | sed -e 's/trans_//g' -e 's/-9_psiR3//g' > designed_pos.name
awk '{print $1}' designed_pos.csv | sed 's/trans_//g' | awk -F"-" '{print $1}' > designed_pos.number
awk '{print $3}' designed_pos.csv > designed_pos.qual
awk '{print $NF}' designed_pos.csv > designed_pos.seq
echo "Sequences,Name,Number,#NinSanger" > designed_pos_format.csv
paste -d "," designed_pos.seq designed_pos.name  designed_pos.number designed_pos.qual >>  designed_pos_format.csv
/bin/rm -f designed_pos.seq designed_pos.name designed_pos.number designed_pos.qual

header=`cat designed_pos.position`
#put in a good format for analysis
for x in $(/bin/ls designed_pos.csv | grep -v char)
do
        awk '{print $1}' ${x%%.*}.csv > ${x%%.*}.name
        awk '{print $3}' ${x%%.*}.csv > ${x%%.*}.qual
        awk -F"\t" '{print $2}' ${x%%.*}.csv | sed 's/\(.\)/\1,/g'  | sed 's/.$//' > ${x%%.*}_seq_char.csv
        #awk -F"\t" '{print $2}' ${x%%.*}.csv | awk NF=NF FS= | sed 's/ /,/g' > ${x%%.*}_seq_char.csv
        awk -F"\t" '{print $2}' ${x%%.*}.csv > ${x%%.*}_seq.csv
        echo "seqname,quality(#N),$header" > ${x%%.*}_aln_char.csv
        paste -d "," ${x%%.*}.name ${x%%.*}.qual ${x%%.*}_seq_char.csv >> ${x%%.*}_aln_char.csv
        paste -d "," ${x%%.*}.name ${x%%.*}.qual ${x%%.*}_seq.csv > ${x%%.*}_aln.csv
        rm -f ${x%%.*}.name ${x%%.*}.qual ${x%%.*}_seq.csv ${x%%.*}_seq_char.csv
done
