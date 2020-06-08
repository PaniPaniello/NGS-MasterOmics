#1. Perform differential expression analysis between brain and liver using the EdgeR. Present results using a heatmap with hierarchical clustering in rows and columns and colored classification of differentially expressed genes (DEGs), i.e. overexpressed in brain versus liver and the other way around.

docker run -i -t -v /home/therea/tutorial:/tutorial -w /tutorial ceciliaklein/teaching:uvic

export PATH=$PATH:/tutorial/teaching-utils/;

#1. Differential expression analysis edgeR

cd analysis

edgeR.analysis.R --input_matrix ../quantifications/encode.mouse.gene.expected_count.idr_NA.tsv \
                 --metadata /tutorial/data/gene.quantifications.index.tsv \
                 --fields tissue \
                 --coefficient 3 \
                 --output brain_X_liver

awk '$NF<0.01 && $2<-10{print $1"\tover_brain_X_liver"}' edgeR.cpm1.n2.brain_X_liver.tsv > edgeR.0.01.over_brain_X_liver.txt
awk '$NF<0.01 && $2>10 {print $1"\tover_liver_X_brain"}' edgeR.cpm1.n2.brain_X_liver.tsv > edgeR.0.01.over_liver_X_brain.txt

wc -l edgeR.0.01.over*.txt

#Heatmap

awk '$3=="gene"{ match($0, /gene_id "([^"]+).+gene_type "([^"]+)/, var); print var[1],var[2] }' OFS="\t" /tutorial/refs/gencode.vM4.gtf \
| join.py --file1 stdin \
          --file2 <(cat edgeR.0.01.over_brain_X_liver.txt edgeR.0.01.over_liver_X_brain.txt) \
| sed '1igene\tedgeR\tgene_type' > gene.edgeR.tsv

cut -f1 gene.edgeR.tsv \
| tail -n+2 \
| selectMatrixRows.sh - ../quantifications/encode.mouse.gene.TPM.idr_NA.tsv \
| ggheatmap.R --width 5 \
              --hclust complete \
              --height 8 \
              --col_metadata /tutorial/data/gene.quantifications.index.tsv \
              --colSide_by tissue \
              --col_labels labExpId \
              --row_metadata gene.edgeR.tsv \
              --merge_row_mdata_on gene \
              --rowSide_by edgeR,gene_type \
              --row_labels none \
              --log \
              --pseudocount 0.1 \
              --col_dendro \
              --row_dendro \
              --matrix_palette /tutorial/palettes/palDiverging.txt \
              --colSide_palette /tutorial/palettes/palTissue.txt \
              --output heatmap.brain_X_liver.pdf

#2. Perform gene ontology enrichment analysis of the two sets of DEGs using the command line wrapper of GOstats R package for biological processes. Plot results using any graphical representation and discuss results.

awk '{split($10,a,/\"|\./); print a[2]}' /tutorial/refs/gencode.vM4.gtf | sort -u > universe.txt

#brainXliver

# BP
awk '{split($1,a,"."); print a[1]}' edgeR.0.01.over_brain_X_liver.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ BP \
                  --output edgeR.over_brain_X_liver \
                  --species mouse

# MF
awk '{split($1,a,"."); print a[1]}' edgeR.0.01.over_brain_X_liver.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ MF \
                  --output edgeR.over_brain_X_liver \
                  --species mouse

# CC
awk '{split($1,a,"."); print a[1]}' edgeR.0.01.over_brain_X_liver.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ CC \
                  --output edgeR.over_brain_X_liver \
                  --species mouse

#liverXbrain

# BP
awk '{split($1,a,"."); print a[1]}' edgeR.0.01.over_liver_X_brain.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ BP \
                  --output edgeR.over_liver_X_brain \
                  --species mouse

# MF
awk '{split($1,a,"."); print a[1]}' edgeR.0.01.over_liver_X_brain.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ MF \
                  --output edgeR.over_liver_X_brain \
                  --species mouse

# CC
awk '{split($1,a,"."); print a[1]}' edgeR.0.01.over_liver_X_brain.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ CC \
                  --output edgeR.over_liver_X_brain \
                  --species mouse

#to check the tsv files:

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_brain_X_liver.BP.tsv

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_brain_X_liver.MF.tsv

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_brain_X_liver.CC.tsv

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_liver_X_brain.BP.tsv

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_liver_X_brain.MF.tsv

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_liver_X_brain.CC.tsv


#3. Analyze differential splicing using SUPPA between brain and liver for skipping exon, intron retention, mutually exclusive exon and alternative first exon. Plot top results using heatmaps. Different thresholds may be chosen for each event type.

cd ../splicing

#Prepare input files

# List of transcript ids
awk '$3=="transcript" && $0~/gene_type "protein_coding"/{ match($0, /transcript_id "([^"]+)/, id); print id[1] }' /tutorial/refs/gencode.vM4.gtf |sort -u > protein_coding_transcript_IDs.txt

# Genome annotation restricted to exon features and filtered by transcript type
cat /tutorial/refs/gencode.vM4.gtf |awk '$3=="exon"' |grep -Ff protein_coding_transcript_IDs.txt > exon-annot.gtf

# Filter transcript TPM matrix
selectMatrixRows.sh protein_coding_transcript_IDs.txt /tutorial/quantifications/encode.mouse.transcript.TPM.idr_NA.tsv > pc-tx.tsv

# Individual transcript expression matrices
for tissue in Brain Liver; do
    selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 pc-tx.tsv > expr.${tissue}.tsv
done
suppa.py generateEvents -i exon-annot.gtf -e SE MX RI FL -o localEvents -f ioe

#Generate local AS events

SUPPA info in github:  https://github.com/comprna/SUPPA/blob/master/README.md#overview 

suppa.py generateEvents -i exon-annot.gtf -e SE SS MX RI FL -o localEvents -f ioe

wc -l localEvents*ioe

#Compute PSI for all events and create it's respective tissue files and generate the heatmaps:

#SE

event=SE; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}

event=SE; for tissue in Brain Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

event=SE; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}

event=SE; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.5 || $2<-0.5) && $3<0.05{print}' DS.${event}.dpsi|cut -f1 > top-examples-SE.txt
selectMatrixRows.sh top-examples-SE.txt DS.SE.psivec > matrix.top-examples-SE.tsv

ggheatmap.R -i matrix.top-examples-SE.tsv -o heatmap_top-examples-SE.pdf --matrix_palette /tutorial/palettes/palSequential.txt --row_dendro  --matrix_fill_limits "0,1" -B 8

#MX

event=MX; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}

event=MX; for tissue in Brain Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

event=MX; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}

event=MX; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.3 || $2<-0.3) && $3<0.05{print}' DS.${event}.dpsi|cut -f1 > top-examples-MX.txt
selectMatrixRows.sh top-examples-MX.txt DS.MX.psivec > matrix.top-examples-MX.tsv

ggheatmap.R -i matrix.top-examples-MX.tsv -o heatmap_top-examples-MX.pdf --matrix_palette /tutorial/palettes/palSequential.txt --row_dendro  --matrix_fill_limits "0,1" -B 8

#RI
event=RI; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}

event=RI; for tissue in Brain Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

event=RI; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}

event=RI; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.3 || $2<-0.3) && $3<0.05{print}' DS.${event}.dpsi|cut -f1 > top-examples-RI.txt
selectMatrixRows.sh top-examples-RI.txt DS.RI.psivec > matrix.top-examples-RI.tsv

ggheatmap.R -i matrix.top-examples-RI.tsv -o heatmap_top-examples-RI.pdf --matrix_palette /tutorial/palettes/palSequential.txt --row_dendro  --matrix_fill_limits "0,1" -B 8

#AF
event=AF; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}

event=AF; for tissue in Brain Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

event=AF; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}

event=AF; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.5 || $2<-0.5) && $3<0.05{print}' DS.${event}.dpsi|cut -f1 > top-examples-AF.txt
selectMatrixRows.sh top-examples-AF.txt DS.AF.psivec > matrix.top-examples-AF.tsv

ggheatmap.R -i matrix.top-examples-AF.tsv -o heatmap_top-examples-AF.pdf --matrix_palette /tutorial/palettes/palSequential.txt --row_dendro --matrix_fill_limits "0,1" -B 8

#4. Find H3K4me3 peaks shared by brain and liver and the ones exclusively found in each tissue using the narrow peaks found in /tutorial/results using bedtools intersect. Show results using a bar plot colored by the color code used during the hands-on. Palette isavailables at /tutorial/palettes/palTissue.txt. Any color may be chosen for shared peaks.

mkdir ../chip-analysis

cd chip-analysis

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoLiver.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Brain_coordinates","Brain_peak","Liver_coordinates","Liver_peak","intersection"}$NF!=0{print $1":"$2"-"$3,$4,$11":"$12"-"$13,$14,$NF}' > common-peaks-Brain-liver.tsv

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoHeart.narrowPeak -wao|awk 'BEGIN{OFS="\t"; print "Brain_coordinates","Brain_peak","Liver_coordinates","Liver_peak","intersection"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14}' > common-peaks-Brain-liver.bed

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoLiver.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Brain_coordinates","Brain_peak"}$NF==0{print $1":"$2"-"$3,$4}' > Brain-specific-peaks.tsv

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoHeart.narrowPeak -wao|awk 'BEGIN{OFS="\t"; print "Brain_coordinates","Brain_peak"}$NF!=0{print $1,$2,$3,$4}' > Brain-specific-peaks.bed

bedtools intersect -a /tutorial/results/CHIPembryoLiver.narrowPeak -b /tutorial/results/CHIPembryoBrain.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Liver_coordinates","Liver_peak"}$NF==0{print $1":"$2"-"$3,$4}' > Liver-specific-peaks.tsv

bedtools intersect -a /tutorial/results/CHIPembryoLiver.narrowPeak -b /tutorial/results/CHIPembryoBrain.narrowPeak -wao|awk 'BEGIN{OFS="\t"; print "Liver_coordinates","Liver_peak"}$NF!=0{print $1,$2,$3,$4}' > Liver-specific-peaks.bed

#barplot:

ls *tsv|while read f; do echo -e $f"\t"$(grep -v nan $f |wc -l);done |ggbarplot.R -i stdin -o number_of_peaks.pdf -f 1 --palette_fill /tutorial/palettes/palTissue.txt --title "H3k4me3 peaks” --y_title "Number of peaks" --x_title "Tissues"

#5. Create a BED file of 200bp up/downstream TSS of genes and overlap DEGs (step 1)with the 3 sets of H3K4me3 peaks classified in the previous step (4). Show three examples in the UCSC genome browser, including RNA-seq, ChIP-seq and ATAC-seqtracks. Ideally, one example of each peak set (i.e. shared peak, peak exclusively calledin brain and peak exclusively called in liver). Discuss the integration of the three datasetsin the TSS of the selected cases.

mkdir /tutorial/tss

# Removing mithocondrial chromosomes:

awk 'BEGIN{FS=OFS="\t"}$1!="chrM" && $3=="transcript" && $0~/gene_type "protein_coding"/ && $7=="+"{ match($0, /transcript_id "([^"]+)/, id); printf ("%s\t20%s\t20%s\t%s\t%s\t\n",$1,$4-200,$4+200,$7,id[1]) }$3=="transcript" && $0~/gene_type "protein_coding"/ && $7=="-"{ match($0, /transcript_id "([^"]+)/, id); printf ("%s\t20%s\t20%s\t%s\t%s\t\n",$1,$5-200,$5+200,$7,id[1]) }' /tutorial/refs/gencode.vM4.gtf |head -n -1 > /tutorial/tss/protein-coding-transcript-200up_downTSS.bed

awk 'BEGIN{FS=OFS="\t"}$3=="gene" && $0~/gene_type "protein_coding"/ && $7=="+"{ match($0, /transcript_id "([^"]+)/, id); print $1,$5-200,$5+200,id[1],$7 }$3=="gene" && $0~/gene_type "protein_coding"/ && $7=="-"{ match($0, /transcript_id "([^"]+)/, id); print $1,$5-200,$5+200,id[1],$7 }' /tutorial/refs/gencode.vM4.gtf |sort -u > /tutorial/tss/protein-coding-genes-200up_downTSS.bed

cut -f1 /tutorial/analysis/edgeR.0.01.over_liver_X_brain.txt |sort -u| grep -Ff - /tutorial/tss/protein-coding-genes-200up_downTSS.bed > /tutorial/tss/degs-over-liver-TSS.bed

cut -f1 /tutorial/analysis/edgeR.0.01.over_brain_X_liver.txt |sort -u| grep -Ff - /tutorial/tss/protein-coding-genes-200up_downTSS.bed > /tutorial/tss/degs-over-brain-TSS.bed

#Intersect

cd tss

cat /tutorial/Brain-specific-peaks.bed | tail -n+2 | sort -u > brain-header-peaks.bed

cat /tutorial/Liver-specific-peaks.bed | tail -n+2 | sort -u > liver-header-peaks.bed

cat /tutorial/common-peaks-Brain-liver.bed | tail -n+2 | sort -u > common-header-peaks.bed

bedtools intersect -a degs-over-brain-TSS.bed -b brain-header-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > brain-peaks-over-brain.bed

bedtools intersect -a degs-over-liver-TSS.bed -b liver-header-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > liver-peaks-over-liver.bed

bedtools intersect -a degs-over-brain-TSS.bed -b liver-header-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > liver-peaks-over-brain.bed

bedtools intersect -a degs-over-liver-TSS.bed -b brain-header-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > brain-peaks-over-liver.bed

bedtools intersect -a degs-over-liver-TSS.bed -b common-header-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > common-peaks-over-liver.bed

bedtools intersect -a degs-over-brain-TSS.bed -b common-header-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > common-peaks-over-brain.bed
