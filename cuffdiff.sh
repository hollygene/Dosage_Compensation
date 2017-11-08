#run cuffdiff to find differential expression between each line and its ancestor
cuffdiff -p $THREADS -o ${SampleName}cuffdiff --labels ${SampleName},anc MG1655.ref.gtf \
SRR5344681_cuffquant/abundances.cxb,SRR5344682_cuffquant/abundances.cxb \
SRR5344683_cuffquant/abundances.cxb,SRR5344684_cuffquant/abundances.cxb
