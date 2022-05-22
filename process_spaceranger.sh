spaceranger_path='./spaceranger-1.2.2'
transcriptome_path='./refdata-gex-mm10-2020-A'
fastq_path='./data/filter'


bsub -q TEST-A -n 16 -e run.err -o run.log -x "$spaceranger_path/spaceranger count --id=E135 --transcriptome=$transcriptome_path  --fastqs=$fastq_path --sample=E136 --image=$fastq_path/b1_he.tif --slide=V10T10-322 --area=B1 --localcores=16 --slidefile $spaceranger_path/V10T10-322.gpr --reorient-images"