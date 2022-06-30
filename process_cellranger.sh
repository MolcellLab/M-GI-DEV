
cellranger_path='opt/cellranger-4.0.0'
transcriptome_path='opt/refdata-gex-mm10-2020-A'
fastq_path='opt/Clean_data'

## sample E95 for example
bsub -q TEST-A -n 16 -e run.err -o run.log "$cellranger_path/cellranger count --id=E95 --transcriptome=$transcriptome_path  --fastqs=$fastq_path --sample=E95"
