wget -O cellranger-4.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-4.0.0.tar.gz?Expires=1598301082&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci00LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE1OTgzMDEwODJ9fX1dfQ__&Signature=hg-ZPP~v8z8uj3mZ3tzn5ZhIqwYSx72QI~6kPrLl~I2h6hKEfRjvyG1UEaZj-95o3Uv4iXVmoA~maJHyFJ6Bk8mnG8TkGP-UVqxZC2bLmiRZlMcQy-TjqagPIBlPXV6PsHtxRWaa0Oh79NVWfN2RDo~YH~cFO8HEuSwJTBmZDmvvzkTgiEm1r5UpeeCQ3tnWEIUCsMURlR2g2fpa~EBRt5YRFgMkHt4haAhrUtWeC~rmEH-0ORcJSbZoLWI3YShQjaO9Od9bJF0TC9VpMx3YHFvo1yHvv-NjeADrL4BOsUV7hLT5yWkbo~3J2hl7H1DXsfWfAE-IQjq7aN67ckInjA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz

tar -xzvf cellranger-4.0.0.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

cellranger_path='.../cellranger-4.0.0'
transcriptome_path='.../refdata-gex-mm10-2020-A'
fastq_path='.../Clean_data'

## sample E95 for example
bsub -q TEST-A -n 16 -e run.err -o run.log "$cellranger_path/cellranger count --id=E95 --transcriptome=$transcriptome_path  --fastqs=$fastq_path --sample=E95"