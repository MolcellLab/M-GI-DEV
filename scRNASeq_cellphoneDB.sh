source activate cpdb

cellphonedb method statistical_analysis meta.txt data.txt --subsampling --subsampling-log true --threads=10 --counts-data=gene_name
## dotplot
cut -f2 ./out/significant_means.txt |sed 1d > ./out/rows.txt
cellphonedb plot dot_plot --means-path ./out/means.txt --pvalues-path ./out/pvalues.txt --output-path ./out --rows ./out/rows.txt --output-name dotsig.pdf
cellphonedb plot heatmap_plot meta.txt --pvalues-path ./out/pvalues.txt --output-path ./out