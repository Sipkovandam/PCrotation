export DYLD_LIBRARY_PATH="/opt/intel/mkl/lib/":"/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib":$DYLD_LIBRARY_PATH
java -jar -Xmx60g /Volumes/Promise_RAID/GeneNetwork/Sipko/CorrelationLargeGenes.jar filename=/Volumes/Promise_RAID/sipko/Tests/PCA/pre_Correlation_Or_Covariance.txt writefile=/Volumes/Promise_RAID/sipko/Tests/PCA/gene_covariance.txt correlation=false
cd /Volumes/Promise_RAID/sipko/Tests/PCA/
/Volumes/Promise_RAID/juha/PCA++/pca evd /Volumes/Promise_RAID/sipko/Tests/PCA/gene_covariance.txt
mv eigenvectors.txt /Volumes/Promise_RAID/sipko/Tests/PCA/GENE.eigenvectors.txt
/Volumes/Promise_RAID/juha/PCA++/pca pc-scores covariance /Volumes/Promise_RAID/sipko/Tests/PCA/pre_Correlation_Or_Covariance.txt GENE.eigenvectors.txt
mv eigenvalues.txt /Volumes/Promise_RAID/sipko/Tests/PCA/GENE.eigenvalues.txt
gzip -f /Volumes/Promise_RAID/sipko/Tests/PCA/gene_covariance.txt
n=0
while read line; do 
if [[ $line == "Cronbach"* ]]; then continue; fi
compare=$(echo $line'<'0.7 | bc)
if [[ compare -eq 1 ]]; then break; fi
((n=$n+1))
done < cronbach.txt
echo $n
cat < GENE.eigenvectors.txt | cut -f1-$n > GENE.eigenvectors0.7.txt
gzip -f GENE.eigenvectors.txt
gzip -f /Volumes/Promise_RAID/sipko/Tests/PCA/pre_Correlation_Or_Covariance.txt
