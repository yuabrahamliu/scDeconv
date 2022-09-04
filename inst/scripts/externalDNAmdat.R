# The Infinium 27K and 450K data on human placentas were obtained from GEO 
# datasets GSE31781, GSE36829, GSE59274, GSE74738, GSE69502, GSE98224, 
# GSE125605, GSE100197, GSE75196, and GSE73375. Then, SeSAMe was used to 
# perform data preprocessing, and then the datasets were merged so that only 
# their overlapping probes shared by the Illumina 27K and 450K platforms were 
# kept. The batch difference was adjusted via ComBat (The GSE98224 methylation 
# data were used as the adjust reference).
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31781
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36829
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59274
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74738
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69502
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98224
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125605
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100197
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75196
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73375
#