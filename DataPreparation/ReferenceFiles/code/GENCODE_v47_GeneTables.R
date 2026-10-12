library(rtracklayer)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# folder with the downloaded raw data (not part of this repository, see Data/RawDataInfo.txt)
rawDataPath <- "path/to/raw_data/"

# data folder
dataPath = paste0( home_directory, "Data/ReferenceFiles/GENCODE_v47/" )


# read the gene annotation of GENCODE release 47
gtf <- import( paste0( rawDataPath, "GENCODE_v47/gencode.v47.annotation.gtf.gz" ) )

# assignment of transcripts to genes (ids with version number), needed by tximport
transcripts = gtf[ gtf$type == "transcript" ]
tx2gene = data.frame( TXNAME = transcripts$transcript_id, GENEID = transcripts$gene_id )
fwrite( tx2gene, file = paste0( dataPath, "GENCODE_v47_tx2gene.txt" ), sep = "\t" )

# protein-coding genes (ids without version number)
allGenes = gtf[ gtf$type == "gene" ]
pcGenes = sort( unique( sub( "[.][0-9]+$", "", allGenes$gene_id[ allGenes$gene_type == "protein_coding" ] ) ) )
writeLines( pcGenes, paste0( dataPath, "GENCODE_v47_ProteinCodingGenes.txt" ) )

# the 13 protein-coding genes of the mitochondrial genome, ordered by the complexes of the respiratory chain (I, III, IV, V)
mtOrder = c( "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "MT-CYB", "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8" )
mtGenes = allGenes[ as.character( seqnames( allGenes ) ) == "chrM" & allGenes$gene_type == "protein_coding" ]
mtTab = data.frame( symbol = mtGenes$gene_name, ensg = sub( "[.][0-9]+$", "", mtGenes$gene_id ) )
mtTab = mtTab[ match( mtOrder, mtTab$symbol ), ]
fwrite( mtTab, file = paste0( dataPath, "GENCODE_v47_MitochondrialProteinCodingGenes.txt" ), sep = "\t" )
