library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
outputPath = paste0( home_directory, "Figure_3/" )
tablePath = paste0( outputPath, "tables/" )


# overview of the donor covariates (age, PMI, sex, brain bank, batches) in the disease groups and the HC group, for each cell state

# table of the clustered OPCs and oligodendrocytes (03_OPC_OL_Clustering_Analysis.R), only cells of the kept clusters
cells = as.data.frame( fread( paste0( tablePath, "OPC_OL_Clustering/OPC_OL_CellTable.txt.gz" ) ) )
cells <- cells[ cells$final_retained, ]

# donor information (columns: sample_id, Diagnosis_author, Diagnosis, Sex, Age, PMI, brain_bank, extract, flowcell, disease_duration, subtype)
# made by DataPreparation/snRNAseq_HC_MSA_PD/code/DonorCovariates_AuthorMetadata.R from the data of Rydbirk et al. (2026)
donor_info = as.data.frame( fread( paste0( dataPath, "snRNAseq_HC_MSA_PD/NatComm_DonorCovariates.txt" ) ) )

states = c( "OPC", paste0( "OL", 1:6 ) )
comparisons = c( MSA_vs_HC = "MSA", PD_vs_HC = "PD" )

# a donor is used for a cell state if it has at least 10 cells of this state, each group needs at least 3 donors
minCells = 10
minDonors = 3

numericVars = data.frame( variable = c( "Age", "PMI" ), label = c( "Age (years)", "PMI (hours)" ), author_precedent = c( TRUE, FALSE ) )
categoricalVars = data.frame( variable = c( "Sex", "brain_bank", "extract", "flowcell" ), label = c( "Sex", "Brain bank", "Extract", "Flowcell" ) )


# standardized mean difference between two groups (second minus first group)
computeSMD <- function( x, g )
{
    x1 = x[ g == levels( g )[ 1 ] ]
    x2 = x[ g == levels( g )[ 2 ] ]
    sdPooled = sqrt( ( sd( x1 )^2 + sd( x2 )^2 ) / 2 )
    if( !is.finite( sdPooled ) || sdPooled == 0 ) return( NA )
    ( mean( x2 ) - mean( x1 ) ) / sdPooled
}

# Cramer's V of a contingency table
computeCramersV <- function( tab )
{
    if( nrow( tab ) < 2 || ncol( tab ) < 2 || sum( tab ) == 0 ) return( NA )
    chi = suppressWarnings( chisq.test( tab, correct = FALSE ) )
    k = min( nrow( tab ) - 1, ncol( tab ) - 1 )
    sqrt( as.numeric( chi$statistic ) / ( sum( tab ) * k ) )
}


supportTab = c()
numericTab = c()
categoricalTab = c()
levelTab = c()

for( state in states )
{
    cellsPerDonor = table( cells$sample_id[ cells$final_cluster_label == state ] )
    donor_info$n_cells = as.numeric( cellsPerDonor[ donor_info$sample_id ] )
    donor_info$n_cells[ is.na( donor_info$n_cells ) ] = 0

    for( comp in names( comparisons ) )
    {
        disease = comparisons[ comp ]
        d = donor_info[ donor_info$Diagnosis %in% c( "HC", disease ) & donor_info$n_cells >= minCells, ]
        d$group = factor( d$Diagnosis, levels = c( "HC", disease ) )
        n_HC = sum( d$group == "HC" )
        n_case = sum( d$group == disease )
        supported = n_HC >= minDonors & n_case >= minDonors

        supportTab = rbind( supportTab, data.frame( natural_label = state, contrast = comp, case = disease, min_cells_per_sample = minCells, min_samples_per_group = minDonors,
                                                    n_HC = n_HC, n_case = n_case, n_samples = nrow( d ), supported_for_DE = supported ) )
        if( !supported ) next

        # age and PMI
        for( i in 1:nrow( numericVars ) )
        {
            x = as.numeric( d[ , numericVars$variable[ i ] ] )
            keep = is.finite( x )
            x = x[ keep ]
            g = droplevels( d$group[ keep ] )

            smd = computeSMD( x, g )
            pval = wilcox.test( x ~ g, exact = FALSE )$p.value

            # classes by the absolute SMD: 0.8 and 0.5
            balance = if( !is.finite( smd ) ) "not_testable" else if( abs( smd ) >= 0.8 ) "major" else if( abs( smd ) >= 0.5 ) "moderate" else "minor"
            action = if( balance == "major" ) "main_model_adjust" else if( balance == "moderate" ) "sensitivity_only" else "do_not_adjust"

            numericTab = rbind( numericTab, data.frame( natural_label = state, contrast = comp, variable = numericVars$variable[ i ], variable_label = numericVars$label[ i ],
                                                        variable_type = "numeric", author_precedent = numericVars$author_precedent[ i ],
                                                        group1 = levels( g )[ 1 ], n1 = sum( g == levels( g )[ 1 ] ), mean1 = mean( x[ g == levels( g )[ 1 ] ] ), median1 = median( x[ g == levels( g )[ 1 ] ] ),
                                                        group2 = levels( g )[ 2 ], n2 = sum( g == levels( g )[ 2 ] ), mean2 = mean( x[ g == levels( g )[ 2 ] ] ), median2 = median( x[ g == levels( g )[ 2 ] ] ),
                                                        smd = smd, abs_smd = abs( smd ), wilcox_p = pval, balance_class = balance, final_action = action ) )
        }

        # sex, brain bank and the two batch variables (RNA extraction and flowcell)
        for( i in 1:nrow( categoricalVars ) )
        {
            v = d[ , categoricalVars$variable[ i ] ]
            keep = !is.na( v )
            tab = table( d$group[ keep ], v[ keep ] )

            pval = tryCatch( fisher.test( tab )$p.value, error = function( e ) NA )
            cramersV = computeCramersV( tab )

            # levels that occur in both groups and levels that occur in one group only
            shared = colnames( tab )[ colSums( tab > 0 ) == nrow( tab ) ]
            exclusive = sapply( 1:nrow( tab ), function( j ) {
              lev = colnames( tab )[ tab[ j, ] > 0 & colSums( tab > 0 ) == 1 ]
              if( length( lev ) > 0 ) paste( lev, collapse = ", " ) else "none" } )

            # classes by the number of shared levels and Cramer's V: 0.5 and 0.3
            balance = if( length( shared ) == 0 ) "structural_confounding" else if( !is.finite( cramersV ) ) "minor" else if( cramersV >= 0.5 ) "major" else if( cramersV >= 0.3 ) "moderate" else "minor"
            action = if( balance == "structural_confounding" ) "structurally_confounded_do_not_adjust" else if( balance %in% c( "major", "moderate" ) ) "sensitivity_only" else "do_not_adjust"

            categoricalTab = rbind( categoricalTab, data.frame( natural_label = state, contrast = comp, variable = categoricalVars$variable[ i ], variable_label = categoricalVars$label[ i ],
                                                                variable_type = "categorical", author_precedent = FALSE, fisher_p = pval, cramers_v = cramersV,
                                                                shared_levels_n = length( shared ), shared_levels = if( length( shared ) > 0 ) paste( shared, collapse = ", " ) else NA,
                                                                exclusive_levels = paste( rownames( tab ), exclusive, sep = ": ", collapse = " | " ),
                                                                balance_class = balance, final_action = action ) )

            lev = as.data.frame( tab )
            colnames( lev ) = c( "Diagnosis", "level", "Freq" )
            lev$natural_label = state
            lev$contrast = comp
            lev$variable = categoricalVars$variable[ i ]
            lev$variable_label = categoricalVars$label[ i ]
            lev$prop = lev$Freq / sapply( lev$Diagnosis, function( dx ) sum( lev$Freq[ lev$Diagnosis == dx ] ) )
            levelTab = rbind( levelTab, lev )
        }
    }
}

# p-values adjusted within each cell state and comparison
numericTab$wilcox_padj = NA
categoricalTab$fisher_padj = NA
for( state in states )
{
    for( comp in names( comparisons ) )
    {
        cand = which( numericTab$natural_label == state & numericTab$contrast == comp )
        numericTab$wilcox_padj[ cand ] = p.adjust( numericTab$wilcox_p[ cand ], method = "BH" )
        cand = which( categoricalTab$natural_label == state & categoricalTab$contrast == comp )
        categoricalTab$fisher_padj[ cand ] = p.adjust( categoricalTab$fisher_p[ cand ], method = "BH" )
    }
}

# short overview of all variables
decisionTab = rbind( data.frame( numericTab[ , c( "natural_label", "contrast", "variable", "variable_label", "variable_type" ) ],
                                 metric_label = paste0( "SMD ", sprintf( "%.2f", numericTab$abs_smd ) ),
                                 numericTab[ , c( "balance_class", "final_action" ) ] ),
                     data.frame( categoricalTab[ , c( "natural_label", "contrast", "variable", "variable_label", "variable_type" ) ],
                                 metric_label = ifelse( categoricalTab$shared_levels_n == 0, "no overlap", paste0( "V ", sprintf( "%.2f", categoricalTab$cramers_v ) ) ),
                                 categoricalTab[ , c( "balance_class", "final_action" ) ] ) )

fwrite( supportTab, file = paste0( tablePath, "DonorCovariates/DonorCovariates_DonorsPerStateAndComparison.txt" ), sep = "\t" )
fwrite( numericTab, file = paste0( tablePath, "DonorCovariates/DonorCovariates_Balance_Age_PMI.txt" ), sep = "\t" )
fwrite( categoricalTab, file = paste0( tablePath, "DonorCovariates/DonorCovariates_Balance_Sex_BrainBank_Batches.txt" ), sep = "\t" )
fwrite( levelTab, file = paste0( tablePath, "DonorCovariates/DonorCovariates_Balance_Sex_BrainBank_Batches_Levels.txt" ), sep = "\t" )
fwrite( decisionTab, file = paste0( tablePath, "DonorCovariates/DonorCovariates_Balance_Overview.txt" ), sep = "\t" )
