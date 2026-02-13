#' Prepare DESeq2 results for MetaVolcanoR
#'
#' Converts DESeq2 results object to MetaVolcanoR format with confidence intervals
#' 
#' @param deseq_result A DESeq2 results object from results() or lfcShrink()
#' @param gene_column Character string specifying which column to use as gene 
#'        identifier. Default is "GeneID" (rownames). Can also be a column name
#'        if gene symbols are in the rowData.
#' @return A data.frame with columns: Symbol, Log2FC, pvalue, CI.L, CI.R
#' @export
#' @examples
#' \dontrun{
#' library(DESeq2)
#' dds <- DESeqDataSet(se, design = ~ condition)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' 
#' # Prepare for MetaVolcanoR
#' deg_table <- prepare_deseq2(res)
#' 
#' # Use in meta-analysis
#' meta_result <- rem_mv(list(study1 = deg_table, study2 = deg_table2))
#' }
prepare_deseq2 <- function(deseq_result, gene_column = "GeneID") {
  
  # Convert to data.frame
  df <- as.data.frame(deseq_result)
  
  # Add gene IDs from rownames
  df <- tibble::rownames_to_column(df, "GeneID")
  
  # Filter out NAs
  df <- dplyr::filter(df,
                     !is.na(log2FoldChange), 
                     !is.na(pvalue),
                     !is.na(lfcSE))
  
  # Calculate 95% confidence intervals and format
  result <- dplyr::mutate(df,
                         Symbol = .data[[gene_column]],
                         Log2FC = log2FoldChange,
                         pvalue = pvalue,
                         CI.L = log2FoldChange - 1.96 * lfcSE,
                         CI.R = log2FoldChange + 1.96 * lfcSE) %>%
           dplyr::select(Symbol, Log2FC, pvalue, CI.L, CI.R)
  
  return(result)
}

#' Prepare limma results for MetaVolcanoR
#'
#' Converts limma topTable results to MetaVolcanoR format
#' 
#' @param limma_result A data.frame from limma::topTable()
#' @param gene_column Column name containing gene identifiers (default: "GeneID")
#' @param fc_column Column name containing log fold changes (default: "logFC")
#' @param pval_column Column name containing p-values (default: "P.Value")
#' @param se_column Column name containing standard errors. If NULL, will try
#'        to calculate from CI columns or t-statistic
#' @param ci_left Column name for left CI limit (default: "CI.L")
#' @param ci_right Column name for right CI limit (default: "CI.R")
#' @return A data.frame with columns: Symbol, Log2FC, pvalue, CI.L, CI.R
#' @export
#' @examples
#' \dontrun{
#' library(limma)
#' fit <- lmFit(eset, design)
#' fit <- eBayes(fit)
#' top <- topTable(fit, number = Inf)
#' 
#' # Prepare for MetaVolcanoR
#' deg_table <- prepare_limma(top)
#' }
prepare_limma <- function(limma_result, 
                         gene_column = "GeneID",
                         fc_column = "logFC",
                         pval_column = "P.Value",
                         se_column = NULL,
                         ci_left = "CI.L",
                         ci_right = "CI.R") {
  
  df <- as.data.frame(limma_result)
  
  # Add rownames as GeneID if not present
  if(!gene_column %in% colnames(df)) {
    df <- tibble::rownames_to_column(df, gene_column)
  }
  
  # Check if CI columns exist
  has_ci <- all(c(ci_left, ci_right) %in% colnames(df))
  
  if(has_ci) {
    # Use existing CI columns
    result <- dplyr::mutate(df,
                           Symbol = .data[[gene_column]],
                           Log2FC = .data[[fc_column]],
                           pvalue = .data[[pval_column]],
                           CI.L = .data[[ci_left]],
                           CI.R = .data[[ci_right]])
  } else if(!is.null(se_column) && se_column %in% colnames(df)) {
    # Calculate CI from SE
    result <- dplyr::mutate(df,
                           Symbol = .data[[gene_column]],
                           Log2FC = .data[[fc_column]],
                           pvalue = .data[[pval_column]],
                           CI.L = .data[[fc_column]] - 1.96 * .data[[se_column]],
                           CI.R = .data[[fc_column]] + 1.96 * .data[[se_column]])
  } else {
    stop("Need either CI columns (", ci_left, ", ", ci_right, 
         ") or SE column (", se_column, ") to calculate confidence intervals.
         For basic meta-analysis without REM, use votecount_mv() or combining_mv() 
         which don't require CI.")
  }
  
  result <- dplyr::select(result, Symbol, Log2FC, pvalue, CI.L, CI.R) %>%
           dplyr::filter(!is.na(Log2FC), !is.na(pvalue))
  
  return(result)
}

#' Prepare edgeR results for MetaVolcanoR
#'
#' Converts edgeR topTags results to MetaVolcanoR format
#' 
#' @param edger_result A data.frame from edgeR::topTags()
#' @param gene_column Column name containing gene identifiers (default: uses rownames)
#' @return A data.frame with columns: Symbol, Log2FC, pvalue, CI.L, CI.R
#' @export
#' @examples
#' \dontrun{
#' library(edgeR)
#' y <- DGEList(counts = counts, group = group)
#' y <- calcNormFactors(y)
#' design <- model.matrix(~group)
#' y <- estimateDisp(y, design)
#' fit <- glmFit(y, design)
#' lrt <- glmLRT(fit)
#' top <- topTags(lrt, n = Inf)
#' 
#' # Prepare for MetaVolcanoR  
#' deg_table <- prepare_edger(top$table)
#' }
prepare_edger <- function(edger_result, gene_column = NULL) {
  
  df <- as.data.frame(edger_result)
  
  # Add gene IDs
  if(is.null(gene_column)) {
    df <- tibble::rownames_to_column(df, "GeneID")
    gene_column <- "GeneID"
  }
  
  # edgeR doesn't provide SE or CI by default
  # Calculate approximate CI from logFC and PValue
  # Using normal approximation: SE ≈ logFC / qnorm(PValue/2)
  df <- dplyr::mutate(df,
                     SE_approx = abs(logFC) / qnorm(PValue/2, lower.tail = FALSE),
                     SE_approx = ifelse(is.finite(SE_approx), SE_approx, NA))
  
  result <- dplyr::mutate(df,
                         Symbol = .data[[gene_column]],
                         Log2FC = logFC,
                         pvalue = PValue,
                         CI.L = logFC - 1.96 * SE_approx,
                         CI.R = logFC + 1.96 * SE_approx) %>%
           dplyr::select(Symbol, Log2FC, pvalue, CI.L, CI.R) %>%
           dplyr::filter(!is.na(Log2FC), !is.na(pvalue), !is.na(CI.L))
  
  return(result)
}

