################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Modifying Objects/-------------------------------#

#' Add counts to vesalius assay
#' Adding custom count matrix to a vesalius assay object.
#' @param vesalius_assay a vesalius assay object 
#' @param counts matrix or sparse matrix containing normalised counts
#' @param raw_counts matrix or sparse matrix containing raw counts
#' @param count_type character string defining the name of the count matrix
#' being added.
#' @param force logical indicating if count matrix provided should also be used
#' as raw count matrix.
#' @param verbose logical indicating if progress messages should be outputed.
#' @details In some case, you might wish to use your own normalisation method.
#' In this case, you can add your own matrix and specify the name you want to 
#' give to that count matrix. 
#' 
#' Differential gene expression tools
#' such as DESeq2 and edgeR require raw counts. As such, we recommend providing
#' raw counts as well. If you do not have raw counts, or do not wish to provide
#' them, you can set the force argument to TRUE. This will force vesalius to
#' generate a copy of your count matrix and use this as "raw" count matrix.
#' 
#' Since coordinates need to be parsed to the vesalius_assay contructor,
#' this function will also compared the your count matrix to the 
#' coordinate data. It will trim the count matrix based on barcodes shared
#' between both. 
#' 
#' @return a vesalius_assay object
#' @export

add_counts <- function(vesalius_assay,
    counts,
    raw_counts = NULL,
    count_type = "custom_counts",
    force = FALSE,
    verbose = TRUE) {
    simple_bar(verbose)
    assay <- get_assay_names(vesalius_assay)
    message_switch("add_counts", verbose, assay = assay)
    #--------------------------------------------------------------------------#
    # First lets check count validity
    # check if raw counts are already present 
    #--------------------------------------------------------------------------#
    counts <- check_counts(counts, assay, verbose)
    raw <-  ifelse(is(
        try(get_counts(vesalius_assay, type = "raw"), silent = TRUE),
        "try-error"),
      yes = FALSE,
      no = TRUE)
    #--------------------------------------------------------------------------#
    # we can do the same thing on raw counts if required
    #--------------------------------------------------------------------------#
    if (is.null(raw_counts) && !force) {
        stop("No raw counts are provided and force is set to FALSE \n
            Please provide raw counts or set force to TRUE")
    } else if (is.null(raw_counts) && force) {
        message_switch("force_count", verbose = verbose)
        raw_counts <- counts
    } else {
        raw_counts <- check_counts(raw_counts, assay, verbose)
    }
    #--------------------------------------------------------------------------#
    # Next we compare the barcodes present in tiles and adjust count
    # matrix if necesary
    #--------------------------------------------------------------------------#
    tiles <- get_tiles(vesalius_assay)
    counts <- adjust_counts(tiles, counts, verbose)
    raw_counts <- adjust_counts(tiles, raw_counts, verbose)
    #--------------------------------------------------------------------------#
    # creating count list and updating vesalius assay object
    #--------------------------------------------------------------------------#
    counts <- list(raw_counts, counts)
    if (raw) {
      message_switch("raw_count", verbose, count_type = count_type)
      names(counts) <- c(paste0("raw_", count_type), count_type)
    } else {
      names(counts) <- c("raw", count_type)
    }
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
      data = counts,
      slot = "counts",
      append = TRUE)
    #--------------------------------------------------------------------------#
    # we can update the comment on the count slot list 
    # this comment will indicate which count matrix is set as default
    #--------------------------------------------------------------------------#
    vesalius_assay <- add_active_count_tag(vesalius_assay, norm = count_type)
    #--------------------------------------------------------------------------#
    # Finally we update the vesalius commit log
    #--------------------------------------------------------------------------#
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(add_counts))
    vesalius_assay <- commit_log(vesalius_assay,
      commit,
      assay)
    simple_bar(verbose)
    return(vesalius_assay)
}


#' Add embeddings
#' Add custom embeddings to vesalius objects
#' @param vesalius_assay a vesalius_assay object
#' @param embeddings a matrix containing embedding values (see details)
#' @param embedding_type character string to be used as embedding name.
#' @param verbose logical indicating if progress message should be outputed.
#' @details Vesalius objects accepts custom embedding values that will be
#' used to generate images.
#' 
#' The embedding matrix should be in the form of a matrix with columns being
#' the latent space dimensiions and rows representing the spatial indices
#' present in the data set.
#' 
#' The intersection between spatial indices present in the tiles and
#' custom embeddings will be used.
#' This intersection will be applied to the custom embeddings and filter
#' out all tiles that are not present in the tiles. The tiles will not
#' be filtered.
#' 
#' Rownames should also be present and should represent the barcode name.
#' These rownames are use to match the latent space embedding to its tile
#' in the image.
#' 
#' @return a vesalius_assay object
#' @export
add_embeddings <- function(vesalius_assay,
    embeddings,
    embedding_type = "custom_embedding",
    verbose = TRUE) {
    simple_bar(verbose)
    assay <- get_assay_names(vesalius_assay)
    message_switch("add_embeds", verbose, assay = assay)
    #--------------------------------------------------------------------------#
    # First we check the embedding matrix to see if it is what is expected
    # filter out any barcode that does not line up with tiles
    # NOTE: Here we don't filter out tiles as we do when adding counts/tiles
    #--------------------------------------------------------------------------#
    embeddings <- check_embeddings(embeddings)
    locs <- check_barcodes(rownames(embeddings),
      unique(get_tiles(vesalius_assay)$barcodes))
    embeddings <- embeddings[locs, ]
    #--------------------------------------------------------------------------#
    # We can update the active slot which hold the default embedding that 
    # should be used 
    #--------------------------------------------------------------------------#
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
      data = embeddings,
      slot = "active",
      append = FALSE)
    #--------------------------------------------------------------------------#
    # We create a list element that will be appended to the embedding list
    # present in the vesalius_assay
    #--------------------------------------------------------------------------#
    embeddings <- list(embeddings)
    names(embeddings) <- embedding_type
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
      data = embeddings,
      slot = "embeddings",
      append = TRUE)
    #--------------------------------------------------------------------------#
    # we can update the comment on the embedding slot list 
    # this comment will indicate which embedding is actually active 
    #--------------------------------------------------------------------------#
    vesalius_assay <- add_active_embedding_tag(vesalius_assay, embedding_type)
    #--------------------------------------------------------------------------#
    # Finally we update the vesalius commit log
    #--------------------------------------------------------------------------#
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(add_embeddings))
    vesalius_assay <- commit_log(vesalius_assay,
      commit,
      assay)
    simple_bar(verbose)
    return(vesalius_assay)
}