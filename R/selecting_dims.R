################################################################################
###############################   Vesalius      ################################
################################################################################

#---------------------/Latent space embeddings/--------------------------------#


#' select latent space dimensions
#' @param vesalius_assay a vesalius_assay object
#' @param method character string describing which method should be use to 
#' select latent space dimensions
#' @param threshold numeric - threshold of image selection
#' @details Since each latent space dimension is represented as an image, 
#' we can compute the information content of each image. If an image drops below
#' a certain threhsold, we assume that it mainly contains noise and will not
#' be uesful for territtory selection.
#' 
#' Currently, we compute the information content using entropy.
#' 
#' Vesalius will always compute information content on the active embedding.
#' Since the active embedding can be subjected to image processing steps,
#' we recommend selecting dimensions after image processing and prior to 
#' image segmentation. 
#' @returns a numeric vector of suggested dimensions to use for analysis.
#' @importFrom future.apply future_sapply
#' @export
#' @examples
#' \dontrun{
#' data(vesalius)
#' # First we build a simple object
#' ves <- build_vesalius_object(coordinates, counts)
#' # We can do a simple run 
#' ves <- generate_embeddings(ves)
#' ves <- smooth_image(ves, iter = 5, sigma = 3)
#' embeds <- select_dimensions(ves)
#'}
select_dimensions <- function(vesalius_assay,
  method = "moran",
  threshold = 0.75,
  verbose = TRUE) {
    simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # check method and convert active
    # this section can be improved to make it faster 
    # you don't need to join to tiles only to origin 
    #--------------------------------------------------------------------------#
    method <- check_dim_selection_method(method)
    tiles <- check_tiles(vesalius_assay)
    active <- get_embeddings(vesalius_assay, active = TRUE)
    message_switch("vtc", verbose)
    active <- future_lapply(seq_len(ncol(active)),
      FUN = future_ves_to_cimg,
      embeddings = active,
      dims = seq_len(ncol(active)),
      tiles = tiles,
      full_image = FALSE,
      future.seed = TRUE)
    active <- future_lapply(active, FUN = simplify_image)
    #--------------------------------------------------------------------------#
    # Calculate information content of image
    # Using a switch for now. I might add more methods later 
    #--------------------------------------------------------------------------#
    message_switch("info", verbose, method = method)
    information_content <- switch(EXPR = method,
      "altieri" = sapply(active, ves_altieri),
      "batty" = sapply(active, ves_batty),
      "contagion" = sapply(active, ves_contagion),
      "karlstrom" = sapply(active, ves_karlstrom),
      "leibovici" = sapply(active, ves_leibovici),
      "oneill" = sapply(active, ves_oneil),
      "parredw" = sapply(active, ves_parredw),
      "shannon" = sapply(active, ves_shannon),
      "shannonZ" = sapply(active, ves_shannonZ),
      "moran" = sapply(active, ves_moran))#, future.seed = TRUE))
    #--------------------------------------------------------------------------#
    # No need to update the object for now 
    # I could add this to meta paramters or something? 
    # Remebering paramters or something along those lines
    #--------------------------------------------------------------------------#
    information_content <- which(information_content >= 
        quantile(information_content, threshold))
    simple_bar(verbose)
    return(information_content)
}

#' simplify images before computing entropy
#' @param image a cimg df
#' @returns a cimg df
#' @importFrom stats quantile
#' 
simplify_image <- function(image) {
    #quant <- stats::quantile(image$value, threshold)
    #image[image$value <= quant, "value"] <- 0
    #image[image$value > quant, "value"] <- 1
    image <- image %>% filter(origin == 1)
    return(image)
}

# #------------------------/ Selecting Dimensions /-----------------------------#
# #' compute Altieri spatial entropy
# #' @param active data frame containing pixel values as cimg df
# #' @returns image entropy
# #' @importFrom SpatEntropy altieri
# #' @importFrom spatstat.geom ppp owin
# ves_altieri <- function(active) {
#   x_range <- c(min(active$x), max(active$x))
#   y_range <- c(min(active$y), max(active$y))
#   ppp_obj <- spatstat.geom::ppp(x =  active$x,
#     y = active$y,
#     window = spatstat.geom::owin(x_range, y_range),
#     marks = active$value)
#   entropy <- SpatEntropy::altieri(ppp_obj)
#   return(entropy$RES)
# }

# #' compute batty spatial entropy
# #' @param active data frame containing pixel values as cimg df
# #' @returns image entropy
# #' @importFrom SpatEntropy batty
# #' @importFrom spatstat.geom ppp owin
# ves_batty <- function(active) {
#   x_range <- c(min(active$x), max(active$x))
#   y_range <- c(min(active$y), max(active$y))
#   ppp_obj <- spatstat.geom::ppp(x =  active$x,
#     y = active$y,
#     window = spatstat.geom::owin(x_range, y_range))
#   entropy <- SpatEntropy::batty(ppp_obj, partition = 100)
#   return(entropy$batty)
# }

# #' compute relative contagion index 
# #' @param active data frame containing pixel values as cimg df
# #' @returns contagion index
# #' @importFrom SpatEntropy contagion
# ves_contagion <- function(active) {
#   dat <- as.matrix(imager::as.cimg(active[, c("x", "y", "value")]))
#   entropy <- SpatEntropy::contagion(dat)
#   return(entropy$contagion)
# }

# #' compute karlstrom entropy
# #' @param active data frame containing pixel values as cimg df
# #' @returns entropy
# #' @importFrom SpatEntropy karlstrom
# ves_karlstrom <- function(active) {
#   x_range <- c(min(active$x), max(active$x))
#   y_range <- c(min(active$y), max(active$y))
#   ppp_obj <- spatstat.geom::ppp(x =  active$x,
#     y = active$y,
#     window = spatstat.geom::owin(x_range, y_range))
#   entropy <- SpatEntropy::karlstrom(ppp_obj, neigh = 6, partition = 100)
#   return(entropy$karlstrom)
# }

# #' compute leibovici entropy
# #' @param active data frame containing pixel values as cimg df
# #' @returns entropy 
# #' @importFrom SpatEntropy leibovici
# ves_leibovici <- function(active) {
#   x_range <- c(min(active$x), max(active$x))
#   y_range <- c(min(active$y), max(active$y))
#   ccdist <- 0.05 * mean(c(x_range[2] - x_range[1],
#     y_range[2] - y_range[1]))
#   ppp_obj <- spatstat.geom::ppp(x =  active$x,
#     y = active$y,
#     window = spatstat.geom::owin(x_range, y_range))
#   entropy <- SpatEntropy::leibovici(ppp_obj, ccdist = ccdist, partition = 100)
#   return(entropy$leib)
# }

# #' compute oneil entropy
# #' @param active data frame containing pixel values as cimg df
# #' @returns image entropy
# #' @importFrom SpatEntropy oneill
# ves_oneil <- function(active) {
#   dat <- as.matrix(imager::as.cimg(active[, c("x", "y", "value")]))
#   entropy <- SpatEntropy::oneill(dat)
#   return(entropy$oneill)
# }

# #' compute parredw entropy
# #' @param active data frame containing pixel values as cimg df
# #' @returns image entropy
# #' @importFrom SpatEntropy parredw
# ves_parredw <- function(active) {
#   dat <- as.matrix(imager::as.cimg(active[, c("x", "y", "value")]))
#   entropy <- SpatEntropy::parredw(dat)
#   return(entropy$parredw) 
# }

# #' compute shannow entropy
# #' @param active data frame containing pixel values as cimg df
# #' @returns image entropy
# #' @importFrom SpatEntropy shannon
# ves_shannon <- function(active) {
#   x_range <- c(min(active$x), max(active$x))
#   y_range <- c(min(active$y), max(active$y))
#   ppp_obj <- spatstat.geom::ppp(x =  active$x,
#     y = active$y,
#     marks = active$value,
#     window = spatstat.geom::owin(x_range, y_range))
#   entropy <- SpatEntropy::shannon(ppp_obj)
#   return(entropy$shann)
# }

# #' compute shannonZ entropy
# #' @param active data frame containing pixel values as cimg df
# #' @returns image entropy
# #' @importFrom SpatEntropy shannonZ
# ves_shannonZ <- function(active) {
#   x_range <- c(min(active$x), max(active$x))
#   y_range <- c(min(active$y), max(active$y))
#   ppp_obj <- spatstat.geom::ppp(x =  active$x,
#     y = active$y,
#     marks = active$value,
#     window = spatstat.geom::owin(x_range, y_range))
#   entropy <- SpatEntropy::shannonZ(ppp_obj)
#   return(entropy$shannZ)
# }


#' compute spatial auto-correlation on image 
#' @param active a cimg data frame 
#' @return autocorrelation
#' @importFrom ape Moran.I
#' 
ves_moran <- function(active) {
    inverse_distance <- active %>% 
        select(x, y) %>%
        dist() %>%
        as.matrix()
    inverse_distance <- 1 / inverse_distance
    diag(inverse_distance) <- 0
    # should return pval as well at least for filtering 
    moran <- Moran.I(active$value, inverse_distance)$observed
    return(moran)
}
