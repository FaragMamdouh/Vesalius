# Running test for Vesalius objects
# Loading data from the packages
data(vesalius)

test_that("Dimension selection works as expected with raw embeds", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "PCA",
        dimensions = 100,
        normalisation = "log_norm")
    # testing before any image
    selected <- select_dimensions(vesalius, method = "moran")
})