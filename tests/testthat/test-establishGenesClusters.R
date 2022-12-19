
test_that("Establish genes clusters", {
  skip("Test not ready yet")

  data("raw.dataset")
  objCOTAN <- COTAN(raw = raw.dataset)
  objCOTAN <- proceedToCoex(objCOTAN, cores = 12, saveObj = FALSE)

  #primaryMarkers <- getGenes(objCOTAN)[sample(getNumGenes(objCOTAN), 10)]
  primaryMarkers <- c("Pcbp2", "Snrpe", "Nfyb", "Prpf40a", "Ergic2",
                      "Ncl", "Cd47", "Macrod2", "Fth1", "Supt16")

  df <- get.gene.coexpression.space(object = as(objCOTAN, "scCOTAN"),
                                    n.genes.for.marker = 11,
                                    primary.markers = primaryMarkers)

})
