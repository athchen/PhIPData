alias <- data.frame(alias = c("EBV", "HIV", "HPV"),
                    pattern = c("Epstein-Barr", "Human immunodeficiency virus", "Human papillomavirus"))

save(list = "alias", file = "inst/extdata/alias.rda")
