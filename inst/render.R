out_name =  paste0("paper_ntop=",50,"_ngene=",2000,".pdf")
    rmarkdown::render("paper/paper.Rmd",knit_root_dir = "..",
                      params = list(ntop = 50,
                                    ngene=2000),
                      output_dir = "paper",
                      output_file = out_name)
    