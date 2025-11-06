ntop=50
ngene=2000
out_name =  paste0("paper_ntop=",ntop,"_ngene=",ngene,".pdf")
rmarkdown::render("paper/paper.Rmd",knit_root_dir = "..",
                  params = list(ntop = ntop,
                                ngene=ngene),
                  output_dir = "paper",
                  output_file = out_name)



trackdown::update_file(
  "paper/paper.Rmd",           # path to your .Rmd / .qmd / .Rnw
  gfile = "paper",   # optional: the name it will have in Google Drive
  gpath = "desurv_paper_trackdown",    # folder in Drive (will create if missing)
  hide_code = FALSE,      # hide code + header from collaborators?
  path_output = "paper/paper_ntop=50_ngene=2000.pdf"     # also upload the rendered PDF or HTML alongside
)

trackdown::update_file(
  "paper/04_results.Rmd",           # path to your .Rmd / .qmd / .Rnw
  gfile = "results",   # optional: the name it will have in Google Drive
  gpath = "desurv_paper_trackdown",    # folder in Drive (will create if missing)
  hide_code = FALSE      # hide code + header from collaborators?
)

trackdown::update_file(
  "paper/03_methods.Rmd",           # path to your .Rmd / .qmd / .Rnw
  gfile = "methods",   # optional: the name it will have in Google Drive
  gpath = "desurv_paper_trackdown",    # folder in Drive (will create if missing)
  hide_code = FALSE      # hide code + header from collaborators?
)

trackdown::upload_file(
  "paper/02_introduction_30102025.Rmd",           # path to your .Rmd / .qmd / .Rnw
  gfile = "introduction",   # optional: the name it will have in Google Drive
  gpath = "desurv_paper_trackdown",    # folder in Drive (will create if missing)
  hide_code = FALSE      # hide code + header from collaborators?
)

trackdown::update_file(
  "paper/05_discussion.Rmd",           # path to your .Rmd / .qmd / .Rnw
  gfile = "discussion",   # optional: the name it will have in Google Drive
  gpath = "desurv_paper_trackdown",    # folder in Drive (will create if missing)
  hide_code = FALSE      # hide code + header from collaborators?
)
