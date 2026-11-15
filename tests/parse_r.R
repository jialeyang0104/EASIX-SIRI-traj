files <- list.files("src",pattern="[.]R$",recursive=TRUE,full.names=TRUE)
for (f in files) parse(f)
cat("PASS R syntax:",length(files),"files\n")
