line78 <- paste(rep("/", 78), collapse = "")

unfolder <- function(txt, rel = "include/barry/") {

    heads <- which(grepl("^\\s*#include \"", txt))

    if (!length(heads))
        return(txt)

    fns   <- gsub("^[^\"]+\"|\"$", "", txt[heads])

    head_start <- sprintf(
        "/*%s\n%1$s//\n\n Start of -%%s-\n\n%s//\n%1$s*/\n\n",
        line78,
        line78
        )

    head_end <- sprintf(
        "/*%s\n%1$s//\n\n End of -%%s-\n\n%s//\n%1$s*/\n\n",
        line78,
        line78
    )

    new_src <- txt
    for (h in rev(1L:length(fns)))
    {

        loc <- heads[h]

        fn <- paste0(rel, fns[h])
        tmp_lines <- readLines(fn, warn = FALSE)

        # Extracting relative path
        sub_rel <- gsub("[^/]+$", "", fns[h])

        tmp_lines <- unfolder(tmp_lines, rel = paste0(rel, "/", sub_rel))

        # Getting the filename
        new_src <- c(
            new_src[1:(loc - 1)],
            sprintf(head_start, fn),
            tmp_lines,
            sprintf(head_end, fn),
            new_src[(loc + 1):length(new_src)]
            )

    }

    return(new_src)

}

rel <- "include/barry/"

# Barry core
src <- readLines(paste0(rel, "barry.hpp"), warn = FALSE)
src_new <- unfolder(src, rel)
writeLines(src_new, "barry.hpp")

# # Geese core
# rel <- "include/barry/models/"
# src <- readLines(paste0(rel, "geese.hpp"), warn = FALSE)
# src_new <- unfolder(src, rel)
# writeLines(src_new, "geese.hpp")

# DEFM core
rel <- "include/barry/models/"
src <- readLines(paste0(rel, "defm.hpp"), warn = FALSE)
src_new <- unfolder(src, rel)
writeLines(src_new, "defm.hpp")
