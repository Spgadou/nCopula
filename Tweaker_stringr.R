assignInNamespace("fix_replacement",
                  function(x) {
                    vapply(x, stringr:::fix_replacement_one, character(1), USE.NAMES = FALSE)
                  },
                  asNamespace("stringr"))

assignInNamespace("fix_replacement_one",
                  function(x){
                    if (is.na(x)) {
                      return(x)
                    }

                    chars <- stringr:::str_split(x, "")[[1]]
                    out <- character(length(chars))
                    escaped <- logical(length(chars))

                    in_escape <- FALSE
                    for (i in seq_along(chars)) {
                      escaped[[i]] <- in_escape
                      char <- chars[[i]]

                      if (in_escape) {
                        # Escape character not printed previously so must include here
                        if (char == "$") {
                          out[[i]] <- "\\\\$"
                        } else if (char >= "0" && char <= "9") {
                          out[[i]] <- paste0("$", char)
                        } else {
                          out[[i]] <- paste0("\\", char)
                        }

                        in_escape <- FALSE
                      } else {
                        if (char == "$") {
                          out[[i]] <- "\\$"
                        } else if (char == "\\") {
                          in_escape <- TRUE
                        } else {
                          out[[i]] <- char
                        }
                      }
                    }

                    # tibble::tibble(chars, out, escaped)
                    paste0(out, collapse = "")
                  },
                  asNamespace("stringr"))



