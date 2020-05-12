
ans <- list()
N <- 3
M <- 4

# S <- t(combn(1:(M*N), 2))
pset <- function(i = 0, m = NULL) {
  
  if (is.null(m)) {
    m <- matrix(0, N, M)
    ans <<- c(ans, list(m))
  }
  
  m1    <- m
  S     <- cbind(i %% N, i %/% N) + 1
  if (m1[S] == 1)
    stop("Should be 0!")
  m1[S] <- 1
  
  ans <<- c(ans, list(m1))
  
  if (i < (N*M - 1)) {
    pset(i + 1, m)
    pset(i + 1, m1)
  }
  
}

pset(0)

# table(sapply(ans, paste, collapse= ""))
