library("EBImage")

stabilizer <- function(input_dir, output_dir) {
  files <- system(paste("ls ", input_dir, sep = ""), intern = TRUE)
  T <- length(files)
  Data <- list()
  Sumnail <- list()
  s.size <- 64
  scale <- 512 / s.size

  for (i in 1:T) {
    Data[[i]] <- readImage(paste(dir, "/", names[i], sep = ""))
    Sumnail[[i]] <- resize(Data[[i]], s.size)
  }

  diff <- matrix(0, nrow = T - 1, ncol = 2)
  X.diff <- rep(0, T - 1)
  Y.diff <- X.diff
  R.diff <- rep(1, T - 1)

  ### Search diff by Sumnails ###
  diff.max <- 9
  x.lo <- rep(c(rep(0, 5), seq(1, 4, 1)), diff.max)
  x.up <- rep(c(seq(-4, -1, 1), rep(0, 5)), diff.max)
  x.diff <- x.lo + x.up
  y.lo <- as.vector(matrix(x.lo, nrow = diff.max, byrow = TRUE))
  y.up <- as.vector(matrix(x.up, nrow = diff.max, byrow = TRUE))
  y.diff <- y.lo + y.up
  for (t in 1:(T - 1)) {
    I1 <- Sumnail[[t]]
    I2 <- Sumnail[[t + 1]]
    R <- 1
    S <- 1
    for (i in 1:diff.max^2) {
      r <- mean((I2[(1 + x.lo[i]):(s.size + x.up[i]), (1 + y.lo[i]):(s.size + y.up[i])]
      - I1[(1 - x.up[i]):(s.size - x.lo[i]), (1 - y.up[i]):(s.size - y.lo[i])])^2)
      if (R > r) {
        R <- r
        S <- i
      }
    }
    X.diff[t] <- x.diff[S]
    Y.diff[t] <- y.diff[S]
  }

  ### Search diff in CENTER ###
  x.lo <- rep(c(rep(0, 5), seq(1, 3, 1)), 8)
  x.up <- rep(c(seq(-4, -1, 1), rep(0, 4)), 8)
  x.diff <- x.lo + x.up
  y.lo <- as.vector(matrix(x.lo, nrow = 8, byrow = TRUE))
  y.up <- as.vector(matrix(x.up, nrow = 8, byrow = TRUE))
  y.diff <- y.lo + y.up
  for (t in 1:(T - 1)) {
    I1 <- Data[[t]][c(96:160) - X.diff[t] * 4, c(96:160) - Y.diff[t] * 4]
    I2 <- Data[[t + 1]][c(96:160) + X.diff[t] * 4, c(96:160) + Y.diff[t] * 4]
    R <- 1
    S <- 1
    for (i in 1:8^2) {
      r <- mean((I2[(1 + x.lo[i]):(s.size + x.up[i]), (1 + y.lo[i]):(s.size + y.up[i])]
      - I1[(1 - x.up[i]):(s.size - x.lo[i]), (1 - y.up[i]):(s.size - y.lo[i])])^2)
      if (R > r) {
        R <- r
        S <- i
      }
    }
    X.diff[t] <- X.diff[t] * 8 + x.diff[S]
    Y.diff[t] <- Y.diff[t] * 8 + y.diff[S]
    R.diff[t] <- R
  }

  cum.s.diff <- matrix(c(0, cumsum(X.diff), 0, cumsum(Y.diff)), ncol = 2)
  plot(c(2:T), R.diff, ylim = c(0, 1), type = "l", xlab = "Frame", ylab = "Difference score")

  Diff.x <- (cum.s.diff[, 1] - min(cum.s.diff[, 1]))
  Diff.y <- (cum.s.diff[, 2] - min(cum.s.diff[, 2]))
  x.ex <- 512 + max(Diff.x)
  y.ex <- 512 + max(Diff.y)

  for (i in 1:T) {
    Im <- matrix(0, nrow = x.ex, ncol = y.ex)
    Im[((1:512) + (max(Diff.x) - Diff.x[i])), ((1:512) + (max(Diff.y) - Diff.y[i]))] <- Data[[i]]
    writeImage(Im, paste(output_dir, "/", names[i], sep = ""))
  }
}
