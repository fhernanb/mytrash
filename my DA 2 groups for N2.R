DA_normal_2p <- function(m1, m2, s1, s2, n1, n2,
                         c12=1, c21=1, p1=0.5, p2=0.5, varcov='equal') {
  # To ensure vectors
  m1 <- as.vector(m1)
  m2 <- as.vector(m2)
  # Argument Verification Using Partial Matching
  varcov <- match.arg(arg=varcov, choices=c("equal","different"))
  # To obtain the log_cost_prob
  log_cost_prob <- log(c12 * p2 / (c21 * p1))
  
  if (varcov == 'equal') {
    Sp <- ((n1-1) * s1 + (n2-1) * s2) / (n1+n2-2)
    a <- matrix(m1 - m2, nrow=1) %*% solve(Sp)
    m <- 0.5 * matrix(m1 - m2, nrow=1) %*% solve(Sp) %*% matrix(m1 + m2, ncol=1)
    m <- as.numeric(m)
    res <- list(Sp=Sp, a=a, m=m, log_cost_prob=log_cost_prob, varcov=varcov)
  }
  else {
    a <- solve(s1) - solve(s2)
    b <- matrix(m1, nrow=1) %*% solve(s1) - matrix(m2, nrow=1) %*% solve(s2)
    k0 <- 0.5 * log(det(s1) / det(s2))
    k1 <- matrix(m1, nrow=1) %*% solve(s1) %*% matrix(m1, ncol=1)
    k2 <- matrix(m2, nrow=1) %*% solve(s2) %*% matrix(m2, ncol=1)
    k <- k0 + 0.5 * (k1 - k2)
    k <- as.numeric(k)
    res <- list(a=a, b=b, k=k, log_cost_prob=log_cost_prob, varcov=varcov)
  }
  
  return(res)
}

predict.DA <- function(object, data, levels=NULL) {
  
  if (is.null(levels)) levels <- 1:2
  
  data <- as.matrix(data)
  if (object$varcov == 'equal') {
    y0 <- as.vector(object$a %*% t(data))
    clasificacion <- ifelse(y0 >= object$m + object$log_cost_prob, 
                            levels[1], levels[2])
  }
  if (object$varcov == 'different') {
    part1 <- as.vector(diag(-0.5 * data %*% object$a %*% t(data)))
    part2 <- as.vector(object$b %*% t(data))
    left_side <- part1 + part2
    right_side <- object$k + object$log_cost_prob
    clasificacion <- ifelse(left_side >= right_side, 
                            levels[1], levels[2])
  }
  return(clasificacion)
}
