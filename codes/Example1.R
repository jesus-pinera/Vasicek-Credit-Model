m <- 6
  set.seed(1)
  pds <- seq(0.4, 0.025, length.out = m)
  cota_inf <- function(pi, pj) {
    aux <- (pmax(pi + pj - 1, 0) - pi * pj) / sqrt(pi * (1 - pi) * pj * (1 -
                                                                           pj))
    return(aux)
  }
  z_lower <- outer(pds, pds, cota_inf)
  p_min <- min(z_lower)
  cota_sup <- function(pi, pj) {
    aux <- (pmin(pi, pj) - pi * pj) / sqrt(pi * (1 - pi) * pj * (1 - pj))
    return(aux)
  }
  z_upper <- outer(pds, pds, cota_sup)
  p_max <- max(z_upper)
  A <- runif(p_min / m, p_max / m, n = m ^ 2)
  dim(A) <- c(m, m)
  M <- t(A) %*% A
  M_2.7 <- M
  M_2.7
  M_cred <- list()
  n <- rep(1e5, m)
  grupos <- paste0("grupo", 1:m)
  colnames(M) <- rownames(M) <- grupos
  for (i in 1:m) {
    M_cred[[i]] <- rgamma(n[i], shape = 2, rate = 2)
  }
  V <- sapply(M_cred, sum)
  IHH <- sapply(M_cred, function(x)
    sum(x ^ 2) / sum(x) ^ 2)
  names(M_cred) <-
    names(pds) <- names(n) <- names(V) <- names(IHH) <- grupos
  start_time <- Sys.time()
  sample_aux <- corrida_unSistema(
    m = m,
    M = M,
    M_cred = M_cred,
    pds = pds,
    n = n,
    numero_muestras = 1e4,
    sim_parallel_function = Modelo_Sim.1_parallel,
    num_ciclos = 10,
    numero_nucleos = 4,
    df = Inf,
    decomp.method = "chol",
    name_modelo.sim = "Modelo_Sim.1"
  )
  end_time <- Sys.time()
  print(end_time - start_time)
  fig_2.6 <- as.data.frame(sample_aux$L_vec) %>%
    rename(
      "Grupo 1" = V1,
      "Grupo 2" = V2,
      "Grupo 3" = V3,
      "Grupo 4" = V4,
      "Grupo 5" = V5,
      "Grupo 6" = V6
    ) %>%
    reshape2::melt() %>%
    ggplot(aes(x = value, group = variable)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 100,
      color = "black",
      fill = "white"
    ) +
    geom_density(color = "red",
                 lty = 2) +
    facet_wrap(. ~ variable,
               scales = "free",
               ncol = 2) +
    theme_clean() +
    theme(
      plot.background = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "NA"
    ) +
    xlab("Pérdida") +
    ylab("") +
    scale_x_continuous(labels = scales::comma)
  fig_2.6
  M * 100
  fig_2.7A <- as.data.frame(sample_aux$L_vec[1:10000, ]) %>%
    ggplot(aes(V1, V2)) +
    geom_point() +
    geom_density2d(color = "white") +
    theme_clean() +
    theme(
      plot.background = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "NA"
    ) +
    xlab("Grupo 1") +
    ylab("Grupo 2") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma)
  fig_2.7B <- as.data.frame(sample_aux$L_vec[1:10000, ]) %>%
    ggplot(aes(V2, V5)) +
    geom_point() +
    geom_density2d(color = "white") +
    theme_clean() +
    theme(
      plot.background = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "NA"
    ) +
    xlab("Grupo 2") +
    ylab("Grupo 5") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma)
  fig_2.7C <- as.data.frame(sample_aux$L_vec[1:10000, ]) %>%
    ggplot(aes(V5, V4)) +
    geom_point() +
    geom_density2d(color = "white") +
    theme_clean() +
    theme(
      plot.background = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "NA"
    ) +
    xlab("Grupo 5") +
    ylab("Grupo 4") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma)
  fig_2.7D <- as.data.frame(sample_aux$L_vec[1:10000, ]) %>%
    ggplot(aes(V4, V2)) +
    geom_point() +
    geom_density2d(color = "white") +
    theme_clean() +
    theme(
      plot.background = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "NA"
    ) +
    xlab("Grupo 4") +
    ylab("Grupo 2") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma)
