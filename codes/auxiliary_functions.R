library(mvtnorm)
library(Matrix)

f_cor <- function(rho, pi,pj){
  ui <- qnorm(pi)
  uj <- qnorm(pj)
  return(mvtnorm::pmvnorm(upper = c(ui,uj), corr = rbind(c(1,rho), c(rho,1)))[1])
}


g_cor <- function(rho, params){
  pi <- params[1]
  pj <- params[2]
  rij <- params[3]
  aux <- f_cor(rho, pi, pj) - (pi * pj + rij*sqrt( pi*(1-pi) * pj*(1-pj)))
  return(aux)
}

g_cor <- Vectorize(g_cor, "rho")



rbinom_V_p <- Vectorize(rbinom, "prob")
sample_V <- Vectorize(dqrng::dqsample.int, "size")


Modelo_Sim.1 <-
  function(nS = 1e4, semilla = 1, pds, n, M, M_cred, tol.raiz = 1e-16, decomp.method = "eig"){
    
  m <- length(pds) ### numero de grupos
  u <- qnorm(as.numeric(pds)) #u's ### Umbral de variables latentes
  rhostildeij <- rep(0,m^2)
  dim(rhostildeij) <- c(m,m)
  for(i in 1:m) for(j in 1:m){ 
        rhostildeij[i,j] <- uniroot(g_cor, 
                                     interval = c(-1,1), 
                                     params = c(pds[i],pds[j],M[i,j]),
                                     tol = tol.raiz)$root
      }
      
     
  
  #rhosijtilde # Matriz de correlaciones extragrupo y diagonal intragrupo calibradas
  
  rhostilde <- as.numeric(diag(rhostildeij))
  rhostilde <- ifelse(abs(rhostilde)>1e-12, rhostilde, rep(0,m))
  
  aux <- matrix(1, nrow = m, ncol = m)
  aux[upper.tri(aux)] <- rhostildeij[upper.tri(rhostildeij)]
  aux[lower.tri(aux)] = t(aux)[lower.tri(aux)]
  diag(aux) <- rhostilde
  
  rhostildeij <- aux
  
  w_aux <-  rhostildeij/sqrt(rhostilde %*% t(rhostilde))  
  
  
  
  
  if(!matrixcalc::is.positive.semi.definite(w_aux)){
    w <- as.matrix(nearPD(w_aux, keepDiag = TRUE, eig.tol = 1e-16,
                        conv.tol = 1e-15,
                        conv.norm.type = "F",
                        maxit = 1000)$mat)
  }else{
    w <- w_aux
  }
  
  Y_s <- mvtnorm::rmvnorm(nS, sigma = w, method = decomp.method)
  
  r_s <- sqrt(rhostilde)
  
  K <- L <- matrix(NA, nrow = nS, ncol = m)
  
  incum_nombres <- list()
  
  
  grupos <- colnames(M)
  
  for (i in 1:length(u)) {
    p_z <- pnorm((u[i] - r_s[i]*Y_s[,i]) / (sqrt(1-r_s[i]^2)) )
    
    K[,i] <- rbinom_V_p(1, n[i], p_z)
    
    indx <- sample_V(n = n[i], size = K[,i], replace = FALSE)
    
    aux1 <- list()
    
    for (j in 1:length(indx)) {
      
      L[j,i] <- sum(M_cred[[grupos[i]]][indx[[j]]])
      
    }
  
  }
  
  return(list(L = L, K = K, Z = Y_s, omega_ori = w_aux, omega_aj = w))
  }


Modelo_Sim.1_parallel <- function(nS) {
  Modelo_Sim.1(nS, pds = pds, M = M, 
          n = n, M_cred = M_cred, 
          tol.raiz = 1e-12, 
          decomp.method = method)
}

Modelo_Sim.2_cop.t <-
  function(nS = 1e4, semilla = 1, pds, n, M, M_cred, tol.raiz = 1e-16, decomp.method = "eig"){
    
    m <- length(pds) ### numero de grupos
    u <- qnorm(as.numeric(pds)) #u's ### Umbral de variables latentes
    rhostildeij <- rep(0,m^2)
    dim(rhostildeij) <- c(m,m)
    for(i in 1:m) for(j in 1:m){ 
      rhostildeij[i,j] <- uniroot(g_cor, 
                                  interval = c(-1,1), 
                                  params = c(pds[i],pds[j],M[i,j]),
                                  tol = tol.raiz)$root
    }
    
    
    
    #rhosijtilde # Matriz de correlaciones extragrupo y diagonal intragrupo calibradas
    
    rhostilde <- as.numeric(diag(rhostildeij))
    aux <- matrix(1, nrow = m, ncol = m)
    aux[upper.tri(aux)] <- rhostildeij[upper.tri(rhostildeij)]
    aux[lower.tri(aux)] = t(aux)[lower.tri(aux)]
    diag(aux) <- diag(rhostildeij)
    rhostildeij <- aux
    
    
    w_aux <-  rhostildeij/sqrt(rhostilde %*% t(rhostilde))
    
    if(!matrixcalc::is.positive.semi.definite(w_aux)){
      w <- as.matrix(nearPD(w_aux, keepDiag = TRUE, eig.tol = 1e-16,
                            conv.tol = 1e-15,
                            conv.norm.type = "F",
                            maxit = 1000)$mat)
    }else{
      w <- w_aux
    }
    
    Y_s <- mvtnorm::rmvt(nS, sigma = w, method = decomp.method,
                         df = df)
    
    r_s <- sqrt(rhostilde)
    
    K <- L <- matrix(NA, nrow = nS, ncol = m)
    
    incum_nombres <- list()
    
    
    grupos <- colnames(M)
    
    for (i in 1:length(u)) {
      p_z <- pnorm((u[i] - r_s[i]*Y_s[,i]) / (sqrt(1-r_s[i]^2)) )
      
      K[,i] <- rbinom_V_p(1, n[i], p_z)
      
      indx <- sample_V(n = n[i], size = K[,i], replace = FALSE)
      
      aux1 <- list()
      
      for (j in 1:length(indx)) {
        
        L[j,i] <- sum(M_cred[[grupos[i]]][indx[[j]]])
        
      }
      
    }
    
    return(list(L = L, K = K, Z = Y_s,omega_ori = w_aux, omega_aj = w))
  }

Modelo_Sim.2_cop.t_un.sector <-
  function(nS = 1e4, semilla = 1, pds, n, r, M_cred, tol.raiz = 1e-16, decomp.method = "eig"){
    
    m <- 1
    u <- qnorm(as.numeric(pds)) #u's ### Umbral de variables latentes
    rhostilde <- uniroot(g_cor, 
                         interval = c(-1,1), 
                         params = c(pds, pds, r),
                         tol = tol.raiz)$root
    
    #rhosijtilde # Matriz de correlaciones extragrupo y diagonal intragrupo calibradas
    
    
    Y_s <- rnorm(nS)
    
    r_s <- sqrt(rhostilde)
    
    K <- L <- matrix(NA, nrow = nS, ncol = m)
    
    incum_nombres <- list()
    
    
    grupos <- colnames(M)
    
    for (i in 1:length(u)) {
      p_z <- pnorm((u[i] - r_s[i]*Y_s[,i]) / (sqrt(1-r_s[i]^2)) )
      
      K[,i] <- rbinom_V_p(1, n[i], p_z)
      
      indx <- sample_V(n = n[i], size = K[,i], replace = FALSE)
      
      aux1 <- list()
      
      for (j in 1:length(indx)) {
        
        L[j,i] <- sum(M_cred[[grupos[i]]][indx[[j]]])
        
      }
      
    }
    
    return(list(L = L, K = K, Z = Y_s,omega_ori = w_aux, omega_aj = w))
  }


Modelo_Sim.2_parallel <- function(nS) {
  Modelo_Sim.2_cop.t(nS, pds = pds, M = M,
                     n = n, M_cred = M_cred, 
                     tol.raiz = 1e-12, 
                     decomp.method = method,
                     df)
}

corrida_unSistema <- function(m,
                              M, M_cred, pds, n,
                              df,
                              decomp.method,
                              numero_muestras,
                              sim_parallel_function,
                              name_modelo.sim,
                              num_ciclos,
                              numero_nucleos = 8){
  
  
  grupos <- names(pds)
  
  Maux <- M * 0
  
  for(i in grupos) for(j in grupos[match(i,grupos):m]){
    k <- match(i,grupos)
    l <- match(j,grupos)
    Maux[k,l] <- ifelse(k==l, V[k]^2*pds[k]*(1-pds[k])*(M[k,k] + (1-M[k,k])*IHH[k]),
                        V[k]*V[l]*M[k,l]*sqrt(pds[k]*(1-pds[k])*pds[l]*(1-pds[l])))
    Maux[l,k] <- Maux[k,l]
  }
  
  Lbar <- pds[grupos]*V[grupos]
  
  method <- decomp.method
  nS <- numero_muestras
  
  cl <- makeCluster(numero_nucleos)  
  
  registerDoParallel(cl) 
  clusterExport(cl=cl, varlist = list("pds", "M", "n", "M_cred",
                                      "df",
                                      "method", "g_cor", "f_cor", 
                                      name_modelo.sim,
                                      "rbinom_V_p",
                                      "sample_V",
                                      "uniroot",
                                      "pmvnorm",
                                      "nearPD"),
                envir=environment())
  
  sample_L <- 
    sample_K <- 
    sample_Z <- list()
  
  num_part <- 8
 
  for(i in 1:num_ciclos){
    aux <- parLapply(cl, X = rep(nS/num_part,num_part) , sim_parallel_function)
    
    
    aux2 <- comprehenr::to_list(for(i in 1:num_part) rbind.data.frame(aux[[i]]$L))
    aux2.1 <- comprehenr::to_list(for(i in 1:num_part) rbind.data.frame(aux[[i]]$K))
    aux2.2 <- comprehenr::to_list(for(i in 1:num_part) rbind.data.frame(aux[[i]]$Z))
    aux3 <- data.table::rbindlist(aux2)
    aux3.1 <- data.table::rbindlist(aux2.1)
    aux3.2 <- data.table::rbindlist(aux2.2)
    sample_L[[i]] <- as.data.frame(aux3)
    sample_K[[i]] <- as.data.frame(aux3.1)
    sample_Z[[i]] <- as.data.frame(aux3.2)
    
    
    rm(aux2, aux2.1, aux2.2, aux3, aux3.1, aux3.2)
    
    print(paste("Ciclo", i, "de", num_ciclos))
  }
  
  omega_ori <- aux[[2]]$omega_ori
  omega_aj <- aux[[2]]$omega_aj
  stopCluster(cl)
  rm(aux)
  
  sample_L_aux <- sample_L
  
  sample_L <- as.data.frame(data.table::rbindlist(sample_L_aux))
  

  sample_K_aux <- sample_K
  
  sample_K <- as.data.frame(data.table::rbindlist(sample_K_aux))
  
  sample_Z_aux <- sample_Z
  
  sample_Z <- as.data.frame(data.table::rbindlist(sample_Z_aux))
  
  
 
  return(list(V = V,
              L_vec = sample_L,
              K_vec = sample_K,
              Z_vec = sample_Z,
              Maux = Maux,
              pds = pds,
              omega_aj = omega_aj,
              omega_ori = omega_ori))
  
}


obtener_cuantiles_torres <- function(u,alfa,h, X){
  u_dir <- as.matrix(u)
  
  norm_u <- norm(u_dir, type = "2")
  u_dir_norm <- u_dir/norm_u
  norm(u_dir_norm, "2")
  
  n <- length(u_dir)
  
  e_n <- as.matrix(1/sqrt(n)*rep(1,n))
  I_n <- diag(n)
  M_u <- cbind(u_dir_norm, sign(u_dir_norm)[-1] * I_n[,-1])
  M_e <- cbind(e_n, I_n[,-1])
  Mu_qr <- qr(M_u, LAPACK = FALSE, tol = 1e-15)
  Q_u <- qr.Q(Mu_qr)
  sgn <- sign(diag(qr.R(Mu_qr)))
  Q_u.new <- Q_u %*% diag(sgn)
  
  
  Me_qr <- qr(M_e, LAPACK = FALSE, tol = 1e-15)
  Q_e <- qr.Q(Me_qr)
  sgn <- sign(diag(qr.R(Me_qr)))
  Q_e.new <- Q_e %*% diag(sgn)
  
  Q_e.new %*% I_n[,1]
  Q_u.new %*% I_n[,1]
  
  R_u <- Q_e.new%*% t(Q_u.new)
  
  for(i in 1:nrow(X)){
    vertex <- X[i,]
    aux <- t(R_u %*% t(sweep(X, 2, vertex) ))
    aux1 <- (aux >= c(0,0))
    prop <- sum(aux1[,1]*aux1[,2]) / nrow(X)
    if(prop-alfa>=h){
      puntos_cuantil[i] <- TRUE
    }else{
      puntos_cuantil[i] <- FALSE
    }
  }
  return(puntos_cuantil = puntos_cuantil)
}


obtener_cuantiles_mltools <- function(alfa,h, X){
  
  X <- data.table::data.table(X)
  aux <- apply(X, 1, function(x) mltools::empirical_cdf(X, ubounds = data.table::data.table(x[1],x[2]))$CDF)
  aux1 <- aux-alfa>=h
  return(puntos_cuantil = aux1)
}

obtener_cuantiles_baseR <- function(alfa, h, X){
  
  data <- as.matrix(as.data.frame(X))
  dt <- data.table::data.table(data)
  
  base_R <- pmap_dbl(dt, .f = function(V1, V2) {
    mean(data[, "V1"] <= V1 & data[, "V2"] <= V2)
  })
  aux <- base_R-alfa>=h
  return(puntos_cuantil = aux)
}



dtrgenbeta <- function(x,shape1, shape2, shape3, scale,a=0,b=100){
  f_x <- dgenbeta(x, 
                  shape1 = shape1,
                  shape2 = shape2,
                  shape3 = shape3,
                  scale = scale)
  f_x / (pgenbeta(b, 
                  shape1 = shape1,
                  shape2 = shape2,
                  shape3 = shape3,
                  scale = scale) - pgenbeta(a, 
                                            shape1 = shape1,
                                            shape2 = shape2,
                                            shape3 = shape3,
                                            scale = scale))
}


ptrgenbeta <- function(q,shape1, shape2, shape3, scale,a=0,b=100){
  aux <- (pgenbeta(q,
                   shape1 = shape1,
                   shape2 = shape2,
                   shape3 = shape3,
                   scale = scale) -
            pgenbeta(a, 
                     shape1 = shape1,
                     shape2 = shape2,
                     shape3 = shape3,
                     scale = scale))/(pgenbeta(b, 
                                               shape1 = shape1,
                                               shape2 = shape2,
                                               shape3 = shape3,
                                               scale = scale) - pgenbeta(a, 
                                                                         shape1 = shape1,
                                                                         shape2 = shape2,
                                                                         shape3 = shape3,
                                                                         scale = scale))
  if(q>b) aux <- 1
  return(aux)
}

ptrgenbeta <- Vectorize(ptrgenbeta,"q")


rtrgenbeta <- function(n, shape1, shape2, shape3, scale,a=0,b=100) {
  
  F.a <- pgenbeta(a, shape1 = shape1,
                  shape2 = shape2,
                  shape3 = shape3,
                  scale = scale)
  F.b <- pgamma(b, shape1 = shape1,
                shape2 = shape2,
                shape3 = shape3,
                scale = scale)
  
  u <- runif(n, min = F.a, max = F.b)
  
  qgenbeta(u, shape1 = shape1,
           shape2 = shape2,
           shape3 = shape3,
           scale = scale)
}

qtrgenbeta <- function(p,shape1, shape2, shape3, scale,a=0,b=100){
  qgenbeta(pgenbeta(a,
                    shape1 = shape1,
                    shape2 = shape2,
                    shape3 = shape3,
                    scale = scale) + p * (pgenbeta(b,
                                                   shape1 = shape1,
                                                   shape2 = shape2,
                                                   shape3 = shape3,
                                                   scale = scale) -pgenbeta(a,
                                                                            shape1 = shape1,
                                                                            shape2 = shape2,
                                                                            shape3 = shape3,
                                                                            scale = scale)),
           shape1 = shape1,
           shape2 = shape2,
           shape3 = shape3,
           scale = scale)
}
