# ============================================================
# 04_bayesiano_jags.R
# Modelo bivariado bayesiano binomial-normal (Chu & Cole) en JAGS,
# con la parametrización de priors documentada para MetaBayesDTA
# (archivo "parametros metabayes.xlsx"):
#   - Intervalo prior 95% para la sensibilidad agrupada: 0,40 - 0,99
#   - Intervalo prior 95% para la especificidad agrupada: 0,80 - 0,995
#     -> se traducen a N(media, sd) en escala logit
#   - Prior de la DE entre estudios: semi-normal(0; 0,75) para logit(Se)
#                                     semi-normal(0; 0,50) para logit(Sp)
#   - 4 cadenas, 2000 de calentamiento, 6000 iteraciones totales, semilla 123
# NOTA: no es el mismo software (MetaBayesDTA usa Stan); es el mismo modelo
# y la misma especificación de priors implementada en JAGS.
# ============================================================

suppressPackageStartupMessages({ library(rjags); library(coda); library(dplyr); library(writexl) })
OUT <- "salidas"
DAT <- readRDS(file.path(OUT, "datasets_dta_campy.rds"))
comparaciones <- c("LFA_Culture","EIA_Culture","LFA_PCR","EIA_PCR")
etiquetas <- c(LFA_Culture="ICA/LFA vs Cultivo", EIA_Culture="EIA microplaca vs Cultivo",
               LFA_PCR="ICA/LFA vs PCR", EIA_PCR="EIA microplaca vs PCR")

logit <- function(p) log(p/(1-p))
prior_normal <- function(lo, hi) {
  m <- (logit(lo) + logit(hi))/2
  s <- (logit(hi) - logit(lo))/(2*qnorm(0.975))
  c(mean = m, sd = s)
}
PR_SE <- prior_normal(0.40, 0.99)
PR_SP <- prior_normal(0.80, 0.995)
SD_SE <- 0.75; SD_SP <- 0.50

modelo <- "
model {
  for (i in 1:K) {
    TP[i] ~ dbin(se[i], n1[i])
    TN[i] ~ dbin(sp[i], n0[i])
    logit(se[i]) <- theta[i,1]
    logit(sp[i]) <- theta[i,2]
    theta[i,1:2] ~ dmnorm.vcov(mu[1:2], Sigma[1:2,1:2])
  }
  mu[1] ~ dnorm(m_se, pow(s_se, -2))
  mu[2] ~ dnorm(m_sp, pow(s_sp, -2))
  tau1 ~ dnorm(0, pow(sd_se, -2)) T(0,)
  tau2 ~ dnorm(0, pow(sd_sp, -2)) T(0,)
  rho ~ dunif(-1, 1)
  Sigma[1,1] <- pow(tau1, 2)
  Sigma[2,2] <- pow(tau2, 2)
  Sigma[1,2] <- rho*tau1*tau2
  Sigma[2,1] <- Sigma[1,2]

  Se <- ilogit(mu[1])
  Sp <- ilogit(mu[2])
  LRpos <- Se/(1-Sp)
  LRneg <- (1-Se)/Sp
  DOR <- (Se/(1-Se))/((1-Sp)/Sp)
  # estudio nuevo (predicción)
  theta_new[1:2] ~ dmnorm.vcov(mu[1:2], Sigma[1:2,1:2])
  Se_new <- ilogit(theta_new[1])
  Sp_new <- ilogit(theta_new[2])
}
"
writeLines(modelo, file.path(OUT, "modelo_bivariado.jags"))

correr_jags <- function(d, semilla = 123) {
  datos <- list(K = nrow(d), TP = d$TP, TN = d$TN,
                n1 = d$TP + d$FN, n0 = d$FP + d$TN,
                m_se = PR_SE[["mean"]], s_se = PR_SE[["sd"]],
                m_sp = PR_SP[["mean"]], s_sp = PR_SP[["sd"]],
                sd_se = SD_SE, sd_sp = SD_SP)
  set.seed(semilla)
  inits <- lapply(1:4, function(ch) list(mu = c(PR_SE[["mean"]], PR_SP[["mean"]]),
                                         tau1 = 0.3, tau2 = 0.3, rho = 0,
                                         .RNG.name = "base::Mersenne-Twister",
                                         .RNG.seed = semilla + ch))
  m <- jags.model(file.path(OUT, "modelo_bivariado.jags"), data = datos,
                  inits = inits, n.chains = 4, quiet = TRUE)
  update(m, 2000)
  # 2000 de calentamiento + 4000 por cadena según el archivo de priors; se
  # extiende a 20000 por cadena porque con k=5 el ESS de tau2 y Sp quedaba < 250.
  coda.samples(m, c("Se","Sp","LRpos","LRneg","DOR","tau1","tau2","rho","Se_new","Sp_new"),
               n.iter = 20000)
}

resumir <- function(s, cmp, k) {
  M <- as.matrix(s)
  q <- t(apply(M, 2, function(z) c(media = mean(z), mediana = median(z),
                                   IC_inf = quantile(z, 0.025), IC_sup = quantile(z, 0.975))))
  rh <- try(gelman.diag(s, multivariate = FALSE)$psrf[,1], silent = TRUE)
  ess <- effectiveSize(s)
  data.frame(Comparacion = cmp, Etiqueta = etiquetas[[cmp]], k = k,
             Parametro = rownames(q), Media = q[,1], Mediana = q[,2],
             IC95_inf = q[,3], IC95_sup = q[,4],
             Rhat = if (inherits(rh, "try-error")) NA_real_ else rh[rownames(q)],
             ESS = ess[rownames(q)], row.names = NULL)
}

bayes_tab <- bind_rows(lapply(comparaciones, function(cmp) {
  d <- DAT$global[[cmp]]
  s <- correr_jags(d)
  resumir(s, cmp, nrow(d))
}))

priors_usados <- data.frame(
  Parametro = c("logit(Se) media prior","logit(Se) DE prior","logit(Sp) media prior",
                "logit(Sp) DE prior","DE entre estudios logit(Se)","DE entre estudios logit(Sp)",
                "Correlación","Cadenas","Calentamiento","Iteraciones post-calentamiento","Semilla"),
  Valor = c(round(PR_SE[["mean"]],4), round(PR_SE[["sd"]],4), round(PR_SP[["mean"]],4),
            round(PR_SP[["sd"]],4), "semi-normal(0; 0,75)", "semi-normal(0; 0,50)",
            "uniforme(-1, 1)", 4, 2000, 4000, 123),
  Origen = c(rep("derivado del intervalo 0,40-0,99 del archivo de priors", 2),
             rep("derivado del intervalo 0,80-0,995 del archivo de priors", 2),
             rep("archivo de priors", 2), "no especificado en el archivo; prior no informativo",
             rep("archivo de priors", 4)))

write_xlsx(list(bayesiano = bayes_tab, priors = priors_usados),
           file.path(OUT, "resultados_bayesianos.xlsx"))
write.csv(bayes_tab, file.path(OUT, "bayes_tab.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(bayes_tab %>% filter(Parametro %in% c("Se","Sp","DOR","tau1","tau2","rho")) %>%
        mutate(across(where(is.numeric), ~round(.x, 4))))
