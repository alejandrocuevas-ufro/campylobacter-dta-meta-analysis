# ============================================================
# 06_bayes_priors.R
# Análisis de sensibilidad a los priors del modelo bivariado bayesiano.
# Tres especificaciones:
#   A. Original del proyecto (revisión de literatura con consensus.app):
#      Se 0,40-0,99 ; Sp 0,80-0,995 ; DE semi-normal(0,75) y (0,50)
#   B. Ampliada (recomendada): mismo centro, intervalos creíbles más anchos
#      Se 0,25-0,995 ; Sp 0,60-0,999 ; DE semi-normal(1,0) en ambos
#   C. Vaga (referencia): mu ~ N(0, 3) en logit ; DE semi-normal(1,5)
# En las tres, correlación ~ uniforme(-1, 1).
# ============================================================

suppressPackageStartupMessages({ library(rjags); library(coda); library(dplyr); library(writexl) })
OUT <- "salidas"
DAT <- readRDS(file.path(OUT, "datasets_dta_campy.rds"))
comparaciones <- c("LFA_Culture","EIA_Culture","LFA_PCR","EIA_PCR")
etiquetas <- c(LFA_Culture="ICA/LFA vs Cultivo", EIA_Culture="EIA microplaca vs Cultivo",
               LFA_PCR="ICA/LFA vs PCR", EIA_PCR="EIA microplaca vs PCR")
logit <- function(p) log(p/(1-p))
de_intervalo <- function(lo, hi) c(mean = (logit(lo)+logit(hi))/2,
                                   sd = (logit(hi)-logit(lo))/(2*qnorm(0.975)))

PRIORS <- list(
  A_original = list(se = de_intervalo(0.40, 0.99),  sp = de_intervalo(0.80, 0.995),
                    sd_se = 0.75, sd_sp = 0.50,
                    desc = "Original del proyecto: Se 0,40-0,99; Sp 0,80-0,995; DE 0,75/0,50"),
  B_ampliado = list(se = de_intervalo(0.25, 0.995), sp = de_intervalo(0.60, 0.999),
                    sd_se = 1.00, sd_sp = 1.00,
                    desc = "Ampliado: Se 0,25-0,995; Sp 0,60-0,999; DE 1,0/1,0"),
  C_vago     = list(se = c(mean = 0, sd = 3), sp = c(mean = 0, sd = 3),
                    sd_se = 1.50, sd_sp = 1.50,
                    desc = "Vago: mu ~ N(0, 3) en logit; DE 1,5/1,5")
)

correr <- function(d, pr, semilla = 123) {
  datos <- list(K = nrow(d), TP = d$TP, TN = d$TN,
                n1 = d$TP + d$FN, n0 = d$FP + d$TN,
                m_se = pr$se[["mean"]], s_se = pr$se[["sd"]],
                m_sp = pr$sp[["mean"]], s_sp = pr$sp[["sd"]],
                sd_se = pr$sd_se, sd_sp = pr$sd_sp)
  inits <- lapply(1:4, function(ch) list(mu = c(pr$se[["mean"]], pr$sp[["mean"]]),
                                         tau1 = 0.3, tau2 = 0.3, rho = 0,
                                         .RNG.name = "base::Mersenne-Twister",
                                         .RNG.seed = semilla + ch))
  m <- jags.model(file.path(OUT, "modelo_bivariado.jags"), data = datos, inits = inits,
                  n.chains = 4, quiet = TRUE)
  update(m, 2000)
  coda.samples(m, c("Se","Sp","DOR","tau1","tau2","rho","Se_new","Sp_new"), n.iter = 20000)
}

res <- bind_rows(lapply(names(PRIORS), function(nm) {
  pr <- PRIORS[[nm]]
  bind_rows(lapply(comparaciones, function(cmp) {
    d <- DAT$global[[cmp]]
    s <- correr(d, pr)
    q <- t(apply(as.matrix(s), 2, function(z) c(mediana = median(z),
                                                lb = quantile(z, 0.025), ub = quantile(z, 0.975))))
    rh <- gelman.diag(s, multivariate = FALSE)$psrf[,1]
    data.frame(Prior = nm, Descripcion = pr$desc, Comparacion = cmp,
               Etiqueta = etiquetas[[cmp]], k = nrow(d), Parametro = rownames(q),
               Mediana = q[,1], ICr_inf = q[,2], ICr_sup = q[,3],
               Rhat = rh[rownames(q)], ESS = effectiveSize(s)[rownames(q)], row.names = NULL)
  }))
}))

priors_tab <- bind_rows(lapply(names(PRIORS), function(nm) {
  pr <- PRIORS[[nm]]
  data.frame(Prior = nm, Descripcion = pr$desc,
             mu_logitSe_media = pr$se[["mean"]], mu_logitSe_DE = pr$se[["sd"]],
             mu_logitSp_media = pr$sp[["mean"]], mu_logitSp_DE = pr$sp[["sd"]],
             DE_entre_estudios_Se = pr$sd_se, DE_entre_estudios_Sp = pr$sd_sp,
             Correlacion = "uniforme(-1, 1)",
             Se_prior_IC95 = paste0(round(plogis(pr$se[["mean"]] - 1.96*pr$se[["sd"]]), 3), " - ",
                                    round(plogis(pr$se[["mean"]] + 1.96*pr$se[["sd"]]), 3)),
             Sp_prior_IC95 = paste0(round(plogis(pr$sp[["mean"]] - 1.96*pr$sp[["sd"]]), 3), " - ",
                                    round(plogis(pr$sp[["mean"]] + 1.96*pr$sp[["sd"]]), 3)),
             row.names = NULL)
}))

write_xlsx(list(sensibilidad_priors = res, especificacion_priors = priors_tab),
           file.path(OUT, "resultados_bayes_priors.xlsx"))
write.csv(res, file.path(OUT, "bayes_priors.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(priors_tab[, c("Prior","Se_prior_IC95","Sp_prior_IC95","DE_entre_estudios_Se")])
print(res %>% filter(Parametro %in% c("Se","Sp")) %>%
        select(Prior, Comparacion, Parametro, Mediana, ICr_inf, ICr_sup) %>%
        mutate(across(where(is.numeric), ~round(.x, 3))) %>% as.data.frame())
cat("Rhat máximo:", max(res$Rhat, na.rm=TRUE), " ESS mínimo:", min(res$ESS, na.rm=TRUE), "\n")
