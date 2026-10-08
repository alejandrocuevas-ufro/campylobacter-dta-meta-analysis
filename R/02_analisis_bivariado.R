# ============================================================
# 02_analisis_bivariado.R
# Metaanálisis DTA de inmunoensayos para Campylobacter spp.
# Modelo bivariado de Reitsma (mada), heterogeneidad y comparación EIA vs LFA
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(mada); library(metafor); library(writexl)
})

OUT <- "salidas"; dir.create(OUT, showWarnings = FALSE)
DAT <- readRDS(file.path(OUT, "datasets_dta_campy.rds"))
comparaciones <- c("LFA_Culture","LFA_PCR","EIA_Culture","EIA_PCR")
etiquetas <- c(LFA_Culture = "ICA/LFA vs Cultivo", LFA_PCR = "ICA/LFA vs PCR",
               EIA_Culture = "EIA microplaca vs Cultivo", EIA_PCR = "EIA microplaca vs PCR")

ci_exact <- function(x, n) {
  if (n == 0) return(c(NA, NA))
  as.numeric(binom.test(x, n)$conf.int)
}

# ---- 1. Tabla de estudios individuales ----
tabla_estudios <- bind_rows(lapply(comparaciones, function(cmp) {
  d <- DAT$global[[cmp]]
  ci_se <- t(mapply(ci_exact, d$TP, d$TP + d$FN))
  ci_sp <- t(mapply(ci_exact, d$TN, d$TN + d$FP))
  data.frame(Comparacion = cmp, Etiqueta = etiquetas[[cmp]],
             Estudio = d$Study_ID, Prueba_indice = d$Index_1, Anio = d$year,
             Pais = if ("country" %in% names(d)) d$country else NA,
             N = d$n, TP = d$TP, FN = d$FN, FP = d$FP, TN = d$TN,
             Prevalencia = round(d$sample_prev_calc, 4),
             Sensibilidad = round(d$TP/(d$TP+d$FN), 4),
             Se_IC_inf = round(ci_se[,1],4), Se_IC_sup = round(ci_se[,2],4),
             Especificidad = round(d$TN/(d$TN+d$FP), 4),
             Sp_IC_inf = round(ci_sp[,1],4), Sp_IC_sup = round(ci_sp[,2],4),
             row.names = NULL)
}))

# ---- 2. Modelo bivariado por comparación ----
fits <- list(); resumen <- list(); heterog <- list()

for (cmp in comparaciones) {
  d <- DAT$global[[cmp]]
  f <- reitsma(d, TP="TP", FN="FN", FP="FP", TN="TN", method = "reml")
  s <- summary(f)
  fits[[cmp]] <- f
  co <- s$coefficients
  set.seed(20260918)
  sp <- summary(SummaryPts(f, n.iter = 1e6))

  # heterogeneidad univariada (logit) con metafor
  es_se <- escalc(measure="PLO", xi=d$TP, ni=d$TP+d$FN)
  es_sp <- escalc(measure="PLO", xi=d$TN, ni=d$TN+d$FP)
  m_se <- rma(yi, vi, data = es_se, method = "REML")
  m_sp <- rma(yi, vi, data = es_sp, method = "REML")

  resumen[[cmp]] <- data.frame(
    Comparacion = cmp, Etiqueta = etiquetas[[cmp]],
    k_estudios = nrow(d), N_total = sum(d$TP+d$FN+d$FP+d$TN),
    Sensibilidad = co["sensitivity","Estimate"],
    Se_IC_inf = co["sensitivity","95%ci.lb"], Se_IC_sup = co["sensitivity","95%ci.ub"],
    Especificidad = 1 - co["false pos. rate","Estimate"],
    Sp_IC_inf = 1 - co["false pos. rate","95%ci.ub"], Sp_IC_sup = 1 - co["false pos. rate","95%ci.lb"],
    LR_pos = sp["posLR","Median"], LRpos_IC_inf = sp["posLR","2.5%"], LRpos_IC_sup = sp["posLR","97.5%"],
    LR_neg = sp["negLR","Median"], LRneg_IC_inf = sp["negLR","2.5%"], LRneg_IC_sup = sp["negLR","97.5%"],
    DOR = sp["DOR","Median"], DOR_IC_inf = sp["DOR","2.5%"], DOR_IC_sup = sp["DOR","97.5%"],
    AUC_SROC = s$AUC$AUC, row.names = NULL)

  heterog[[cmp]] <- data.frame(
    Comparacion = cmp, Etiqueta = etiquetas[[cmp]], k_estudios = nrow(d),
    tau2_logitSe_biv = f$Psi[1,1], tau2_logitFPR_biv = f$Psi[2,2],
    cor_aleatoria_biv = s$corRandom[1,2],
    I2_logitSe = m_se$I2, tau2_logitSe_univ = m_se$tau2,
    Q_Se = m_se$QE, p_Q_Se = m_se$QEp,
    I2_logitSp = m_sp$I2, tau2_logitSp_univ = m_sp$tau2,
    Q_Sp = m_sp$QE, p_Q_Sp = m_sp$QEp, row.names = NULL)
}
resumen <- bind_rows(resumen); heterog <- bind_rows(heterog)

# ---- 3. Comparación EIA vs LFA (metarregresión bivariada por estándar de referencia) ----
comparar_plataforma <- function(ref) {
  d <- bind_rows(DAT$global[[paste0("LFA_", ref)]], DAT$global[[paste0("EIA_", ref)]])
  d$platform <- factor(ifelse(d$Index_2 == "microplate assay", "EIA", "LFA"), levels = c("LFA","EIA"))
  f0 <- reitsma(d, TP="TP", FN="FN", FP="FP", TN="TN", method = "ml")
  f1 <- reitsma(d, formula = cbind(tsens, tfpr) ~ platform,
                TP="TP", FN="FN", FP="FP", TN="TN", method = "ml")
  lrt <- anova(f1, f0)
  s1 <- summary(f1)$coefficients
  list(ref = ref, d = d, f0 = f0, f1 = f1, lrt = lrt, coef = s1,
       n_estudios_compartidos = sum(duplicated(d$Study_ID)))
}
cmp_culture <- comparar_plataforma("Culture")
cmp_pcr     <- comparar_plataforma("PCR")

extraer_lrt <- function(x) {
  st <- x$lrt$statistic
  data.frame(Referencia = x$ref,
             k_total = nrow(x$d),
             estudios_con_ambas_plataformas = x$n_estudios_compartidos,
             chi2_LRT = unname(st["Chi-squared"]),
             gl = unname(st["Df"]),
             p_valor = unname(st["Pr(>Chi-squared)"]),
             row.names = NULL)
}
tabla_comparacion <- bind_rows(lapply(list(cmp_culture, cmp_pcr), extraer_lrt))

coef_comparacion <- bind_rows(lapply(list(cmp_culture, cmp_pcr), function(x) {
  cf <- as.data.frame(x$coef); cf$Parametro <- rownames(cf); cf$Referencia <- x$ref
  cf[, c("Referencia","Parametro","Estimate","Std. Error","z","Pr(>|z|)","95%ci.lb","95%ci.ub")]
}))

# ---- 4. Sesgo de publicación: prueba de Deeks (sólo k >= 10) ----
deeks <- function(d) {
  tp <- d$TP + 0.5; fn <- d$FN + 0.5; fp <- d$FP + 0.5; tn <- d$TN + 0.5
  lnDOR <- log((tp*tn)/(fp*fn))
  n1 <- d$TP + d$FN; n0 <- d$FP + d$TN
  ess <- 4*n1*n0/(n1+n0)
  fit <- lm(lnDOR ~ I(1/sqrt(ess)), weights = ess)
  s <- summary(fit)$coefficients
  data.frame(pendiente = s[2,1], EE = s[2,2], t = s[2,3], p_valor = s[2,4], row.names = NULL)
}
tabla_deeks <- bind_rows(lapply(comparaciones, function(cmp) {
  d <- DAT$global[[cmp]]
  if (nrow(d) < 10) return(data.frame(Comparacion = cmp, Etiqueta = etiquetas[[cmp]],
                                      k = nrow(d), pendiente = NA, EE = NA, t = NA,
                                      p_valor = NA, nota = "k<10: no evaluado"))
  cbind(data.frame(Comparacion = cmp, Etiqueta = etiquetas[[cmp]], k = nrow(d)),
        deeks(d), nota = "")
}))

# ---- 5. Exportar ----
write_xlsx(list(estudios = tabla_estudios, resumen_bivariado = resumen,
                heterogeneidad = heterog, comparacion_EIA_LFA = tabla_comparacion,
                coef_metarregresion = coef_comparacion, deeks = tabla_deeks),
           file.path(OUT, "resultados_dta_campylobacter.xlsx"))
write.csv(resumen, file.path(OUT, "resumen_bivariado.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(heterog, file.path(OUT, "heterogeneidad.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tabla_estudios, file.path(OUT, "estudios_incluidos.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(fits = fits, cmp_culture = cmp_culture, cmp_pcr = cmp_pcr,
             resumen = resumen, heterog = heterog, tabla_estudios = tabla_estudios,
             tabla_comparacion = tabla_comparacion, coef_comparacion = coef_comparacion,
             tabla_deeks = tabla_deeks),
        file.path(OUT, "modelos_dta_campy.rds"))

print(resumen); cat("\n"); print(heterog); cat("\n"); print(tabla_comparacion); cat("\n"); print(tabla_deeks)
