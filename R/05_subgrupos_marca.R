# ============================================================
# 05_subgrupos_marca.R
# (a) Análisis por marca comercial usando TODAS las entradas
#     (la deduplicación por ranking de mercado no aplica cuando el
#      análisis se estratifica por marca: cada estudio aporta una
#      sola entrada por marca y estándar de referencia)
# (b) Prueba de Deeks forzada con k < 10 (exploratoria)
# (c) Metarregresión multivariable exploratoria sobre logDOR
# (d) Tabla 1: características de los estudios de la revisión
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(mada); library(metafor); library(writexl); library(tidyr)
})
OUT <- "salidas"
DAT <- readRDS(file.path(OUT, "datasets_dta_campy.rds"))
qd_cols <- c("rob_PS","rob_IT","rob_RS","rob_FT","ac_PS","ac_IT","ac_RS")

fit_biv <- function(d) try(reitsma(d, TP="TP", FN="FN", FP="FP", TN="TN", method="reml"), silent=TRUE)
est_biv <- function(f) {
  if (inherits(f,"try-error")) return(rep(NA_real_, 7))
  s <- summary(f); co <- s$coefficients
  c(co["sensitivity","Estimate"], co["sensitivity","95%ci.lb"], co["sensitivity","95%ci.ub"],
    1-co["false pos. rate","Estimate"], 1-co["false pos. rate","95%ci.ub"], 1-co["false pos. rate","95%ci.lb"],
    if (!is.null(s$AUC$AUC)) s$AUC$AUC else NA_real_)
}
calc_logdor <- function(d, cc = 0.5) {
  tp<-d$TP+cc; fn<-d$FN+cc; fp<-d$FP+cc; tn<-d$TN+cc
  d$.logDOR <- log((tp*tn)/(fp*fn)); d$.var <- 1/tp+1/fn+1/fp+1/tn; d
}

# ---------- (a) Subgrupos por marca ----------
# Base: todas las entradas de las 4 comparaciones centrales (sin deduplicar).
todas <- bind_rows(DAT$all)
todas$Referencia <- ifelse(todas$Reference_2 == "Culture", "Cultivo", "PCR")
todas$Plataforma <- ifelse(todas$Index_2 == "microplate assay", "EIA microplaca", "ICA/LFA")

# Entradas de la revisión que no entran en las 4 comparaciones (marcas de
# plataformas no comparables), sólo para descripción por marca.
dc <- DAT$datos_campy
fuera <- dc %>%
  filter(!Record_ID %in% DAT$datos_comparaciones$Record_ID,
         is.na(eligibility) | eligibility != 0,
         Reference_2 %in% c("Culture","PCR")) %>%
  mutate(Referencia = ifelse(Reference_2 == "Culture", "Cultivo", "PCR"),
         Plataforma = case_when(Index_2 == "rapid membrane EIA" ~ "EIA de membrana rápida",
                                Index_2 == "chemiluminescent assay" ~ "Quimioluminiscente",
                                TRUE ~ Index_2))

marca_resumen <- function(d, en_analisis) {
  d %>% group_by(Marca = Index_1, Plataforma, Referencia) %>%
    group_modify(function(g, key) {
      f <- if (nrow(g) >= 4) fit_biv(g) else NULL
      e <- if (is.null(f)) rep(NA_real_, 7) else est_biv(f)
      data.frame(k = nrow(g),
                 estudios = paste(sort(unique(g$Study_ID)), collapse = "; "),
                 N = sum(g$TP+g$FN+g$FP+g$TN),
                 Se_bruta = sum(g$TP)/sum(g$TP+g$FN),
                 Sp_bruta = sum(g$TN)/sum(g$TN+g$FP),
                 Se_min = min(g$TP/(g$TP+g$FN)), Se_max = max(g$TP/(g$TP+g$FN)),
                 Sp_min = min(g$TN/(g$TN+g$FP)), Sp_max = max(g$TN/(g$TN+g$FP)),
                 Se_biv = e[1], Se_lb = e[2], Se_ub = e[3],
                 Sp_biv = e[4], Sp_lb = e[5], Sp_ub = e[6], AUC = e[7],
                 modelo = if (nrow(g) >= 4) "bivariado de Reitsma" else "sólo descriptivo (k<4)")
    }) %>% ungroup() %>% mutate(en_analisis_conjunto = en_analisis)
}

subgrupos_marca <- bind_rows(
  marca_resumen(todas, "sí"),
  marca_resumen(fuera, "no (plataforma no comparable)")
) %>% arrange(desc(en_analisis_conjunto == "sí"), Referencia, Plataforma, desc(k))

# Metarregresión por marca dentro de cada comparación (marcas con k >= 2)
marca_metareg <- bind_rows(lapply(names(DAT$all), function(nm) {
  d <- DAT$all[[nm]]
  tb <- table(d$Index_1); keep <- names(tb)[tb >= 2]
  d2 <- d[d$Index_1 %in% keep, ]
  if (length(keep) < 2 || nrow(d2) < 5) {
    return(data.frame(dataset = nm, k = nrow(d2), marcas = paste(keep, collapse="; "),
                      chi2 = NA_real_, gl = NA_real_, p = NA_real_,
                      estado = "insuficiente (menos de 2 marcas con k>=2 o k<5)"))
  }
  d2$.mod <- droplevels(factor(d2$Index_1))
  f0 <- try(reitsma(d2, TP="TP",FN="FN",FP="FP",TN="TN", method="ml"), silent=TRUE)
  f1 <- try(reitsma(d2, formula = cbind(tsens,tfpr) ~ .mod, TP="TP",FN="FN",FP="FP",TN="TN",
                    method="ml"), silent=TRUE)
  if (inherits(f0,"try-error") || inherits(f1,"try-error"))
    return(data.frame(dataset = nm, k = nrow(d2), marcas = paste(keep, collapse="; "),
                      chi2 = NA_real_, gl = NA_real_, p = NA_real_, estado = "no convergió"))
  a <- anova(f1, f0)
  data.frame(dataset = nm, k = nrow(d2), marcas = paste(keep, collapse="; "),
             chi2 = unname(a$statistic["Chi-squared"]), gl = unname(a$statistic["Df"]),
             p = unname(a$statistic["Pr(>Chi-squared)"]), estado = "estimado")
}))

# ---------- (b) Deeks forzado con k < 10 ----------
deeks_calc <- function(d) {
  d <- calc_logdor(d)
  n1 <- d$TP+d$FN; n0 <- d$FP+d$TN; ess <- 4*n1*n0/(n1+n0)
  fit <- lm(d$.logDOR ~ I(1/sqrt(ess)), weights = ess)
  s <- summary(fit)$coefficients
  data.frame(k = nrow(d), intercepto = s[1,1], pendiente = s[2,1], EE = s[2,2],
             t = s[2,3], p = s[2,4], gl_residual = fit$df.residual)
}
deeks_forzado <- bind_rows(lapply(names(DAT$global), function(cmp) {
  cbind(data.frame(Comparacion = cmp, conjunto = "global deduplicado"),
        deeks_calc(DAT$global[[cmp]]))
})) %>%
  bind_rows(bind_rows(lapply(names(DAT$all), function(nm) {
    cbind(data.frame(Comparacion = nm, conjunto = "todas las entradas"),
          deeks_calc(DAT$all[[nm]]))
  }))) %>%
  mutate(interpretacion = ifelse(k >= 10, "prueba estándar (k>=10)",
                                 "exploratoria: potencia insuficiente, no interpretar el p"))

# ---------- (c) Metarregresión multivariable exploratoria (logDOR) ----------
multivar <- bind_rows(lapply(c("LFA_Culture","EIA_Culture"), function(cmp) {
  d <- calc_logdor(DAT$global[[cmp]])
  d$rob_RS_f <- factor(ifelse(d$rob_RS == 1, "bajo", "incierto"), levels = c("bajo","incierto"))
  d$prosp_f  <- factor(ifelse(d$prospective == 1, "prospectivo", "no prospectivo"),
                       levels = c("prospectivo","no prospectivo"))
  m1 <- rma(yi = .logDOR, vi = .var, mods = ~ rob_RS_f, data = d, method = "REML")
  m2 <- rma(yi = .logDOR, vi = .var, mods = ~ prosp_f, data = d, method = "REML")
  m3 <- rma(yi = .logDOR, vi = .var, mods = ~ rob_RS_f + prosp_f, data = d, method = "REML")
  filas <- function(m, etq) {
    cf <- coef(summary(m))
    data.frame(Comparacion = cmp, modelo = etq, k = m$k,
               parametros = m$p + 1, termino = rownames(cf),
               beta = cf$estimate, EE = cf$se, p = cf$pval,
               QM = m$QM, QM_p = m$QMp, tau2 = m$tau2, I2 = m$I2,
               AICc = fitstats(m)["AICc",1], row.names = NULL)
  }
  bind_rows(filas(m1, "univariable: rob_RS"), filas(m2, "univariable: prospectivo"),
            filas(m3, "multivariable: rob_RS + prospectivo"))
}))

# ---------- (d) Tabla 1: características de los estudios ----------
en_comp <- DAT$datos_comparaciones %>%
  mutate(comp = paste0(ifelse(Index_2 == "microplate assay", "EIA", "LFA"), "-",
                       ifelse(Reference_2 == "Culture", "Cultivo", "PCR"))) %>%
  group_by(Study_ID) %>% summarise(comparaciones = paste(sort(unique(comp)), collapse = ", "),
                                   .groups = "drop")

tabla1 <- dc %>%
  group_by(Study_ID) %>%
  summarise(
    Anio = suppressWarnings(min(year, na.rm = TRUE)),
    Pais = paste(unique(na.omit(country)), collapse = "/"),
    n_muestras = suppressWarnings(max(n, na.rm = TRUE)),
    Prevalencia_pct = round(100*max(sample_prev_calc, na.rm = TRUE), 1),
    Pruebas_indice = paste(sort(unique(Index_1)), collapse = "; "),
    Plataformas = paste(sort(unique(Index_2)), collapse = "; "),
    Estandares = paste(sort(unique(Reference_2)), collapse = "; "),
    Diseno_prospectivo = case_when(all(prospective == 1, na.rm = TRUE) ~ "Sí",
                                   all(prospective == 0, na.rm = TRUE) ~ "No",
                                   TRUE ~ "Mixto/NI"),
    Muestras_preseleccionadas = case_when(all(preselected == 1, na.rm = TRUE) ~ "Sí",
                                          all(preselected == 0, na.rm = TRUE) ~ "No",
                                          TRUE ~ "Mixto/NI"),
    Lectura = paste(unique(na.omit(result_reading)), collapse = "/"),
    QUADAS_dominios_inciertos = max(rowSums(across(all_of(qd_cols), ~ .x == 3), na.rm = TRUE)),
    QUADAS_dominios_alto = max(rowSums(across(all_of(qd_cols), ~ .x == 2), na.rm = TRUE)),
    Elegibilidad = paste(sort(unique(eligibility)), collapse = "/"),
    Entradas_2x2 = n(),
    .groups = "drop") %>%
  left_join(en_comp, by = "Study_ID") %>%
  mutate(comparaciones = ifelse(is.na(comparaciones), "sólo descripción narrativa", comparaciones)) %>%
  arrange(Study_ID)

tabla1_entradas <- dc %>%
  transmute(Estudio = Study_ID, Anio = year, Pais = country,
            Prueba_indice = Index_1, Plataforma = Index_2, Estandar = Reference_2,
            n, TP, FN, FP, TN,
            Prevalencia_pct = round(100*sample_prev_calc, 1),
            Sensibilidad = round(sensitivity, 3), Especificidad = round(specificity, 3),
            Elegibilidad = eligibility,
            En_analisis = ifelse(Record_ID %in% DAT$datos_comparaciones$Record_ID, "sí", "no")) %>%
  arrange(Estudio, Plataforma, Estandar, Prueba_indice)

write_xlsx(list(Tabla1_estudios = tabla1, Tabla1_entradas_2x2 = tabla1_entradas,
                subgrupos_marca = subgrupos_marca, marca_metarregresion = marca_metareg,
                deeks_forzado = deeks_forzado, metarreg_multivariable = multivar),
           file.path(OUT, "resultados_subgrupos_marca.xlsx"))
write.csv(subgrupos_marca, file.path(OUT, "subgrupos_marca.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tabla1, file.path(OUT, "tabla1_caracteristicas.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tabla1_entradas, file.path(OUT, "tabla1_entradas_2x2.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("== subgrupos por marca ==\n")
print(as.data.frame(subgrupos_marca[, c("Marca","Plataforma","Referencia","k","N","Se_bruta","Sp_bruta","Se_biv","Sp_biv","modelo")]))
cat("\n== metarregresión por marca ==\n"); print(marca_metareg)
cat("\n== Deeks forzado ==\n"); print(deeks_forzado[, c("Comparacion","conjunto","k","pendiente","p","interpretacion")])
cat("\n== multivariable ==\n"); print(multivar[, c("Comparacion","modelo","parametros","termino","beta","p","QM_p","I2","AICc")])
