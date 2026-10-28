# ============================================================
# 03_extensiones.R
# (a) Heterogeneidad bivariada: I2 de Zhou-Dendukuri y de Holling
# (b) QUADAS-2 por separado
# (c) Metarregresión (bivariada de Reitsma + cribado sobre logDOR)
# (d) Análisis de sensibilidad
# (e) Insumos del funnel plot de Deeks
# (f) Brazos incluidos en la revisión pero no en el metaanálisis conjunto
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(mada); library(metafor); library(writexl)
})
OUT <- "salidas"
DAT <- readRDS(file.path(OUT, "datasets_dta_campy.rds"))
comparaciones <- c("LFA_Culture","EIA_Culture","LFA_PCR","EIA_PCR")
etiquetas <- c(LFA_Culture="ICA/LFA vs Cultivo", EIA_Culture="EIA microplaca vs Cultivo",
               LFA_PCR="ICA/LFA vs PCR", EIA_PCR="EIA microplaca vs PCR")
qd_cols <- c("rob_PS","rob_IT","rob_RS","rob_FT","ac_PS","ac_IT","ac_RS")

fit_biv <- function(d, method = "reml") {
  try(reitsma(d, TP="TP", FN="FN", FP="FP", TN="TN", method = method), silent = TRUE)
}
resumen_biv <- function(f) {
  if (inherits(f, "try-error")) return(rep(NA_real_, 6))
  s <- summary(f); co <- s$coefficients
  c(co["sensitivity","Estimate"], co["sensitivity","95%ci.lb"], co["sensitivity","95%ci.ub"],
    1 - co["false pos. rate","Estimate"], 1 - co["false pos. rate","95%ci.ub"],
    1 - co["false pos. rate","95%ci.lb"])
}

# ---------- (a) I2 bivariados ----------
i2_tab <- bind_rows(lapply(comparaciones, function(cmp) {
  d <- DAT$global[[cmp]]; f <- fit_biv(d); s <- summary(f)
  data.frame(Comparacion = cmp, Etiqueta = etiquetas[[cmp]], k = nrow(d),
             tau2_logitSe = f$Psi[1,1], tau2_logitFPR = f$Psi[2,2],
             cor_aleatoria = s$corRandom[1,2],
             # Valor tal como lo devuelve mada 0.5.12. NO es el I2 de Zhou-Dendukuri
             # que se informa en el manuscrito: mada calcula el componente entre
             # estudios con los errores estándar de los efectos fijos y usa TP+FP /
             # TN+FN como tamaños de grupo. El recalculado está en 03b_i2_zhou.R
             # (columna I2_Zhou_recalc_bivariado de i2_bivariado_completo.csv).
             I2_Zhou_mada_0_5_12 = s$i2$Zhou,
             I2_Holling_noAjust_1 = s$i2$HollingUnadjusted1,
             I2_Holling_noAjust_2 = s$i2$HollingUnadjusted2,
             I2_Holling_noAjust_3 = s$i2$HollingUnadjusted3,
             I2_Holling_ajust_1 = s$i2$HollingAdjusted1,
             I2_Holling_ajust_2 = s$i2$HollingAdjusted2,
             I2_Holling_ajust_3 = s$i2$HollingAdjusted3,
             row.names = NULL)
}))

# ---------- (b) QUADAS-2 ----------
# Codificación de la planilla: 1 = bajo/baja, 2 = alto/alta, 3 = incierto/incierta
niv <- c("1"="Bajo", "2"="Alto", "3"="Incierto")
dom_lab <- c(rob_PS="RS: selección de pacientes", rob_IT="RS: prueba índice",
             rob_RS="RS: estándar de referencia", rob_FT="RS: flujo y tiempo",
             ac_PS="AP: selección de pacientes", ac_IT="AP: prueba índice",
             ac_RS="AP: estándar de referencia")

quadas_estudio <- bind_rows(lapply(comparaciones, function(cmp) {
  d <- DAT$global[[cmp]]
  out <- d[, c("Study_ID","Index_1", qd_cols)]
  out$Comparacion <- cmp; out$Etiqueta <- etiquetas[[cmp]]
  for (q in qd_cols) out[[q]] <- unname(niv[as.character(out[[q]])])
  out[, c("Comparacion","Etiqueta","Study_ID","Index_1", qd_cols)]
}))

quadas_resumen <- bind_rows(lapply(comparaciones, function(cmp) {
  dd <- quadas_estudio %>% filter(Comparacion == cmp)
  bind_rows(lapply(qd_cols, function(q) {
    v <- factor(dd[[q]], levels = c("Bajo","Alto","Incierto"))
    tb <- table(v, useNA = "no")
    data.frame(Comparacion = cmp, Etiqueta = etiquetas[[cmp]],
               Dominio = q, Descripcion = unname(dom_lab[q]),
               k = sum(!is.na(dd[[q]])),
               n_Bajo = as.integer(tb["Bajo"]), n_Alto = as.integer(tb["Alto"]),
               n_Incierto = as.integer(tb["Incierto"]),
               pct_Bajo = round(100*as.integer(tb["Bajo"])/sum(tb), 1),
               pct_Alto = round(100*as.integer(tb["Alto"])/sum(tb), 1),
               pct_Incierto = round(100*as.integer(tb["Incierto"])/sum(tb), 1),
               row.names = NULL)
  }))
}))

# ---------- (c) Metarregresión ----------
mods_cand <- c("preselected","sampling_type","epidemiol_context","clinical_setting",
               "study_population","study_design","prospective","rob_PS","rob_IT","rob_RS",
               "rob_FT","ac_PS","ac_IT","ac_RS","adequate_blinding","blinding_index",
               "high_coi","coi_mfr","region","europe","northamerica","developing_cntry",
               "country_income","result_reading","cutoff_defined","target_species",
               "genus_reactive","delayed_proc","sample_state_EIA","frozen_thawed",
               "simultaneous_fresh","transport_medium","partial_verif","partial_gs",
               "incubation_48h","culture_incubation_h","pcr_type","prevalence",
               "sample_prev_calc","eligibility","year","Index_1")
continuas <- c("culture_incubation_h","sample_prev_calc","year")

calc_logdor <- function(d, cc = 0.5) {
  tp <- d$TP+cc; fn <- d$FN+cc; fp <- d$FP+cc; tn <- d$TN+cc
  d$.logDOR <- log((tp*tn)/(fp*fn))
  d$.var <- 1/tp + 1/fn + 1/fp + 1/tn
  d
}

metareg_una <- function(d, cmp, m) {
  base <- data.frame(Comparacion = cmp, Etiqueta = etiquetas[[cmp]], Moderador = m,
                     k = NA_integer_, niveles = NA_character_, tipo = NA_character_,
                     biv_chi2 = NA_real_, biv_gl = NA_real_, biv_p = NA_real_,
                     dor_QM = NA_real_, dor_gl = NA_real_, dor_p = NA_real_,
                     dor_tau2 = NA_real_, dor_I2 = NA_real_,
                     estado = "no evaluable", row.names = NULL)
  if (!m %in% names(d)) { base$estado <- "variable ausente"; return(base) }
  dm <- d[!is.na(d[[m]]), ]
  if (nrow(dm) < 4) { base$k <- nrow(dm); base$estado <- "k<4 con dato"; return(base) }
  x <- dm[[m]]
  es_cont <- m %in% continuas && is.numeric(x) && length(unique(x)) > 3
  if (es_cont) {
    dm$.mod <- as.numeric(x); tipo <- "continua"; niveles <- "—"
  } else {
    fx <- droplevels(factor(as.character(x)))
    if (nlevels(fx) < 2) { base$k <- nrow(dm); base$niveles <- paste(levels(fx), collapse="/")
                           base$estado <- "sin variación"; return(base) }
    tb <- table(fx)
    if (min(tb) < 2 || nlevels(fx) > 3) {
      base$k <- nrow(dm); base$tipo <- "categórica"
      base$niveles <- paste0(names(tb), "=", as.integer(tb), collapse="; ")
      base$estado <- if (nlevels(fx) > 3) "demasiados niveles para k disponible" else "nivel con k<2"
      return(base)
    }
    dm$.mod <- fx; tipo <- "categórica"
    niveles <- paste0(names(tb), "=", as.integer(tb), collapse="; ")
  }
  base$k <- nrow(dm); base$tipo <- tipo; base$niveles <- niveles

  f0 <- try(reitsma(dm, TP="TP", FN="FN", FP="FP", TN="TN", method="ml"), silent=TRUE)
  f1 <- try(reitsma(dm, formula = cbind(tsens, tfpr) ~ .mod,
                    TP="TP", FN="FN", FP="FP", TN="TN", method="ml"), silent=TRUE)
  if (!inherits(f0,"try-error") && !inherits(f1,"try-error")) {
    a <- try(anova(f1, f0), silent=TRUE)
    if (!inherits(a,"try-error")) {
      base$biv_chi2 <- unname(a$statistic["Chi-squared"])
      base$biv_gl   <- unname(a$statistic["Df"])
      base$biv_p    <- unname(a$statistic["Pr(>Chi-squared)"])
    }
  }
  dd <- calc_logdor(dm)
  rr <- try(rma(yi = .logDOR, vi = .var, mods = ~ .mod, data = dd, method = "REML"), silent=TRUE)
  if (!inherits(rr,"try-error")) {
    base$dor_QM <- rr$QM; base$dor_gl <- rr$m; base$dor_p <- rr$QMp
    base$dor_tau2 <- rr$tau2; base$dor_I2 <- rr$I2
  }
  base$estado <- if (is.na(base$biv_p) && is.na(base$dor_p)) "no convergió" else "estimado"
  base
}

metareg_tab <- bind_rows(lapply(comparaciones, function(cmp) {
  d <- DAT$global[[cmp]]
  bind_rows(lapply(mods_cand, function(m) metareg_una(d, cmp, m)))
}))

# Grupos de moderadores colineales: dentro de una comparación, dos moderadores
# que inducen la misma partición de estudios producen exactamente el mismo modelo.
# La firma incluye el TIPO de modelo además de la partición: un moderador continuo
# y su versión categorizada inducen la misma partición cuando todos los valores son
# distintos, pero no son el mismo modelo (1 gl frente a k-1 gl), de modo que no deben
# agruparse como colineales.
firma_particion <- function(d, m) {
  if (!m %in% names(d)) return(NA_character_)
  x <- d[[m]]
  if (all(is.na(x))) return(NA_character_)
  tipo <- if (m %in% continuas && is.numeric(x) && length(unique(na.omit(x))) > 3) "C" else "F"
  g <- as.integer(factor(as.character(x)))
  g[is.na(x)] <- NA_integer_
  # canonizar: reetiquetar por orden de aparición para que 1/2 y 2/1 coincidan
  lv <- unique(na.omit(g)); paste0(tipo, ":", paste(match(g, lv), collapse = "-"))
}
metareg_tab$grupo_colineal <- NA_character_
for (cmp in comparaciones) {
  d <- DAT$global[[cmp]]
  idx <- which(metareg_tab$Comparacion == cmp)
  fir <- vapply(metareg_tab$Moderador[idx], function(m) {
    z <- firma_particion(d, m); if (is.na(z)) NA_character_ else z
  }, character(1))
  metareg_tab$grupo_colineal[idx] <- ifelse(is.na(fir), NA_character_,
                                            paste0(cmp, "_g", match(fir, unique(fir))))
}
# Representante del grupo: se prefiere un modelo estimable; si hay varios, el primero
# por orden alfabético. Antes se tomaba el primero sin más, de modo que un moderador
# no estimable podía desplazar al modelo estimado de su propio grupo.
metareg_tab <- metareg_tab %>%
  group_by(grupo_colineal) %>%
  mutate(moderadores_equivalentes = if (all(is.na(grupo_colineal))) NA_character_
         else paste(Moderador, collapse = ", "),
         representante = Moderador == Moderador[order(estado != "estimado", Moderador)][1]) %>%
  ungroup()

# Control de falsos descubrimientos (Benjamini-Hochberg). Familia: los moderadores
# estimados y representantes de su grupo DENTRO de cada comparación; el número de
# pruebas de la familia se informa en n_familia_BH porque difiere entre comparaciones.
metareg_tab$biv_q_BH <- NA_real_
metareg_tab$dor_q_BH <- NA_real_
metareg_tab$n_familia_BH <- NA_integer_
for (cmp in comparaciones) {
  sel <- which(metareg_tab$Comparacion == cmp & metareg_tab$estado == "estimado" &
               metareg_tab$representante)
  metareg_tab$n_familia_BH[sel] <- length(sel)
  for (par in c("biv", "dor")) {
    p <- metareg_tab[[paste0(par, "_p")]][sel]
    metareg_tab[[paste0(par, "_q_BH")]][sel] <- p.adjust(p, method = "BH")
  }
}
metarreg_representantes <- metareg_tab %>%
  filter(representante) %>%
  arrange(Comparacion, biv_q_BH, biv_p)
write.csv(metarreg_representantes, file.path(OUT, "metarregresion_representantes.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Subgrupos: Se/Sp por nivel de los moderadores estimados
subgrupos_tab <- bind_rows(lapply(comparaciones, function(cmp) {
  d <- DAT$global[[cmp]]
  ms <- metareg_tab %>% filter(Comparacion == cmp, estado == "estimado", tipo == "categórica")
  bind_rows(lapply(ms$Moderador, function(m) {
    dm <- d[!is.na(d[[m]]), ]; fx <- factor(as.character(dm[[m]]))
    bind_rows(lapply(levels(fx), function(lv) {
      ds <- dm[fx == lv, ]
      est <- if (nrow(ds) >= 4) resumen_biv(fit_biv(ds)) else rep(NA_real_, 6)
      data.frame(Comparacion = cmp, Moderador = m, Nivel = lv, k = nrow(ds),
                 Se_bruta = sum(ds$TP)/sum(ds$TP+ds$FN),
                 Sp_bruta = sum(ds$TN)/sum(ds$TN+ds$FP),
                 Se_biv = est[1], Se_lb = est[2], Se_ub = est[3],
                 Sp_biv = est[4], Sp_lb = est[5], Sp_ub = est[6], row.names = NULL)
    }))
  }))
}))

# ---------- (d) Sensibilidad ----------
escenarios <- function(cmp) {
  d <- DAT$global[[cmp]]; dall <- DAT$all[[paste0(cmp, "_all")]]
  rob4 <- c("rob_PS","rob_IT","rob_RS","rob_FT")
  rob_bajo_todos <- apply(d[, rob4], 1, function(z) all(z == 1, na.rm = FALSE))
  lst <- list(
    "Base (dedup. por ranking de mercado)" = d,
    "Sólo eligibility = 1" = d[!is.na(d$eligibility) & d$eligibility == 1, ],
    "n >= 300" = d[d$n >= 300, ],
    "Riesgo de sesgo bajo en los 4 dominios" = d[rob_bajo_todos, ],
    "Estándar de referencia con RS bajo (rob_RS = 1)" = d[d$rob_RS == 1, ],
    "Todas las entradas (sin deduplicar)" = dall
  )
  # Evidencia sólo accesible por resumen: Tachikawa 2017 (2x2 reconstruido)
  if (!is.null(DAT$sensibilidad) && nrow(DAT$sensibilidad)) {
    ds <- DAT$sensibilidad[DAT$sensibilidad$Comparison == cmp, , drop = FALSE]
    if (nrow(ds)) {
      comunes <- intersect(names(d), names(ds))
      lst[["Añadiendo evidencia reconstruida del resumen (Tachikawa 2017)"]] <-
        rbind(d[, comunes, drop = FALSE], ds[, comunes, drop = FALSE])
    }
  }
  bind_rows(lapply(names(lst), function(nm) {
    ds <- lst[[nm]]
    est <- if (nrow(ds) >= 4) resumen_biv(fit_biv(ds)) else rep(NA_real_, 6)
    data.frame(Comparacion = cmp, Etiqueta = etiquetas[[cmp]], Escenario = nm,
               k = nrow(ds), Se = est[1], Se_lb = est[2], Se_ub = est[3],
               Sp = est[4], Sp_lb = est[5], Sp_ub = est[6], row.names = NULL)
  }))
}
sens_tab <- bind_rows(lapply(comparaciones, escenarios))

loo_tab <- bind_rows(lapply(comparaciones, function(cmp) {
  d <- DAT$global[[cmp]]
  bind_rows(lapply(seq_len(nrow(d)), function(i) {
    est <- resumen_biv(fit_biv(d[-i, ]))
    data.frame(Comparacion = cmp, Etiqueta = etiquetas[[cmp]],
               Excluido = d$Study_ID[i], Prueba = d$Index_1[i], k = nrow(d)-1,
               Se = est[1], Se_lb = est[2], Se_ub = est[3],
               Sp = est[4], Sp_lb = est[5], Sp_ub = est[6], row.names = NULL)
  }))
}))

# ---------- (e) Insumos del funnel plot de Deeks ----------
deeks_datos <- bind_rows(lapply(comparaciones, function(cmp) {
  d <- calc_logdor(DAT$global[[cmp]])
  n1 <- d$TP + d$FN; n0 <- d$FP + d$TN
  data.frame(Comparacion = cmp, Etiqueta = etiquetas[[cmp]],
             Estudio = d$Study_ID, Prueba = d$Index_1,
             lnDOR = d$.logDOR, DOR = exp(d$.logDOR),
             ESS = 4*n1*n0/(n1+n0), inv_sqrt_ESS = 1/sqrt(4*n1*n0/(n1+n0)),
             row.names = NULL)
}))
deeks_fit <- bind_rows(lapply(comparaciones, function(cmp) {
  dd <- deeks_datos %>% filter(Comparacion == cmp)
  if (nrow(dd) < 10) return(data.frame(Comparacion = cmp, intercepto = NA_real_,
                                       pendiente = NA_real_, EE = NA_real_, t = NA_real_,
                                       p = NA_real_, nota = "k<10: prueba no realizada"))
  fit <- lm(lnDOR ~ inv_sqrt_ESS, data = dd, weights = dd$ESS)
  s <- summary(fit)$coefficients
  data.frame(Comparacion = cmp, intercepto = s[1,1], pendiente = s[2,1], EE = s[2,2],
             t = s[2,3], p = s[2,4], nota = "", row.names = NULL)
}))

# ---------- (f) Brazos de la revisión no incluidos en el metaanálisis ----------
dc <- DAT$datos_campy
en_meta <- dc$Record_ID %in% DAT$datos_comparaciones$Record_ID
motivo <- ifelse(dc$solo_sensibilidad == 1,
                 "sólo análisis de sensibilidad: 2x2 reconstruido del resumen",
          ifelse(!is.na(dc$eligibility) & dc$eligibility == 0, "eligibility = 0 (excluido por el equipo)",
          ifelse(!dc$Index_2 %in% c("lateral flow immunoassay","microplate assay"),
                 paste0("plataforma no comparable: ", dc$Index_2),
          ifelse(!dc$Reference_2 %in% c("Culture","PCR"),
                 paste0("estándar de referencia no comparable: ", dc$Reference_2), "otro"))))
revision_no_agrupada <- dc[!en_meta, ] %>%
  transmute(Estudio = Study_ID, Prueba_indice = Index_1, Plataforma = Index_2,
            Referencia = Reference_2, n, TP, FN, FP, TN,
            Sensibilidad = round(sensitivity, 4), Especificidad = round(specificity, 4),
            Motivo = motivo[!en_meta]) %>%
  arrange(Motivo, Estudio)

# ---------- Exportar ----------
write_xlsx(list(I2_bivariado = i2_tab, QUADAS_por_estudio = quadas_estudio,
                QUADAS_resumen = quadas_resumen, metarregresion = metareg_tab,
                subgrupos = subgrupos_tab, sensibilidad = sens_tab,
                leave_one_out = loo_tab, deeks_datos = deeks_datos, deeks_prueba = deeks_fit,
                revision_no_agrupada = revision_no_agrupada),
           file.path(OUT, "resultados_extensiones.xlsx"))
for (nm in c("i2_tab","quadas_estudio","quadas_resumen","metareg_tab","subgrupos_tab",
             "sens_tab","loo_tab","deeks_datos","deeks_fit","revision_no_agrupada")) {
  write.csv(get(nm), file.path(OUT, paste0(nm, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
}
cat("listo\n")
print(i2_tab)
