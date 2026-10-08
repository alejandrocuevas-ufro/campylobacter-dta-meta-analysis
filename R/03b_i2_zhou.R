# ============================================================
# 03b_i2_zhou.R
# I2 de Zhou & Dendukuri (2014) recalculado a partir del modelo bivariado.
#
# Motivo: en mada 0.5.12, summary.reitsma()$i2$Zhou calcula el componente
# "entre estudios" con los errores estándar de los efectos fijos
# (tabfixed[,"Std. Error"]^2) en lugar de las varianzas entre estudios (Psi),
# y usa como tamaños de grupo TP+FP y TN+FN en lugar de TP+FN y FP+TN.
# Aquí se reproduce la misma estructura de la fórmula usando Psi y los
# tamaños de grupo de enfermos / no enfermos. Se informan ambas versiones.
#
# Estructura: I2 = sqrt(det(Sigma_b)) / ( sqrt(det(Sigma_b)) + sqrt(det(Sigma_w)) )
#   det(Sigma_b) = tau_A^2 * tau_B^2 * (1 - rho^2)
#   E[sigma^2_iA] = (exp(tau_A^2/2 + theta_A) + exp(tau_A^2/2 - theta_A) + 2) / n_iA
#   Sigma_w diagonal con el promedio de E[sigma^2_i.] sobre los estudios
# ============================================================

suppressPackageStartupMessages({ library(dplyr); library(mada); library(writexl) })
OUT <- "salidas"
DAT <- readRDS(file.path(OUT, "datasets_dta_campy.rds"))
comparaciones <- c("LFA_Culture","EIA_Culture","LFA_PCR","EIA_PCR")
etiquetas <- c(LFA_Culture="ICA/LFA vs Cultivo", EIA_Culture="EIA microplaca vs Cultivo",
               LFA_PCR="ICA/LFA vs PCR", EIA_PCR="EIA microplaca vs PCR")

exp_within <- function(tau2, theta, n) {
  (exp(tau2/2 + theta) + exp(tau2/2 - theta) + 2) / n
}

i2_zhou_recalc <- function(f, d) {
  Psi <- f$Psi; mu <- as.numeric(coef(f))
  tauA <- Psi[1,1]; tauB <- Psi[2,2]
  rho  <- Psi[1,2]/sqrt(tauA*tauB)
  thetaA <- mu[1]; thetaB <- -mu[2]
  nA <- d$TP + d$FN; nB <- d$FP + d$TN
  wA <- mean(exp_within(tauA, thetaA, nA))
  wB <- mean(exp_within(tauB, thetaB, nB))
  det_b <- tauA*tauB*(1 - rho^2)
  c(I2_biv = sqrt(det_b)/(sqrt(det_b) + sqrt(wA*wB)),
    I2_Se  = tauA/(tauA + wA),
    I2_Sp  = tauB/(tauB + wB),
    E_within_Se = wA, E_within_Sp = wB, rho = rho)
}

i2_comp <- bind_rows(lapply(comparaciones, function(cmp) {
  d <- DAT$global[[cmp]]
  f <- reitsma(d, TP="TP", FN="FN", FP="FP", TN="TN", method = "reml")
  s <- summary(f); z <- i2_zhou_recalc(f, d)
  data.frame(Comparacion = cmp, Etiqueta = etiquetas[[cmp]], k = nrow(d),
             tau2_logitSe = f$Psi[1,1], tau2_logitFPR = f$Psi[2,2],
             cor_aleatoria = z[["rho"]],
             I2_Zhou_recalc_bivariado = z[["I2_biv"]],
             I2_Zhou_recalc_Se = z[["I2_Se"]],
             I2_Zhou_recalc_Sp = z[["I2_Sp"]],
             I2_Zhou_mada = s$i2$Zhou,
             I2_Holling_mada_noAjust_1 = s$i2$HollingUnadjusted1,
             I2_Holling_mada_noAjust_2 = s$i2$HollingUnadjusted2,
             I2_Holling_mada_noAjust_3 = s$i2$HollingUnadjusted3,
             I2_Holling_mada_ajust_1 = s$i2$HollingAdjusted1,
             I2_Holling_mada_ajust_2 = s$i2$HollingAdjusted2,
             I2_Holling_mada_ajust_3 = s$i2$HollingAdjusted3,
             row.names = NULL)
}))

write.csv(i2_comp, file.path(OUT, "i2_bivariado_completo.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(i2_comp[, c("Comparacion","k","cor_aleatoria","I2_Zhou_recalc_bivariado",
                  "I2_Zhou_recalc_Se","I2_Zhou_recalc_Sp","I2_Zhou_mada")])
