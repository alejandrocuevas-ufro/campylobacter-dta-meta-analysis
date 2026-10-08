#!/usr/bin/env bash
# Reproduce los cálculos del metaanálisis. Requiere R con los paquetes listados en
# sessionInfo.txt y JAGS para el paso bayesiano. Los resultados se escriben en salidas/.
set -euo pipefail
mkdir -p salidas

for f in R/01_prep_datos.R R/02_analisis_bivariado.R R/03_extensiones.R R/03b_i2_zhou.R \
         R/05_subgrupos_marca.R R/04_bayesiano_jags.R R/06_bayes_priors.R; do
  echo ">>> $f"
  Rscript -e "eval(parse('$f', encoding='UTF-8'))"
done

echo "listo: resultados en salidas/"
