#!/bin/bash
# Esegue i 18 script uno alla volta, poi il confronto. Log in logs/, riepilogo in logs/status.txt
SS=/mnt/data3/appierdo/AutoimmuneLandscape/Analysis/sex_specific_disease_signature
cd /mnt/data3/appierdo/AutoimmuneLandscape   # root del progetto: .Rprofile attiva renv
: > $SS/logs/status.txt
for f in $SS/scripts/*/*.R; do
  n=$(basename "$f" .R)
  echo "$(date "+%F %T") START $n" >> $SS/logs/status.txt
  Rscript -e "source(\"$f\")" > $SS/logs/$n.log 2>&1
  echo "$(date "+%F %T") END   $n exit=$?" >> $SS/logs/status.txt
done
echo "$(date "+%F %T") START compare" >> $SS/logs/status.txt
Rscript -e "source(\"$SS/compare_sex_specific_disease_signature.R\")" > $SS/logs/compare.log 2>&1
echo "$(date "+%F %T") END   compare exit=$?" >> $SS/logs/status.txt
echo ALL_DONE >> $SS/logs/status.txt
