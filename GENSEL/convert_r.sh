#!/bin/bash
#awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' run.mrkRes1 | sed -e 's/_/ /g' | sed '1 s/-/ /2' | sed '1 s/map/chr /' | sed '1 s/posn/position/' | awk '{print $11,$12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > 1_run.mrkRes1_pp14pi99

awk '{print $1,$4,$5,$6,$7,$8,$9}' run.winQTL1 | sed -e 's/_/ /g' | sed '1 s/-/ /2' | sed '1 s/map/chr /' | sed '1 s/pos0/position/' |  sed '1 s/%//g' |  sed '1 s/>//g' |  sed '1 s/#//' | awk '{ print $7,$8,$2,$3,$4,$5,$6}' | sort -n -k1 -k2 > 1_run.winQTL_pp14pi99

#awk '{print $1,$4,$5,$6,$7,$8,$9}' run.winQTL2 | sed -e 's/_/ /g' | sed '1 s/-/ /2' | sed '1 s/map/chr /' | sed '1 s/pos0/position/' |  sed '1 s/%//g' |  sed '1 s/>//g' |  sed '1 s/#//' | awk '{ print $7,$8,$2,$3,$4,$5,$6}' | sort -n -k1 -k2 > 1_run2.winQTL_pp14pi99
