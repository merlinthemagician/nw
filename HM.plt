#!/bin/sh
if [ $1 ]; then
    TF=$1
else
    TF=25
fi

if [ $2 ]; then
    SCL=$2
else
    SCL=1.2
fi

gnuplot <<EOF
set term aqua enhanced

title="TF = ${TF} pM"

if($3) {
title = sprintf("%s, with TFPI", title)
}
if($4) {
title = sprintf("%s, with ATIII", title)
}

set title title

set grid

set xlabel "time [s]"
set ylabel "Concentrations [{/Symbol m}M]"

set yrange [:1800]

plot "< ./HockinMann_JBC2002 $1 $3 $4" u 1:8 t column(8) w l, "" u 1:26 t column(26) w l, "" u 1:(column(8)+${SCL}*column(26)) t "IIa + ${SCL}mIIa" w l
EOF
