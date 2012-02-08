#!/bin/sh
# Antiprism Example File - http://www.antiprism.com
# This file may be copied, modified and redistributed
#
# Make a dual pair - Roger Kaufman

if [ $# -eq 0 ]
then
        rm -f tmp_base.off
        while read data; do
            echo "$data" >> tmp_base.off
        done
else
        cp $1 tmp_base.off
fi

pol_recip tmp_base.off | off_color -f royalblue > tmp_dual.off
off_util tmp_base.off tmp_dual.off > tmp_base_dual.off

cat tmp_base_dual.off

rm -f tmp_base.off tmp_dual.off tmp_base_dual.off


