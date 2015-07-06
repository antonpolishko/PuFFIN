#!/bin/sh
samtools view -f 67 -F 1804 $1 |  awk '{if ($9>0) print $3, $4, $4+$9, $5; if ($9<0) print $3, $8, $8-$9, $5}'