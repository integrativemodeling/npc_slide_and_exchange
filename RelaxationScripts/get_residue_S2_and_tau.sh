#!/bin/bash
echo "#resid   Tu[ns] Su2      Tf[ns] Su2*Sf2  Ts[ns] Su2*Sf2*Ss2";  
for i in `seq 252 375`; do
    TMPA=`mktemp`
    TMPB=`mktemp`
    tail -n 2 $i.txt.*.fit7.log  | awk '(NR%4)==2' > $TMPA;
    tail -n 2 $i.txt.*.fit7.log  | awk '(NR%4)==3' > $TMPB;
    paste $TMPA $TMPB |  awk '{n=n+1; tu=tu+$1; tf=tf+$2; ts=ts+$3; Su=Su+$4; Sf=Sf+$5; Ss=Ss+$6; }END{printf "'$i'     %6.2f %5.3f    %6.2f %5.3f    %6.2f %5.3f\n", tu/n,Su/n,tf/n,Su*Sf/n/n,ts/n,Su*Sf*Ss/n/n/n}'; 
done
