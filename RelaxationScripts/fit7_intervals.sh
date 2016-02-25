#!/bin/bash
MATLAB='/Applications/MATLAB_R2013a.app/bin/matlab'
f1=$2
f2=$3
f3=$4
f4=$5
MEAN="((\$$f1+\$$f2+\$$f3+\$$f4)/4.0)"
FITFILE=$1.$f1.$f2.$f3.$f4.fit7;
PSFILE=$1.$f1.$f2.$f3.$f4.fit7.ps
R1R2FILE=$1.$f1.$f2.$f3.$f4.fit7.R1R2
if [ -e $FITFILE ] && [ -e $R1R2FILE ] ; then
    echo $FITFILE and $R1R2FILE already exist
    exit 0;
fi

echo Fitting
# if [ ! -e $FITFILE ] || [ 1 ] ; then
if [ ! -e $FITFILE ]  ; then
    gnuplot -p >& $1.$f1.$f2.$f3.$f4.fit7.log <<EOF
# GT and BT01 are proxies for scanning a range of values of x that produe
# images between x0 and 1000.0 (for GT) or between 0.0 and 1.0 (for BT01)
GT(x,x0) = (1000.0-x0)/pi*(atan(x)+pi/2)+x0
BT01(x) = (atan(x)+pi/2)/pi

# See eq. 4 in Pfeifer et al., 2001 ("Simulated and NMR-derived..."):  
# f(t) is sum of three decaying exponentials with three different time-scales 
# of the form: 
# f(t)=S_MD^2 + Cu*exp(-t/tau_u) + Cf*exp(-t/tau_f) + Cs*exp(-t/tau_s) 
# The scaling coefficients of each exponent are also defined:
# Cu=(1-Su^2); Cf=Su^2*(1-Sf^2); Cs=Su^2*Sf^2)*(1-Ss^2)
# Su^2, St^2, Ss^2 and S_MD^2 are the squaredorder parameters.
# The parameters satisfy: 
#    S_MD^2=Su^2*St^2*Ss^2 
#    0 < tau_u < tau_f < tau_s
#    0 < Si^2 < 1 for all order parameters S_i
f(t) = BT01(s1)*BT01(s2)*BT01(s3) + (1-BT01(s1))*exp(-t/GT(t1,0)) + BT01(s1)*(1-BT01(s2))*exp(-t/GT(t2,GT(t1,0))) + BT01(s1)*BT01(s2)*(1-BT01(s3))*exp(-t/GT(t3,GT(t2,GT(t1,0))))

t1=-0.1; t2=0.01; t3=0.1;
s1=0.1; s2=0.01; s3=0.1;
fit f(x) '$1' u 1:($MEAN) via s1,t1,s2,t2,s3,t3

# Print tau1 (ultrafast), tau2 (fast) and tau3 (slow)
pr GT(t1,0),GT(t2,GT(t1,0)),GT(t3,GT(t2,GT(t1,0)))
# Print Sigma-1^2 (ultrafast), Sigma-2^2 (fast) and Sigma-3^2 (slow)
pr BT01(s1),BT01(s2),BT01(s3)

#set terminal postscript color
#set output '$PSFILE'
#plot [0:300][]f(x) w l, '$1' u 1:($MEAN) w p pt 6;
#unset output
#`ps2pdf $PSFILE`

set samples 300001
set table '$FITFILE'
plot [0:300]f(x)
unset table

EOF
BASE=`basename $0`
TMPFILE=`mktemp -t $BASE`;
awk '{if (NR>=5 && NF>=2) {print $1,$2;}}' $FITFILE > $TMPFILE
mv $TMPFILE $FITFILE
fi
echo Done fitting


echo Matlabbing
tmux has-session -t matlab;
if [ $? == 1 ]; then
    tmux new -d -s matlab "$MATLAB -nodesktop -nojvm -nosplash -nodisplay" 
fi
PWD=`pwd`
bash mx_matlab.sh <<EOF
  fit7_R1R2_part ('$PWD', '$FITFILE', '$R1R2FILE');
EOF


