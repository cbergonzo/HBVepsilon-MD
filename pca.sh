#!/bin/bash

PARM='../nowat.noions.parm7'
REF='../../0BUILD/nowat.noions.01.rst7'
TRAJ0='~/WORK/HBV_rna/2PROD/analysis/nowat.nc.0'
TRAJ1='~/WORK/HBV_rna/2PROD/analysis/nowat.nc.1'
TRAJ2='~/WORK/HBV_rna/2PROD/analysis/nowat.nc.2'
TRAJ3='~/WORK/HBV_rna/2PROD/analysis/nowat.nc.3'
TRAJ4='~/WORK/HBV_rna/2PROD/analysis/nowat.nc.4'
TRAJ5='~/WORK/HBV_rna/2PROD/analysis/nowat.nc.5'
TRAJ6='~/WORK/HBV_rna/2PROD/analysis/nowat.nc.6'
TRAJ7='~/WORK/HBV_rna/2PROD/analysis/nowat.nc.7'
TRAJ8='~/WORK/HBV_rna/2PROD/analysis/nowat.nc.8'
TRAJ9='~/WORK/HBV_rna/2PROD/analysis/nowat.nc.9'

TRAJ10='../combined.NMR03.01.nc'
TRAJ11='../combined.NMR03.02.nc'
TRAJ12='../combined.NMR03.03.nc'

cpptraj <<EOF
parm $PARM
trajin $TRAJ0 1 last 20
trajin $TRAJ1 1 last 20
trajin $TRAJ2 1 last 20
trajin $TRAJ3 1 last 20
trajin $TRAJ4 1 last 20
trajin $TRAJ5 1 last 20
trajin $TRAJ6 1 last 20
trajin $TRAJ7 1 last 20
trajin $TRAJ8 1 last 20
trajin $TRAJ9 1 last 20
trajin $TRAJ10 1 last 2
trajin $TRAJ11 1 last 2
trajin $TRAJ12 1 last 2


# RMS-fit to first frame
rms first :1-61&!@H=

# Create an average structure
average Avg.fullRNA.rst7 ncrestart

# Save coordinates as 'crd1'
createcrd crd1
run

# Fit to average structure
reference Avg.fullRNA.rst7 [avg]

# RMS-fit to average structure
crdaction crd1 rms ref [avg] :1-61&!@H=

# Calculate coordinate covariance matrix
crdaction crd1 matrix covar :1-61&!@H= name gaccCovar

# Diagonalize coordinate covariance matrix, first 15 E.vecs
runanalysis diagmatrix gaccCovar name gaccEvecs out evecs.fullRNA.dat vecs 15

# Now create separate projections for loop region for each trajectory
#
crdaction crd1 projection P1 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 1,5000
crdaction crd1 projection P2 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 5001,10000
crdaction crd1 projection P3 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 10001,15000
crdaction crd1 projection P4 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 15001,20000
crdaction crd1 projection P5 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 20001,25000
crdaction crd1 projection P6 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 25001,30000
crdaction crd1 projection P7 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 30001,35000
crdaction crd1 projection P8 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 35001,40000
crdaction crd1 projection P9 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 40001,45000
crdaction crd1 projection P10 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 45001,50000
crdaction crd1 projection P11 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 50001,55000
crdaction crd1 projection P12 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 55001,60000
crdaction crd1 projection P13 modes gaccEvecs \
beg 1 end 15 :1-61&!@H= crdframes 60001,65000

dataset makexy P1:1 P1:2 name pca112
writedata plot.run1.12.agr pca112
dataset makexy P2:1 P2:2 name pca212
writedata plot.run2.12.agr pca212
dataset makexy P3:1 P3:2 name pca312
writedata plot.run3.12.agr pca312
dataset makexy P4:1 P4:2 name pca412
writedata plot.run4.12.agr pca412
dataset makexy P5:1 P5:2 name pca512
writedata plot.run5.12.agr pca512
dataset makexy P6:1 P6:2 name pca612
writedata plot.run6.12.agr pca612
dataset makexy P7:1 P7:2 name pca712
writedata plot.run7.12.agr pca712
dataset makexy P8:1 P8:2 name pca812
writedata plot.run8.12.agr pca812
dataset makexy P9:1 P9:2 name pca912
writedata plot.run9.12.agr pca912
dataset makexy P10:1 P10:2 name pca1012
writedata plot.run10.12.agr pca1012
dataset makexy P11:1 P11:2 name pca1112
writedata plot.run11.12.agr pca1112
dataset makexy P12:1 P12:2 name pca1212
writedata plot.run12.12.agr pca1212
dataset makexy P13:1 P13:2 name pca1312
writedata plot.run13.12.agr pca1312


# Now histogram first 5 projections for each
#
hist P1:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P1-1
hist P1:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P1-2
hist P1:3,*,*,*,100 out pca.hist.fullRNA.agr norm name P1-3
hist P2:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P2-1
hist P2:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P2-2
hist P2:3,*,*,*,100 out pca.hist.fullRNA.agr norm name P2-3
hist P3:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P3-1
hist P3:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P3-2
hist P3:3,*,*,*,100 out pca.hist.fullRNA.agr norm name P3-3
hist P4:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P4-1
hist P4:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P4-2
hist P5:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P5-1
hist P5:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P5-2
hist P5:3,*,*,*,100 out pca.hist.fullRNA.agr norm name P5-3
hist P6:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P6-1
hist P6:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P6-2
hist P6:3,*,*,*,100 out pca.hist.fullRNA.agr norm name P6-3
hist P7:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P7-1
hist P7:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P7-2
hist P7:3,*,*,*,100 out pca.hist.fullRNA.agr norm name P7-3
hist P8:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P8-1
hist P8:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P8-2
hist P8:3,*,*,*,100 out pca.hist.fullRNA.agr norm name P8-3
hist P9:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P9-1
hist P9:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P9-2
hist P9:3,*,*,*,100 out pca.hist.fullRNA.agr norm name P9-3
hist P10:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P10-1
hist P10:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P10-2
hist P10:3,*,*,*,100 out pca.hist.fullRNA.agr norm name P10-3
hist P11:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P11-1
hist P11:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P11-2
hist P11:3,*,*,*,100 out pca.hist.fullRNA.agr norm name P11-3
hist P12:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P12-1
hist P12:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P12-2
hist P12:3,*,*,*,100 out pca.hist.fullRNA.agr norm name P12-3
hist P13:1,*,*,*,100 out pca.hist.fullRNA.agr norm name P13-1
hist P13:2,*,*,*,100 out pca.hist.fullRNA.agr norm name P13-2
hist P13:3,*,*,*,100 out pca.hist.fullRNA.agr norm name P13-3

EOF

