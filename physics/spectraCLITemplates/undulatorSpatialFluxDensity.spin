TOTAL_PROCESSES 1

@ProcessNo 0 7 -1 0

[OUTPUTDATA]
Name	<outputRoot>

[ACCELERATOR]
eGeV	<beamEnergy-GeV>
imA	<beamCurrent-mA>
emitt	<ex0>
coupl	<kappa>
espread	<Sdelta0>
betax	<betax>
alphax	<alphax>
betay	<betay>
alphay	<alphay>
eta	<etax>
deta	<etaxp>
etay	<etay>
detay	<etayp>

[SOURCE]
lu	<period-cm>
length	<length-m>
ky	<Ky>
kx	<Kx>
type    <sourceType>

[CALCULATION]
Name	Untitled
slit_dist	<pinholeDistance>
fixep	<EpFixed-eV>
min_x	0
max_x	<xMax-mm>
meshx	<nx>
min_y	0
max_y	<yMax-mm>
meshy	<ny>
tgtharm	1
slitvarx	2
slitvary	2
fin_dist	1
auto	1
autoratio	0.5
epdiv	-1
accuracy	1
energy_conv	0
zeroemitt	0
zeroespread	0
header	1
unitlabel	1
suffix_sd_mesh	spout

