#!/bin/bash
shmip=F
# ogag or ogag_2
id=ogag
#F1 	-6.0
#F2 	-3.0
#F3 	 0.0
#F4 	+3.0
#F5 	+6.0

dt=( "-6.0" "-3.0" "0.0" "3.0" "6.0" )

for ((i=1 ; i<=5 ; i++))
do
  num=$(printf "%01d\n" $i)
  run=$id"_"$shmip
  name=$shmip$num"_"$id
  echo $name

  WorkPath=/r510/work1/gagliar/SHMIP/$run
  #WorkPath=/Users/ogagliardini/Simulations/ElmerAppli/SHMIP/$run
  echo $WorkPath

  sifName=$WorkPath/$name.sif
  echo $sifName
  ScketchPath=$WorkPath/
  scketch=$ScketchPath"/shmip_"$shmip"_ogag.sif_model"

  cd $WorkPath
  echo ${dt[$i-1]}
  dti=${dt[$i-1]}

  cat $scketch | sed -e "s#<DT>#$dti#g" \
                     -e "s#<name>#$name#g" \
                     -e "s#<num>#$num#g" > $sifName

  echo $name.sif 
  nohup  seqrun -hosts lachouf1 ElmerSolver $name.sif > output$num.out &
done
