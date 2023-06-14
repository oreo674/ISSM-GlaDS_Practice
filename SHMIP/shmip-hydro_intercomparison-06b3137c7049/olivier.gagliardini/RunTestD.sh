#!/bin/bash
shmip=D
id=ogag
#D1 	-4.0
#D2 	-2.0
#D3 	 0.0
#D4 	+2.0
#D5 	+4.0

dt=( "-4.0" "-2.0" "0.0" "2.0" "4.0" )

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
  nohup  seqrun -hosts lachouf2 ElmerSolver $name.sif > output$num.out &
done
