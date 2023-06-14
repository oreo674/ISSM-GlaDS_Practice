#!/bin/bash
shmip=B
id=ogag
#B1 	1 	B1_M.csv
#B2 	10 	B2_M.csv
#B3 	20 	B3_M.csv
#B4 	50 	B4_M.csv
#B5 	100 	B5_M.csv

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

  cat $scketch | sed -e "s#<name>#$name#g" \
                     -e "s#<num>#$num#g" > $sifName

  echo $name.sif 
  nohup  seqrun -hosts lachouf2 ElmerSolver $name.sif > output$num.out &
done
