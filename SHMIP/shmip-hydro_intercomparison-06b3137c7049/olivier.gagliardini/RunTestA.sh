#!/bin/bash
shmip=A
id=ogag
#A1 	7.93e-11
#A2 	1.59e-09
#A3 	5.79e-09
#A4 	2.5e-08
#A5 	4.5e-08
#A6 	5.79e-07

sources=( "7.93e-11" "1.59e-09" "5.79e-9" "2.5e-8" "4.5e-8" "5.79e-7" )

for ((i=1 ; i<=6 ; i++))
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
  echo ${sources[$i-1]}
  source=${sources[$i-1]}

  cat $scketch | sed -e "s#<source>#$source#g" \
                     -e "s#<name>#$name#g" \
                     -e "s#<num>#$num#g" > $sifName

  echo $name.sif 
   nohup  seqrun -hosts lachouf2 ElmerSolver $name.sif > output$num.out &
done
