#!/bin/bash
shmip=C
id=ogag
#C1 	0.25
#C2 	0.50
#C3 	1.00
#C4 	2.00

ra=( "0.25" "0.50" "1.00" "2.00" )

for ((i=1 ; i<=4 ; i++))
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
  echo ${ra[$i-1]}
  rai=${ra[$i-1]}

  cat $scketch | sed -e "s#<ra>#$rai#g" \
                     -e "s#<name>#$name#g" \
                     -e "s#<num>#$num#g" > $sifName

  echo $name.sif 
   nohup  seqrun -hosts lachouf2 ElmerSolver $name.sif > output$num.out &
done
