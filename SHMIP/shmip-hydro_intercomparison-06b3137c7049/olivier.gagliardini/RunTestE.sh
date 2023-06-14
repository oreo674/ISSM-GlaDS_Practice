#!/bin/bash
shmip=E
# ogag or ogag_2
id=ogag
#E1      0.05 	Bench Glacier reference geometry
#E2 	0 	 
#E3 	-0.1 	Starting to have an overdeepening
#E4 	-0.5 	Overdeepening around supercooling threshold
#E5 	-0.7

paras=( "0.05" "0.0" "-0.1" "-0.5" "-0.7" )

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
  echo ${paras[$i-1]}
  para=${paras[$i-1]}

  cat $scketch | sed -e "s#<para>#$para#g" \
                 -e "s#<name>#$name#g" \
		 -e "s#<num>#$num#g" > $sifName

  echo $name.sif 
# echo $name.sif > ELMERSOLVER_STARTINFO
  nohup  seqrun -hosts lachouf1 ElmerSolver $name.sif > output$num.out &
done
