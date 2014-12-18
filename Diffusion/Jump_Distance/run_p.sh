#!/bin/bash

executable=plot
data=../../../../Data

cd Codes/
gfortran $executable.f90 -o $executable
cd ..

#Read from screen which group of data should we analyse
echo "Do you want to analyse (1) OH- (2) H+ or (3)all? (Reply a number)"
read n

#
#

if [ $n = 1 ] || [ $n = 3 ]; then
 
  ion=OH-
  
  #Create input.dat from input.template
  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  for prefix in  \
        PBE  \
        PBE_vdW  \
        PBE0_vdW ;
  do 

    start=
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    ./$executable <$input
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi

if [[ $n = 2 || $n = 3 ]]; then
  

  ion=H+

  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  for prefix in  \
        PBE  \
        PBE_vdW   \
        PBE0_vdW  ;
  do 
  
  #prefix=PBE

    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    ./$executable <$input
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi




echo "Analyse completes."
