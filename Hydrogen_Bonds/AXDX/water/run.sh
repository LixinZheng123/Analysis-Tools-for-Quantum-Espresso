#!/bin/bash

executable=hydrogen_bond_water
data=../../../../../Data

cd Codes/
gfortran $executable.f90 -o $executable
cd ..

#
#
ion=water

start=1
stop=-1


#if [ $n = 1 ] || [ $n = 3 ]; then
 
  cs=1
  
  cd Data
  
  #Create input.dat from input.template
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  sed -i -e "s/xxCSxx/$cs/" temp
  sed -i -e "s/xxN0xx/$start/" temp
  sed -i -e "s/xxN1xx/$stop/" temp

  for prefix in  \
        PBE_vdW  \
        PBE0_vdW ;
  do 
  
    ln -s $data/$ion/$prefix/*.pos .

    name=$ion"_"$prefix
    input=input-$name.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    #ln -s $data/Water/$name.pos
    ./$executable <$input

    cp $name.bond $data/$ion/$prefix

    rm -f *.pos
    rm $input


  done

  rm temp
  rm $executable

  cd ..

