#!/bin/bash

executable1=new_bond_length
executable2=new_bond_length_H+
data=../../../../Data/One_Thermostat_330K_CPMD_Zhaofeng

cd Codes/
gfortran $executable1.f90 -o $executable1
gfortran $executable2.f90 -o $executable2
cd ..

#Read from screen which group of data should we analyse
echo "Do you want to analyse (1) OH- (2) H+ or (3)all? (Reply a number)"
read n

#
#

start=1
stop=-1


if [ $n = 1 ] || [ $n = 3 ]; then
 
  ion=OH-
  cs=2
  
  #Create input.dat from input.template
  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable1 .

  sed -i -e "s/xxCSxx/$cs/" temp
  sed -i -e "s/xxN0xx/$start/" temp
  sed -i -e "s/xxN1xx/$stop/" temp

  for prefix in  \
        PBE  \
        PBE_vdW  \
        PBE0_vdW ;
  do 
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    ln -s $data/$ion/$prefix/*.fss .
    ln -s $data/$ion/$prefix/*.trans .
    ./$executable1 <$input

    cp $name.sigma $data/$ion/$prefix

    rm -f *.fss
    rm -f *.trans
    rm $input

  done

  rm temp
  rm $executable1

  cd ..

fi

if [[ $n = 2 || $n = 3 ]]; then
  

  ion=H+
  cs=3

  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable2 .

  sed -i -e "s/xxCSxx/$cs/" temp
  sed -i -e "s/xxN0xx/$start/" temp
  sed -i -e "s/xxN1xx/$stop/" temp
 
 
  for prefix in  \
        PBE  \
        PBE_vdW  \
        PBE0_vdW ;
  do 
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    ln -s $data/$ion/$prefix/*.fss .
    ln -s $data/$ion/$prefix/*.trans .
    ./$executable2 <$input
    
    cp $name.sigma $data/$ion/$prefix
    
    rm -f *.fss
    rm -f *.trans
    rm $input

  done

  rm temp
  rm $executable2

  cd ..

fi




