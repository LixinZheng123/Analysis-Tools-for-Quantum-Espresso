#!/bin/bash

executable=wfc_rdf
data=../../../../Calculation/330K_CPMD

cd Codes/
gfortran $executable.f90 -o $executable
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
  nsp1=64
  nsp2=127
  
  #Create input.dat from input.template
  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  sed -i -e "s/xxCSxx/$cs/" temp
  sed -i -e "s/xxN0xx/$start/" temp
  sed -i -e "s/xxN1xx/$stop/" temp
  sed -i -e "s/xxNSP1xx/$nsp1/" temp
  sed -i -e "s/xxNSP2xx/$nsp2/" temp

  for prefix in  \
        PBE  \
        PBE_vdW  \
        PBE0_vdW ;
  do 
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    ln -s $data/$ion/$prefix/*.pos .
    ./$executable <$input
    rm -f *.pos
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi

if [[ $n = 2 || $n = 3 ]]; then
  

  ion=H+
  cs=3
  nsp1=63
  nsp2=127

  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  sed -i -e "s/xxCSxx/$cs/" temp
  sed -i -e "s/xxN0xx/$start/" temp
  sed -i -e "s/xxN1xx/$stop/" temp
  sed -i -e "s/xxNSP1xx/$nsp1/" temp
  sed -i -e "s/xxNSP2xx/$nsp2/" temp
 
 
  for prefix in  \
        PBE  \
        PBE_vdW  \
        PBE0_vdW ;
  do 
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    ln -s $data/$ion/$prefix/*.pos .
    ./$executable <$input
    rm -f *.pos
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi




echo "Analyse completes."
