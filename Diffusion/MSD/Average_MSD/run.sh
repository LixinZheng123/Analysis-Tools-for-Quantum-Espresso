#!/bin/bash

executable=average_msd
data=../../../../../Data

cd Codes/
gfortran $executable.f90 -o $executable
cd ..

#Read from screen which group of data should we analyse
echo "Do you want to analyse (1) OH- (2) H+ or (3)all? (Reply a number)"
read n

#
#

start=1000

if [ $n = 1 ] || [ $n = 3 ]; then
 
  ion=OH-
  nsp1=64
  nsp2=127
  nat=$nsp1+$nsp2
  
  #Create input.dat from input.template
  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  sed -i -e "s/xxNSP1xx/$nsp1/" temp
  sed -i -e "s/xxNSP2xx/$nsp2/" temp
  sed -i -e "s/xxN0xx/$start/" temp

  for prefix in  \
        PBE  \
        PBE_vdW  \
        PBE0_vdW ;
  do 
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    file=$ion"_"$prefix.pos

    cp temp $input
    ln -s $data/$ion/$prefix/*.pos .
    
    nstep=`wc -l ${file} | awk '{print $1/'$nat'}'`

    sed -i -e "s/prefix/$name/" $input
    sed -i -e "s/xxN1xx/$nstep/" $input


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
  nsp1=63
  nsp2=127
  nat=$nsp1+$nsp2
  
  #Create input.dat from input.template
  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  sed -i -e "s/xxNSP1xx/$nsp1/" temp
  sed -i -e "s/xxNSP2xx/$nsp2/" temp
  sed -i -e "s/xxN0xx/$start/" temp

  for prefix in  \
        PBE  \
        PBE_vdW  \
        PBE0_vdW ;
  do 
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    file=$ion"_"$prefix.pos

    cp temp $input
    ln -s $data/$ion/$prefix/*.pos .
    
    nstep=`wc -l ${file} | awk '{print $1/'$nat'}'`

    sed -i -e "s/prefix/$name/" $input
    sed -i -e "s/xxN1xx/$nstep/" $input


    ./$executable <$input

    rm -f *.pos
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi




echo "Analyse completes."
