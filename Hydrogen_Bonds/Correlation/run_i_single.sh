#!/bin/bash

executable=correlation_i
data=../../../Data

cd Codes/
gfortran $executable.f90 -o $executable
cd ..

n=2

#Read from screen which group of data should we analyse
#echo "Do you want to analyse (1) OH- (2) H+ or (3)all? (Reply a number)"
#read n

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

    dstep=10
    dk=15000
    if [ $prefix = PBE0_vdW ]; then
      dstep=5
      dk=30000
    fi
    
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    sed -i -e "s/xxDSTEPxx/$dstep/" $input
    sed -i -e "s/xxDKxx/$dk/" $input
    ln -s $data/$ion/$prefix/*.trans .
    ./$executable <$input
    #cp $name.corr $data/$ion/$prefix/
    rm -f *.trans
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

  #for prefix in  \
  #      PBE  \
  #      PBE_vdW   \
  #      PBE0_vdW  ;
  #do 
  
  prefix=PBE
  
    dstep=10
    dk=15000
    if [ $prefix = PBE0_vdW ]; then
      dstep=5
      dk=30000
    fi
    
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    sed -i -e "s/xxDSTEPxx/$dstep/" $input
    sed -i -e "s/xxDKxx/$dk/" $input
    #ln -s $data/$ion/$prefix/*.trans .
    ln -s $data/$ion/$prefix/new/*.trans .
    ./$executable <$input
    #cp $name.corr $data/$ion/$prefix/
    rm -f *.trans
    rm $input

  #done

  rm temp
  rm $executable

  cd ..

fi




echo "Analyse completes."
