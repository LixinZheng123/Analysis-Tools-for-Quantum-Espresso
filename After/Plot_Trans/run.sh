#!/bin/bash

executable=plot
data=../../../Data

cd Codes/
gfortran $executable.f90 -o $executable
cd ..


#Read from screen which group of data should we analyse
echo "Do you want to analyse (1) OH- (2) H+ or (3)all? (Reply a number)"
read n

#
#

start=1

if [ $n = 1 ] || [ $n = 3 ]; then
 
  ion=OH-
  cs=2
  
  #Create input.dat from input.template
  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  sed -i -e "s/xxCSxx/$cs/" temp

  for prefix in  \
        PBE  \
        PBE_vdW  \
        PBE0_vdW ;
  do 
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    ln -s $data/$ion/$prefix/*.trans .
    ./$executable <$input
    cp trans_$name.dat $data/$ion/$prefix/
    rm -f *.trans
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi

if [[ $n = 2 || $n = 3 ]]; then
  

  ion=H+
  cs=3

  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  sed -i -e "s/xxCSxx/$cs/" temp
 
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
    ln -s $data/$ion/$prefix/*.trans .
    #ln -s $data/$ion/$prefix/new/*.pos .
    ./$executable <$input
    cp trans_$name.dat $data/$ion/$prefix/
    rm -f *.trans
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi




echo "Analyse completes."
