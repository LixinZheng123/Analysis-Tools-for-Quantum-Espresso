#!/bin/bash

executable=find_transfer
data=../../../../../Data

cd Codes/
gfortran $executable.f90 -o $executable
cd ..

#n=2

#Read from screen which group of data should we analyse
echo "Do you want to analyse (1) OH- (2) H+ or (3)all? (Reply a number)"
read n

#
#

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
    cp $name.fss $data/$ion/$prefix/
    cp $name.trans $data/$ion/$prefix/
    cp $name.trans_no_r $data/$ion/$prefix/
    cp $name.ion $data/$ion/$prefix/
    cp $name.jump $data/$ion/$prefix/
    if [ $(stat -c %s $name.error) = 0 ]; then
      rm $name.error
    fi
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
  sed -i -e "s/xxNSP1xx/$nsp1/" temp
  sed -i -e "s/xxNSP2xx/$nsp2/" temp
 
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
    ln -s $data/$ion/$prefix/*.pos .
    #ln -s $data/$ion/$prefix/new/*.pos .
    ./$executable <$input
    cp $name.fss $data/$ion/$prefix/
    cp $name.trans $data/$ion/$prefix/
    cp $name.trans_no_r $data/$ion/$prefix/
    cp $name.ion $data/$ion/$prefix/
    cp $name.jump $data/$ion/$prefix/
    if [ $(stat -c %s $name.error) = 0 ]; then
      rm $name.error
    fi
    rm -f *.pos
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi




echo "Analyse completes."
