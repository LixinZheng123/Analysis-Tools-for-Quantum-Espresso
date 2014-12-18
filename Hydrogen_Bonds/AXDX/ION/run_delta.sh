#!/bin/bash

executable=hydrogen_bond_different_state
data=../../../../../Data/330K_CPMD

cd Codes/
gfortran $executable.f90 -o $executable
cd ..

#Read from screen which group of data should we analyse
echo "Do you want to analyse (1) OH- (2) H+ or (3)all? (Reply a number)"
read n

#
#

stop=-1

rm bond_report.dat

if [ $n = 1 ] || [ $n = 3 ]; then
 
  ion=OH-
  cs=2
  
  #Create input.dat from input.template
  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  sed -i -e "s/xxCSxx/$cs/" temp
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
    ln -s $data/$ion/$prefix/*.sigma .
    ln -s $data/$ion/$prefix/*.trans .
    ./$executable <$input

    cp $name.hbcase $data/$ion/$prefix

    rm -f *.fss
    rm -f *.sigma
    rm -f *.trans
    rm $input

    echo $name >> ../bond_report.dat
    paste $name.bond  >> ../bond_report.dat
    echo ====== >> ../bond_report.dat
    echo >> ../bond_report.dat

    rm *.bond

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
    ln -s $data/$ion/$prefix/*.sigma .
    ln -s $data/$ion/$prefix/*.trans .
    ./$executable <$input
    
    cp $name.hbcase $data/$ion/$prefix
    
    rm -f *.fss
    rm -f *.sigma
    rm -f *.trans
    rm $input

    echo $name >> ../bond_report.dat
    paste $name.bond  >> ../bond_report.dat
    echo ====== >> ../bond_report.dat
    echo >> ../bond_report.dat
    
    rm *.bond

  done

  rm temp
  rm $executable

  cd ..

fi

cp bond_report.dat ../../../Data


