#!/bin/bash

executable=sdf
data=../../../../Data/330K_CPMD

cd Codes/
gfortran $executable.f90 -o $executable
cd ..

#Read from screen which group of data should we analyse
echo "Do you want to analyse (1) OH- (2) H+ or (3)all? (Reply a number)"
read n

start=1
stop=-1

if [ $n = 1 ] || [ $n = 3 ]; then
 
  echo "What hb number case do you want to analyse? [0 for all, 3 or 4 for specification]"
  read hbnum

  ion=OH-
  cs=2
  
  #Create input.dat from input.template
  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  sed -i -e "s/xxCSxx/$cs/" temp
  sed -i -e "s/xxSTATExx/$state/" temp
  sed -i -e "s/xxN0xx/$start/" temp
  sed -i -e "s/xxN1xx/$stop/" temp
  sed -i -e "s/xxHBxx/$hbnum/" temp

  #for prefix in PBE PBE_vdW PBE0_vdW ;
  for prefix in PBE0_vdW ;
  do 
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    ln -s $data/$ion/$prefix/*.fss .
    ln -s $data/$ion/$prefix/*.hbcase .
    ./$executable <$input
    cp $name.angle $data/$ion/$prefix
    rm -f *.fss
    rm -f *.hbcase
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi

if [[ $n = 2 || $n = 3 ]]; then
  

  hbnum=0
  ion=H+
  cs=3

  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  sed -i -e "s/xxCSxx/$cs/" temp
  sed -i -e "s/xxSTATExx/$state/" temp
  sed -i -e "s/xxN0xx/$start/" temp
  sed -i -e "s/xxN1xx/$stop/" temp
  sed -i -e "s/xxHBxx/$hbnum/" temp
 
 
  for prefix in PBE PBE_vdW PBE0_vdW ;
  do 
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    ln -s $data/$ion/$prefix/*.fss .
    ln -s $data/$ion/$prefix/*.hbcase .
    ./$executable <$input
    cp $name.angle $data/$ion/$prefix
    rm -f *.fss
    rm -f *.hbcase
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi




echo "Analyse completes."
