#!/bin/bash

executable=ion_wfc_rdf
data=../../../../Data/330K_CPMD

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
 
  echo "What hb number case do you want to analyse? [0 for all, 3 or 4 for specification]"
  read hbnum

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
  sed -i -e "s/xxHBxx/$hbnum/" temp

  #for prefix in PBE ;
  for prefix in PBE PBE_vdW ;
  do 
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    ln -s $data/$ion/$prefix/*.ion .
    ln -s $data/$ion/$prefix/*.wfc .
    ln -s $data/$ion/$prefix/*.hbcase .
    ./$executable <$input
    rm -f *.ion
    rm -f *.wfc
    rm -f *.hbcase
    rm $input

  done

  for prefix in PBE0_vdW ;
  do 
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    ln -s $data/$ion/$prefix/*.ion .
    ln -s $data/$ion/$prefix/*.wfc.copy $name.wfc
    ln -s $data/$ion/$prefix/*.hbcase .
    ./$executable <$input
    rm -f *.ion
    rm -f *.wfc
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
  sed -i -e "s/xxHBxx/$hbnum/" temp
 
 
  for prefix in PBE PBE_vdW PBE0_vdW ;
  do 
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    ln -s $data/$ion/$prefix/*.ion .
    ln -s $data/$ion/$prefix/*.wfc .
    ln -s $data/$ion/$prefix/*.hbcase .
    ./$executable <$input
    rm -f *.ion
    rm -f *.wfc
    rm -f *.hbcase
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi




echo "Analyse completes."
