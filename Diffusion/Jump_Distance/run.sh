#!/bin/bash

executable=jump_distance
data=../../../../Data

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
  
  #Create input.dat from input.template
  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  for prefix in  \
        PBE  \
        PBE_vdW  \
        PBE0_vdW ;
  do 

    ln -s $data/$ion/$prefix/*.ion .
    ln -s $data/$ion/$prefix/*.trans .
    ln -s $data/$ion/$prefix/*.trans_no_r .
    name=$ion"_"$prefix
    total=$(wc -l <$name.trans)
    total2=$(wc -l <$name.trans_no_r)
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    sed -i -e "s/xxSTOPxx/$total/" $input
    sed -i -e "s/xxSTOP2xx/$total2/" $input
    ./$executable <$input
    rm -f *.ion
    rm -f *.trans
    rm -f *.trans_no_r
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

  for prefix in  \
        PBE  \
        PBE_vdW   \
        PBE0_vdW  ;
  do 
  
  #prefix=PBE

    ln -s $data/$ion/$prefix/*.ion .
    ln -s $data/$ion/$prefix/*.trans .
    ln -s $data/$ion/$prefix/*.trans_no_r .
    name=$ion"_"$prefix
    total=$(wc -l <$name.trans)
    total2=$(wc -l <$name.trans_no_r)
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    sed -i -e "s/xxSTOPxx/$total/" $input
    sed -i -e "s/xxSTOP2xx/$total2/" $input
    ./$executable <$input
    #cp $name.corr $data/$ion/$prefix/
    rm -f *.ion
    rm -f *.trans
    rm -f *.trans_no_r
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi




echo "Analyse completes."
