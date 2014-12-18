#!/bin/bash

executable=correlation_c_no_rattling
data=../../../Data

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

    dstep=10
    dk=6000
    if [ $prefix = PBE0_vdW ]; then
      dstep=5
      dk=12000
    fi
    
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    sed -i -e "s/xxDSTEPxx/$dstep/" $input
    sed -i -e "s/xxDKxx/$dk/" $input
    ln -s $data/$ion/$prefix/*.trans .
    ./$executable <$input
    cp $name.corr_c $data/$ion/$prefix/
    cp $name.trans_no_r $data/$ion/$prefix/
    mv $name.corr_c Continous_Rattling_Excluded
    mv $name.trans_no_r Continous_Rattling_Excluded
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

  for prefix in  \
        PBE  \
        PBE_vdW   \
        PBE0_vdW  ;
  do 
  
  #prefix=PBE
  
    dstep=10
    dk=2500
    if [ $prefix = PBE0_vdW ]; then
      dstep=5
      dk=5800
    fi
    
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    sed -i -e "s/xxDSTEPxx/$dstep/" $input
    sed -i -e "s/xxDKxx/$dk/" $input
    ln -s $data/$ion/$prefix/*.trans .
    ./$executable <$input
    cp $name.corr_c $data/$ion/$prefix/
    cp $name.trans_no_r $data/$ion/$prefix/
    mv $name.corr_c Continous_Rattling_Excluded
    mv $name.trans_no_r Continous_Rattling_Excluded
    rm -f *.trans
    rm $input

  done

  rm temp
  rm $executable

  cd ..

fi




echo "Analyse completes."
