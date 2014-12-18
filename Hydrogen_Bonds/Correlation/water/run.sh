#!/bin/bash

executable=correlation
data=../../../../Data

cd Codes/
gfortran $executable.f90 -o $executable
cd ..

#
#
ion=water
  
  cd Data

  #Create input.dat from input.template
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  for prefix in  \
        PBE_vdW  \
        PBE0_vdW ;
  do

    dstep=10
    dk=7000
    if [ $prefix = PBE0_vdW ]; then
      dstep=5
      dk=10000
    fi

    name=$ion"_"$prefix
    input=input-$name.dat

    ln -s $data/$ion/$prefix/*.bond .

    Total=$(wc -l <$name.bond)
    stop=$(echo "$Total/$(echo "64+1"|bc)"|bc)
    
    cp temp $input
    sed -i -e "s/prefix/$name/" $input
    sed -i -e "s/xxDSTEPxx/$dstep/" $input
    sed -i -e "s/xxDKxx/$dk/" $input
    sed -i -e "s/xxSTOPxx/$stop/" $input
    ./$executable <$input
    rm -f *.bond
    #rm $input

  done
 
  rm temp
  rm $executable

  cd ..

echo "Analyse completes."
