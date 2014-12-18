#!/bin/bash

executable=exp_O_msd
data=../../../../../Data/One_Thermostat_330K_CPMD

cd Codes/
gfortran $executable.f90 -o $executable
cd ..

#Read from screen which group of data should we analyse
echo "Do you want to analyse (1) OH- (2) H+ or (3)all? (Reply a number)"
read n

#
#

start=2000

if [ $n = 1 ] || [ $n = 3 ]; then
 
  ion=OH-
  nsp1=64
  nsp2=127
  cs=2
  
  #Create input.dat from input.template
  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  sed -i -e "s/xxN0xx/$start/" temp
  sed -i -e "s/xxCSxx/$cs/" temp

  for prefix in PBE PBE_vdW ;
  #for prefix in PBE ;
  do 


    dstep=10
    dk=10000
    if [ $prefix = PBE0_vdW ]; then
      dstep=5
      dk=40000
    fi
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    file=$ion"_"$prefix.ion

    cp temp $input
    ln -s $data/$ion/$prefix/*.ion .
    

    sed -i -e "s/prefix/$name/" $input
    sed -i -e "s/xxdstepxx/$dstep/" $input
    sed -i -e "s/xxDKxx/$dk/" $input


    ./$executable <$input

    cp $name.msd $data/$ion/$prefix/

    if [ $(stat -c %s $name.msd_error) = 0 ]; then
      rm $name.msd_error
    fi

    rm -f *.ion
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
  
  #Create input.dat from input.template
  cd $ion
  cp ../INPUT.template temp
  cp ../Codes/$executable .

  sed -i -e "s/xxN0xx/$start/" temp
  sed -i -e "s/xxCSxx/$cs/" temp

  for prefix in PBE PBE_vdW ;
  do 
  
  #prefix=PBE

    dstep=10
    dk=10000
    if [ $prefix = PBE0_vdW ]; then
      dstep=5
      dk=40000
    fi
  
    name=$ion"_"$prefix
    input=input-$ion"_"$prefix.dat
    file=$ion"_"$prefix.ion

    cp temp $input
    ln -s $data/$ion/$prefix/*.ion .
    #ln -s $data/$ion/$prefix/new/*.ion .
    #ln -s $data/$ion/$prefix/new/*.jump .
    

    sed -i -e "s/prefix/$name/" $input
    sed -i -e "s/xxdstepxx/$dstep/" $input
    sed -i -e "s/xxDKxx/$dk/" $input


    ./$executable <$input
    
    cp $name.msd $data/$ion/$prefix/
    #cp $name.msd $data/$ion/$prefix/new

    if [ $(stat -c %s $name.msd_error) = 0 ]; then
      rm $name.msd_error
    fi

    rm -f *.ion
    rm $input

  done


  rm temp
  rm $executable

  cd ..

fi




echo "Analyse completes."
