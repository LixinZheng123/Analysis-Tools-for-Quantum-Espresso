#!/bin/bash

prefix1=PI_PBE_H+_vdw_
prefix2=PI_PBE_OH-_vdw_
surfix1=xyz
#OH- Names
ion1=OH-
cs1=2
#H+ Names
ion2=H+
cs2=3
#Basic parameters
input_temp=INPUT.template
atom11=64
atom12=127
atom21=63
atom22=127
Data_Dir=Data/BOMD/PI_PBE_vdW

#Read from screen which group of data should we analyse
echo "Do you want to analyse (1) OH- (2) H+ or (3)all? (Reply a number)"
read n
#echo "Copying..."

#if [ $n = 1 ] || [ $n = 3 ]; then
  start=1
  #Copy Data
  #cp $GSCRATCH/$Data_Dir/$ion1/$prefix1/*.$surfix1 .
  echo "Copy of OH- completes."
  #Create input.dat from input.template
  #First round
  name=$ion1'_'$prefix1
  Total=$(wc -l <$name.$surfix1)
  stop=$(echo "$Total/$(echo "$atom11+$atom12+2"|bc)"|bc)
  input=input-$name.dat
  cp $input_temp $input
  sed -i -e "s/xxN0xx/$start/" $input
  sed -i -e "s/xxN1xx/$stop/" $input
  sed -i -e "s/prefix/$name/" $input
  sed -i -e "s/xxNSP1xx/$atom11/" $input
  sed -i -e "s/xxNSP2xx/$atom12/" $input
  #sed -i -e "s/xxCSxx/$cs1/" $input
  #Second round
  name=$ion1'_'$prefix2
  Total=$(wc -l <$name.$surfix1)
  stop=$(echo "$Total/$(echo "$atom11+$atom12+1"|bc)"|bc)
  input=input-$name.dat
  cp $input_temp $input
  sed -i -e "s/xxN0xx/$start/" $input
  sed -i -e "s/xxN1xx/$stop/" $input
  sed -i -e "s/prefix/$name/" $input
  sed -i -e "s/xxNSP1xx/$atom11/" $input
  sed -i -e "s/xxNSP2xx/$atom12/" $input
  #sed -i -e "s/xxCSxx/$cs1/" $input
  #Third round
  name=$ion1'_'$prefix3
  Total=$(wc -l <$name.$surfix1)
  stop=$(echo "$Total/$(echo "$atom11+$atom12+1"|bc)"|bc)
  input=input-$name.dat
  cp $input_temp $input
  sed -i -e "s/xxN0xx/$start/" $input
  sed -i -e "s/xxN1xx/$stop/" $input
  sed -i -e "s/prefix/$name/" $input
  sed -i -e "s/xxNSP1xx/$atom11/" $input
  sed -i -e "s/xxNSP2xx/$atom12/" $input
  #sed -i -e "s/xxCSxx/$cs1/" $input
fi

if [[ $n = 2 || $n = 3 ]]; then
  start=1
  cp $GSCRATCH/$Data_Dir/$ion2/$prefix1/*.$surfix1 .
  cp $GSCRATCH/$Data_Dir/$ion2/$prefix1/*.$surfix2 .
  cp $GSCRATCH/$Data_Dir/$ion2/$prefix2/*.$surfix1 .
  cp $GSCRATCH/$Data_Dir/$ion2/$prefix2/*.$surfix2 .
  cp $GSCRATCH/$Data_Dir/$ion2/$prefix3/*.$surfix1 .
  cp $GSCRATCH/$Data_Dir/$ion2/$prefix3/*.$surfix2 .
  echo "Copy of H+ completes."
  #First round
  name=$ion2'_'$prefix2
  Total=$(wc -l <$name.$surfix1)
  stop=$(echo "$Total/$(echo "$atom21+$atom22+1"|bc)"|bc)
  input=input-$name.dat
  cp $input_temp $input
  sed -i -e "s/xxN0xx/$start/" $input
  sed -i -e "s/xxN1xx/$stop/" $input
  sed -i -e "s/prefix/$name/" $input
  sed -i -e "s/xxNSP1xx/$atom21/" $input
  sed -i -e "s/xxNSP2xx/$atom22/" $input
  #sed -i -e "s/xxCSxx/$cs2/" $input
  #Second round
  name=$ion2'_'$prefix1
  Total=$(wc -l <$name.$surfix1)
  stop=$(echo "$Total/$(echo "$atom21+$atom22+1"|bc)"|bc)
  input=input-$name.dat
  cp $input_temp $input
  sed -i -e "s/xxN0xx/$start/" $input
  sed -i -e "s/xxN1xx/$stop/" $input
  sed -i -e "s/prefix/$name/" $input
  sed -i -e "s/xxNSP1xx/$atom21/" $input
  sed -i -e "s/xxNSP2xx/$atom22/" $input
  #sed -i -e "s/xxCSxx/$cs2/" $input
  #Third round
  name=$ion2'_'$prefix3
  Total=$(wc -l <$name.$surfix1)
  stop=$(echo "$Total/$(echo "$atom21+$atom22+1"|bc)"|bc)
  input=input-$name.dat
  cp $input_temp $input
  sed -i -e "s/xxN0xx/$start/" $input
  sed -i -e "s/xxN1xx/$stop/" $input
  sed -i -e "s/prefix/$name/" $input
  sed -i -e "s/xxNSP1xx/$atom21/" $input
  sed -i -e "s/xxNSP2xx/$atom22/" $input
  #sed -i -e "s/xxCSxx/$cs2/" $input
fi
echo "Preparation Completes."
