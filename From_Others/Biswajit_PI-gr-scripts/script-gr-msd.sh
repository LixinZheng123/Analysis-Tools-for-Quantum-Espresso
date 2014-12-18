#!/bin/bash

# compile the code 
 ifort -O3 compute_gOOr.f90 -o compute_gOOr.out 
 ifort -O3 compute_gOHr.f90 -o compute_gOHr.out 
 ifort -O3 compute_gHHr.f90 -o compute_gHHr.out 
 ifort -O3 compute_msd.f90 -o compute_msd.out


#for dir in 64W-PI-VDWPBE0-BO-330K \
for dir in  \
          64W-PI-PBE-BO-300K \
          64W-PI-PBE-BO-330K \
          64W-PI-VDWPBE-BO-300K \
          64W-PI-VDWPBE-BO-330K ;
do

  if [ ! -d $dir ]; then
    mkdir $dir
  fi
  
  cd $dir/

  nskip=0 # number of configurations to skip from pos files ...
  
  nat=194 # number of lines for each step ...
  
  for ((i=0; i<=7; i+=1)) ; do
  
   ln -s ../../$dir/data.pos_${i}.xyz .
  
   file=data.pos_${i}.xyz
  
   nstep=`wc -l ${file} | awk '{print $1/'$nat'}'`
   echo $nstep
  
   #
   # computing g(r)s ....
   #
   echo $file > input 
   echo $nstep >> input 
   echo $nskip >> input 
   ../compute_gOOr.out < input 
   ../compute_gOHr.out < input 
   ../compute_gHHr.out < input 
   rm -f input
  
   mv gOO.dat gOO.pos_${i}.dat
   mv gOH.dat gOH.pos_${i}.dat
   mv gHH.dat gHH.pos_${i}.dat

   #
   # computing msds ....
   #
   echo $file > input 
   echo $nstep >> input 
   ../compute_msd.out < input 
   rm -f input

   mv msd.dat msd.pos_${i}.dat
  
   # remove the soft link
   rm -f $file
   echo Finished $file 
  done
  
  # computing mean g(r)s ....
  paste gOO.pos_0.dat gOO.pos_1.dat gOO.pos_2.dat gOO.pos_3.dat gOO.pos_4.dat gOO.pos_5.dat gOO.pos_6.dat gOO.pos_7.dat | awk '{printf "%15.8f %15.8f\n", $1,($2+$4+$6+$8+$10+$12+$14+$16)/8}' > gOO-mean.dat
  paste gOH.pos_0.dat gOH.pos_1.dat gOH.pos_2.dat gOH.pos_3.dat gOH.pos_4.dat gOH.pos_5.dat gOH.pos_6.dat gOH.pos_7.dat | awk '{printf "%15.8f %15.8f\n", $1,($2+$4+$6+$8+$10+$12+$14+$16)/8}' > gOH-mean.dat
  paste gHH.pos_0.dat gHH.pos_1.dat gHH.pos_2.dat gHH.pos_3.dat gHH.pos_4.dat gHH.pos_5.dat gHH.pos_6.dat gHH.pos_7.dat | awk '{printf "%15.8f %15.8f\n", $1,($2+$4+$6+$8+$10+$12+$14+$16)/8}' > gHH-mean.dat

  # computing mean msds ....
  paste msd.pos_0.dat msd.pos_1.dat msd.pos_2.dat msd.pos_3.dat msd.pos_4.dat msd.pos_5.dat msd.pos_6.dat msd.pos_7.dat | awk '{printf "%15.8f %15.8f\n", $1,($4+$8+$12+$16+$20+$24+$28+$32)/8}' > msd-mean.dat

  cd ../ 

done  
