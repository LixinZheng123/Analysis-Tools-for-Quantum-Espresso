#!/bin/bash

#OH- Names
OH_prefix1=PBE0_OH-_vdw
OH_prefix2=PBE_OH-_vdw
OH_prefix3=PBE_OH-_novdw_flexible
OH_Dir=OH-_Data
#H+ Names
H_prefix1=PBE0_H+_vdw
H_prefix2=PBE_H+_vdw
H_Dir=H+_Data

#Read from screen which group of data should we analyse
echo "Do you want to clean (1) OH- (2) H+ or (3)all? (Reply a number)"
read n

if [[ $n = 1 || $n = 3 ]]; then
  mv *$OH_prefix1.* $OH_Dir
  mv *$OH_prefix2.* $OH_Dir
  mv *$OH_prefix3.* $OH_Dir
fi

if [[ $n = 2 || $n = 3 ]]; then
  mv *$H_prefix1.* $H_Dir
  mv *$H_prefix2.* $H_Dir
fi
