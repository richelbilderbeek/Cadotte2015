#! /bin/bash

if [ ! -e EXP1_MLtree_codes.txt ]
then
  wget http://datadryad.org/bitstream/handle/10255/dryad.63823/EXP1_MLtree_codes.txt
fi

if [ ! -e planted.csv ]
then
  wget http://datadryad.org/bitstream/handle/10255/dryad.63820/planted.csv
fi

if [ ! -e anal_dat.csv ]
then

  wget http://datadryad.org/bitstream/handle/10255/dryad.63824/anal_dat.csv
fi

