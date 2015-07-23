#!/usr/bin/env bash

for i in $(seq 1 20);
do

rm -f orbit_$i*.txt
rm -f perturbed/orbit_$i*.txt
rm -f unperturbed/orbit_$i*.txt

done

rm -f plots_unpert_vs_pert/*
rm -f orbit*.txt
rm -f perturbed/*
rm -f unperturbed/*
rm -f jP_*
rm -f *.png
