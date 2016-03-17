#!/bin/bash

rm -f snapshot/*
rm -f snapshot_subhalo/*
idl make_sgr.idl
f77 frog.f -o frog.out
./frog.out
idl sgr_plot.idl
