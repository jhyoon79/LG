#!/bin/bash
idl make_pal5.idl
g77 frog.f -o frog.out
./frog.out
idl figure1.idl
#idl Pal5_plot.idl
