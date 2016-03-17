#!/bin/sh

alias idl="/Applications/exelis/idl82/bin/idl"

idl make_initial_Pal5.idl
gfortran frog_initial_Pal5.f
./a.out


gfortran frog_central.f
./a.out

idl make_pal5.idl
gfortran frog.f -o frog.out
./frog.out
idl figure1.idl
