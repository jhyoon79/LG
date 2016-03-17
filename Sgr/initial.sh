#!/bin/bash

idl make_initial_sgr.idl
f77 frog_initial_sgr.f
./a.out 

