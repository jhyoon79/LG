#!/bin/bash

idl make_encounter_xy.idl
f77 frog_xy.f -o frog_xy.out
./frog_xy.out
idl various_encounter.idl
