#!/bin/bash

idl make_evol_encounter_xy.idl
f77 frog_xy_evol.f -o frog_xy_evol.out
./frog_xy_evol.out
