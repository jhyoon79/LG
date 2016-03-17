#!/bin/bash
idl make_encounter_xy.idl
f77 frog_xy.f
./a.out
idl one_encounter_snap_movie.idl
