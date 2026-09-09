#! /bin/bash
  python3 clean.py
  mpirun -np 4 eqdyna
  python3 plotRuptureDynamics
