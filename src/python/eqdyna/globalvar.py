"""globalvar.f90 -- the module-level constants the solver shares.

ONLY constants live here, exactly as in src/globalvar.f90. The `fric` array
is 0-indexed in Python and 1-indexed in Fortran, so every slot below is
`FRIC_SLOT_* - 1`. Naming them is not cosmetic: the ports this replaces
addressed these columns as bare integers (`fric[:, 70]`, `fric[:, 5]`), and
a column read from the wrong slot is a silent wrong answer, not a crash.
"""

# --- fric(:) slot map, globalvar.f90:22-72 (Fortran slot minus 1) ----------
SW_FS = 0              # slip-weakening static friction coefficient
SW_FD = 1              # slip-weakening dynamic friction coefficient
SW_D0 = 2              # slip-weakening critical slip distance D0, m
COHESION = 3           # fault cohesion, Pa
TW_T0 = 4              # time-weakening rupture-time offset t0, s
INIT_NORM = 6          # initial effective normal stress, Pa
INIT_STRIKE_SHEAR = 7  # initial strike-direction shear stress, Pa
RSF_A = 8
RSF_B = 9
RSF_DC = 10
RSF_V0 = 11
RSF_R0 = 12
RSF_FW = 13
RSF_VW = 14
TP_A_HY = 15            # thermal pressurization hydraulic diffusivity
TP_A_TH = 16            # thermal pressurization thermal diffusivity
TP_ROUC = 17            # rho*c, heat capacity per volume
TP_LAMBDA = 18          # pore-pressure/temperature coupling
STATE = 19             # RSF state variable theta
THETA_PC = 22          # normal-stress-evolution state variable (Shi & Day 2013)
THETA_PC_DOT = 23
VINI_N = 24            # background/creep slip-rate offset, normal
VINI_X = 25            # ... strike (x)
VINI_Z = 26            # ... dip (z)
VEL_MASTER_X = 30      # 30,31,32 = master x,y,z
VEL_SLAVE_X = 33       # 33,34,35 = slave  x,y,z
TP_H = 39
TP_TINI = 40
TP_PINI = 41
CREEP_VMIN = 45        # creeping/initial slip-rate lower bound, m/s
PEAK_SLIPRATE = 46
SHEAR_MAG = 47
INIT_DIP_SHEAR = 48
TP_NORM_TP = 50        # TP-derived pore-pressure contribution to normal traction
TP_TEMP = 51
SLIP_STRIKE = 70
SLIP_DIP = 71
SLIP_NORM = 72
SLIPRATE_STRIKE = 73
SLIPRATE_DIP = 74
SLIPRATE_MAX = 75
CUM_SLIP = 76
TRACT_NORM = 77        # 77,78,79 = n,s,d -- written as a contiguous triple
TRACT_STRIKE = 78
TRACT_DIP = 79
NUC_DTAU0 = 80

# --- nucleation constants, globalvar.f90:79-81 ----------------------------
# TRUNCATED AS WRITTEN IN THE FORTRAN. Do not substitute "more correct"
# values: the port must reproduce the reference's arithmetic, not improve it.
NUC_VS_FIXED = 3464.0    # fixed shear-wave speed, m/s
NUC_TAPER_COEF = 0.081   # rupture-time taper coefficient
NUC_VR_TO_VS = 0.7       # nucleation rupture-speed-to-Vs ratio

# --- normal-stress caps, globalvar.f90:153-154 ---------------------------
# NOT case parameters. min_norm > max_norm, so faulting.f90's if/elseif
# branches are mutually exclusive and values strictly between max_norm and
# min_norm pass through untouched -- that gap is deliberate, reproduce it.
MAX_NORM = -40.0e6
MIN_NORM = -10.0e6

# eqdyna3d.f90:140 -- fnft's "never ruptured" sentinel, and the >5000 test
# storeRuptureTime uses against it (faulting.f90:312).
FNFT_SENTINEL = 99999.0
FNFT_UNRUPTURED_ABOVE = 5000.0
