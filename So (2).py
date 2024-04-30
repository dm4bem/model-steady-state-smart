# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 14:38:11 2024

@author: matth
"""

import numpy as np

np.set_printoptions(precision=1)

# Data
# ====
# dimensions
L1, L2, c1, c2 = 4, 6, 3, 4  # m
H = 3        # m
H_vitre = 1   #m
L_vitre = 1  #m
W_mur_ext = 0.20,   # m
w_mur_int = 0.10,   # m

# B prend en compte w_int !!!!

# surfaces
S_A_mur_ext = (L1+(c1-L_vitre))*H + L_vitre*(H-H_vitre) + W_mur_ext*2*H

S_B_mur_ext = ((L2-L_vitre)+c1+2*W_mur_ext+w_mur_int)*H + L_vitre*(H-H_vitre)

S_C_mur_ext = (L1+L2+w_mur_int-2*L_vitre+2*c2+4*W_mur_ext)*H + 2*L_vitre*(H-H_vitre)


# thermo-physical propertites
λ = 1.7             # W/(m K) wall thermal conductivity
ρ, c = 1.2, 1000    # kg/m3, J/(kg K) density, specific heat air
hi, ho = 8, 25      # W/(m2 K) convection coefficients in, out

# short-wave solar radiation absorbed by each wall
E = 200             # W/m2

# outdoor temperature
To = 0              # °C

# ventilation rate (air-changes per hour)
ACH = 1             # volume/h

VA_dot =L1 * c1 * H * ACH / 3600  # volumetric air flow rate
VB_dot = L2 * c1 * H * ACH / 3600  # volumetric air flow rate
VC_dot = (L1+L2+w_mur_int) * c2 * H * ACH / 3600  # volumetric air flow rate
mA_dot = ρ * VA_dot               # mass air flow rate
mB_dot = ρ * VB_dot 
mC_dot = ρ * VC_dot 

nq, nθ = 20, 8  # number of flow-rates branches and of temperaure nodes

# Incidence matrix
# ================
A = np.zeros([nq, nθ])

# q0 ... q2 Convection extérieure
A[0, 0] = 1
A[1, 3] = 1
A[2, 5] = 1

# q3 ... q5 Conduction/convection extérieure/intérieure
A[3, 1], A[3, 0] = 1, -1
A[4, 3], A[4, 2] = 1, -1
A[5, 4], A[5, 5] = 1, -1


# q6 ... q8 vitre
A[6, 1], A[6, 6] = 1,-1
A[7, 2], A[7, 7] = 1, -1
A[8, 4], A[8, 8] = 1, -1


# q9 ... q11 Conduction/convection intérieure/intérieure
A[9, 1],A[9, 2] = -1,1
A[10, 2], A[10, 4] = 1, -1
A[11, 1],A[11, 4] = 1,-1

# q12 ... q14 Ventillation extérieure
A[12, 2] = 1
A[13, 3] = 1
A[14, 4] = 1

# q15, q16 Ventillation intérieure
A[15, 2],A[15, 4] = 1,-1
A[16, 1],A[15, 4] = 1,-1

# q17 ... q19 Régulateur de température
A[17, 1] = 1
A[18, 2] = 1
A[19, 4] = 1


# Conductance matrix
# ==================
G = np.zeros(A.shape[0])

# G0 ... G3 (cyan branches): outdoor convection
#L4 = 2 * l + 3 * L + 2 * w  	# length outdoor wall room 4
Lc = L1 + L2 + w_mur_int + c1 + c2 + w_mur_int_    # length outdoor wall room 4
#So = np.array([L, L + l, L + l, L4]) * H    # outdoor surfaces
So = np.array([L, L + l, L + l, L4]) * H
G[0:4] = ho * So

# G4 ... G7 (blue branches): conduction, indoor convection
G[4:8] = 1 / (w / λ + 1 / hi) * So

# G8 ... G12 (yellow branches): indoor walls
#    indoor convection, conduction, indoor convection
Si = np.array([l, l, L, L, L]) * H
G[8:13] = 1 / (1 / hi + w / λ + 1 / hi) * Si

# G13 ... G15 (green branches): advection by ventilation
G[13:16] = np.zeros(3)

# G16 ... G19 (red branches): gains of proportional controllers
G[16:20] = np.zeros(4)

# Vector of temperature sources
# =============================
b = np.zeros(A.shape[0])

b[0:4] = To         # cyan branches: outdoor temperature for walls
b[[13, 15]] = To    # green branches: outdoor temperature for ventilation
b[[16, 18]] = 20, 22    # red branches: setpoints room 1 & 3c

# Vector of flow-rate sources
# =============================
f = np.zeros(A.shape[1])

# Indexes of outputs
# ==================
indoor_air = [2, 3, 4, 6]   # indoor air temperature nodes
controller = range(16, 20)  # controller branches

# Question 1 : all rooms are controlled
# ==================
print(f"Maximum value of conductance: {max(G):.0f} W/K")

b[controller] = 20, 20, 22, 18  # °C setpoint temperature of the rooms
G[controller] = 1e9             # P-controller gain

θ = np.linalg.inv(A.T @ np.diag(G) @ A) @ (A.T @ np.diag(G) @ b + f)
q = np.diag(G) @ (-A @ θ + b)
print("1. All 4 rooms controlled")
print("θ:", θ[indoor_air], "°C")
print("q:", q[controller], "W")

# Question 2
# ===================
# Zone 2 & 4 free-running; solar rad; without ventilation
G[[17, 19]] = 0     # controller gains for room 2 & 4

# Solar radiation
exterior_wall = [0, 1, 5, 7]
f[exterior_wall] = E * So

θ = np.linalg.inv(A.T @ np.diag(G) @ A) @ (A.T @ np.diag(G) @ b + f)
q = np.diag(G) @ (-A @ θ + b)
print("2. 2 & 4 free-run w/o ventilation")
print("θ:", θ[indoor_air], "°C")
print("q:", q[controller], "W")

# Question 3
# ===================
# Zone 2 & 4 free-running; solar rad;
# Ventilation outdoor -> room 2 -> room 4 -> outdoor
ventilation = range(13, 16)
G[ventilation] = m_dot * c, m_dot * c, 0

θ = np.linalg.inv(A.T @ np.diag(G) @ A) @ (A.T @ np.diag(G) @ b + f)
q = np.diag(G) @ (-A @ θ + b)
print("3. 2 & 4 free-run, ventilation out -> 2 -> 4 -> out")
print("θ:", θ[indoor_air], "°C")
print("q:", q[controller], "W")

# Question 4
# ===================
# Zone 2 & 4 free-running; solar rad;
# Ventilation outdoor -> room 4 -> room 2 -> outdoor
G[ventilation] = 0, m_dot * c, m_dot * c

θ = np.linalg.inv(A.T @ np.diag(G) @ A) @ (A.T @ np.diag(G) @ b + f)
q = np.diag(G) @ (-A @ θ + b)
print("4. 2 & 4 free-run, ventilation out -> 4 -> 2 -> out")
print("θ:", θ[indoor_air], "°C")
print("q:", q[controller], "W")
