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
H_vitre = 1   # m
L_vitre = 1  # m
W_mur_ext = 0.30   # m
w_mur_int_beton = 0.20   # m
w_mur_int_insulation = 0.1 # m
w_mur_int = 0.3 # m
w_vitre = 0.02 # m

# Pièce B prend en compte w_int !!!!

# surfaces
S_A_mur_ext = (L1+(c1-L_vitre))*H + L_vitre*(H-H_vitre) + W_mur_ext*2*H
S_A_mur_int = S_A_mur_ext - W_mur_ext*2*H
S_B_mur_ext = ((L2-L_vitre)+c1+2*W_mur_ext+w_mur_int)*H + L_vitre*(H-H_vitre)
S_B_mur_int = S_B_mur_ext - 2*W_mur_ext+w_mur_int)*H
S_C_mur_ext = (L1+L2+w_mur_int-2*L_vitre+2*c2+4*W_mur_ext)*H + 2*L_vitre*(H-H_vitre)
S_C_mur_int = S_C_mur_ext - 4*W_mur_ext*H
S_fenetre = H_vitre*L_vitre
S_AB = c1*H
S_AC = L1*H
S_BC = L2*H


# thermo-physical propertites
λ_concrete = 1.4             # W/(m K) concrete thermal conductivity
λ_insulation = 0.04          # W/(m K) insulation thermal conductivity
λ_window = 1.4               # W/(m K) insulation thermal conductivity
ρ, c = 1.2, 1000    # kg/m3, J/(kg K) density, specific heat air
hi, ho = 8, 25      # W/(m2 K) convection coefficients in, out

# short-wave solar radiation absorbed by each wall
E = 200             # W/m2

# outdoor temperature
To = 0              # °C

# ventilation rate (air-changes per hour)
ACH = 1             # volume/h

VA_dot = L1 * c1 * H * ACH / 3600  # volumetric air flow rate
VB_dot = L2 * c1 * H * ACH / 3600  # volumetric air flow rate
VC_dot = (L1+L2+w_mur_int) * c2 * H * ACH / 3600  # volumetric air flow rate
mA_dot = ρ * VA_dot               # mass air flow rate
mB_dot = ρ * VB_dot               # mass air flow rate
mC_dot = ρ * VC_dot               # mass air flow rate



# radiative properties
ε_wLW = 0.85    # long wave emmisivity: wall surface (concrete)
ε_gLW = 0.90    # long wave emmisivity: glass pyrex
α_wSW = 0.25    # short wave absortivity: white smooth surface
α_gSW = 0.38    # short wave absortivity: reflective blue glass
τ_gSW = 0.30    # short wave transmitance: reflective blue glass
σ = 5.67e-8     # W/(m²⋅K⁴) Stefan-Bolzmann constant


nq, nθ = 23, 8  # number of flow-rates branches and of temperaure nodes

# Incidence matrix
# ================
A = np.zeros([nq, nθ])

# q0 ... q5 Convection extérieure
A[0, 3] = 1
A[1, 7] = 1
A[2, 11] = 1
A[3, 5] = 1
A[4, 9] = 1
A[5, 12] = 1

# q6 q10 q14  Conduction murs extérieurs
A[6, 4], A[6, 3] = 1, -1
A[10, 7], A[10, 8] = -1, 1
A[14, 13], A[14, 11] = 1, -1

# q7 q11 q15 Conduction vitres extérieures
A[7, 6], A[7, 5] = 1, -1
A[11, 10], A[11, 9] = 1, -1
A[15, 12], A[15,14] = -1, 1

# q8 q12 q16 Convection murs intérieurs
A[8, 0], A[8, 4] = 1, -1
A[12, 8], A[12, 1] = 1, -1
A[16, 2], A[16,13] = 1, -1

# q9 q13 q17 Convection vitres intérieures
A[9, 0], A[9, 6] = 1, -1
A[13, 1], A[13, 10] = 1, -1
A[17, 2], A[17,14] = 1, -1

# q18 ... q20 radiation entre vitres et murs
A[18, 4], A[18, 6] = 1, -1
A[19, 8], A[19, 10] = 1, -1
A[20, 13], A[20,14] = 1, -1


# q21 ... q23 Conduction/convection intérieure/intérieure
A[21, 1],A[21, 0] = 1,-1
A[22, 1], A[22, 2] = 1, -1
A[23, 0],A[23, 2] = 1,-1

# q24 ... q26 Ventillation extérieure
A[24, 0] = 1
A[25, 1] = 1
A[26, 2] = 1

# q27, q28 Ventillation intérieure
A[27, 0],A[27, 2] = 1,-1
A[28, 1],A[28, 2] = 1,-1

# q29 ... q31 Régulateur de température
A[29, 0] = 1
A[30, 1] = 1
A[31, 2] = 1


# Conductance matrix
# ==================
G = np.zeros(A.shape[0])

# G0 ... G2 : outdoor convection
G[0] = ho * S_A_mur_ext
G[1] = ho * S_B_mur_ext
G[2] = ho * S_C_mur_ext

# G3 ... G5 : conduction, indoor convection for the wall
G[3] = 1 / (w_mur_int_beton / λ_concrete + w_mur_int_insulation / λ_insulation + 1 / hi) * S_A_mur_ext
G[4] = 1 / (w_mur_int_beton / λ_concrete + w_mur_int_insulation / λ_insulation + 1 / hi) * S_B_mur_ext
G[5] = 1 / (w_mur_int_beton / λ_concrete + w_mur_int_insulation / λ_insulation + 1 / hi) * S_C_mur_ext


# G6 ... G8 : outdoor convection, conduction, indoor convection for the window
G[6:7] = 1 / (w_vitre / λ_window + 1/ho + 1 / hi) * S_fenetre
G[8] = 1 / (w_vitre / λ_window + 1/ho + 1 / hi) * 2*S_fenetre



# G9 ... G11 : solar radiation

# long wave radiation
Tm = 20 + 273   # K, mean temp for radiative exchange

F_A = S_fenetre / S_A_mur_int
GLW_mur_A = 4 * σ * Tm**3 * ε_wLW / (1 - ε_wLW) * S_A_mur_int
GLW_mur_vitre_A = 4 * σ * Tm**3 * F_A * S_A_mur_int
GLW_vitre_A = 4 * σ * Tm**3 * ε_gLW / (1 - ε_gLW) * S_fenetre
G[9] = 1 / (1 / GLW_mur_A + 1 / GLW_mur_vitre_A + 1 / GLW_vitre_A)

F_B = S_fenetre / S_B_mur_int
GLW_mur_B = 4 * σ * Tm**3 * ε_wLW / (1 - ε_wLW) * S_B_mur_int
GLW_mur_vitre_B = 4 * σ * Tm**3 * F_B * S_B_mur_int
GLW_vitre_B = 4 * σ * Tm**3 * ε_gLW / (1 - ε_gLW) * S_fenetre
G[10] = 1 / (1 / GLW_mur_B + 1 / GLW_mur_vitre_B + 1 / GLW_vitre_B)

F_C = 2*S_fenetre / S_C_mur_int
GLW_mur_C = 4 * σ * Tm**3 * ε_wLW / (1 - ε_wLW) * S_C_mur_int
GLW_mur_vitre_C = 4 * σ * Tm**3 * F_C * S_C_mur_int
GLW_vitre_C = 4 * σ * Tm**3 * ε_gLW / (1 - ε_gLW) * 2*S_fenetre
G[11] = 1 / (1 / GLW_mur_C + 1 / GLW_mur_vitre_C + 1 / GLW_vitre_C)

# G12 ... G14 : internal convection, conduction, internal conduction
G[12] = 1 / (w_mur_int_beton / λ_concrete + 2 / hi) * S_AB
G[13] = 1 / (w_mur_int_beton / λ_concrete + 2 / hi) * S_BC
G[14] = 1 / (w_mur_int_beton / λ_concrete + 2 / hi) * S_AC

# G15 ... G19 : Ventilation
G[15] = mA_dot*c
G[16] = mB_dot*c
G[17] = mC_dot*c
G[18]
G[19]

# G20 ... G22 : gains of proportional controllers
G[20:22] = np.zeros(3)

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
