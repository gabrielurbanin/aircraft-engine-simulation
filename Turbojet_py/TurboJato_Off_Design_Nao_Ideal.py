from ambiance import Atmosphere
from TurboJato_Nao_Ideal_function import TurboJato_Nao_Ideal_function
# Tipo do motor
# 1: Rolls-Royce Nene
# 2: VT80
tipo = 2

if tipo == 1:
    #  Inputs do Rolls-Royce Nene
    #  Escolhas
    H = 1000 # m
    # [T0, a0, P0] = atmosisa(H)
    M0 = 1.5
    T0 = 229.8 # K
    P0 = 30.8*10**3 # Pa
    Tt4 = 1670 # K
    P0_P9 = 0.955
    # Constantes
    gamma_c = 1.4
    cp_c = 1004 # J/(kg.K)
    gamma_t = 1.3
    cp_t = 1239 # J/(kg.K)
    hpr = 42.8*10**6 # J/kg
    pi_d_max = 0.95
    pi_b = 0.94
    pi_t = 0.3746
    pi_n = 0.96
    tau_t = 0.8155
    eta_c = 0.8641
    eta_b = 0.98
    eta_m = 0.99
    # Condições de referência
    [M0_R,T0_R,P0_R,tau_r_R,pi_r_R,Tt4_R,pi_d_R,pi_c_R,tau_c_R,Pt9_P9_R] = TurboJato_Nao_Ideal_function(tipo)
    #  Inputs extras do Rolls-Royce Nene
    m0_R = 50 # kg/s (?)
    Tt2_R = T0_R*tau_r_R

if tipo == 2:
    #  #  Inputs do VT80
    #  Escolhas
    H = 1000 # m
    # [T0, a0, P0] = atmosisa(H)
    atmosphere = Atmosphere(H)
    T0 = atmosphere.temperature
    P0 = atmosphere.pressure
    a0 = atmosphere.speed_of_sound
    M0 = 0.5
#      T0 = 281.65 # K
#      P0 = 89.87*10**3 # Pa
    Tt4 = 914.95 # K
    P0_P9 = 0.8388
    # Constantes
    gamma_c = 1.4
    cp_c = 1004 # J/(kg.K)
    gamma_t = 1.3
    cp_t = 1239 # J/(kg.K)
    hpr = 42.8*10**6 # J/kg
    pi_d_max = 0.95
    pi_b = 0.95
    pi_t = 0.6550
    pi_n = 0.9872
    tau_t = 0.9230
    eta_c = 0.8879
    eta_b = 0.99
    eta_m = 0.99
    # Condições de referência
    # [M0_R,T0_R,P0_R,tau_r_R,pi_r_R,Tt4_R,pi_d_R,pi_c_R,tau_c_R,Pt9_P9_R] = TurboJato_Nao_Ideal_function(tipo)
    #  Inputs extras do VT80
    M0_R = 0.03
    T0_R = 298.15 # K
    P0_R = 101325 # Pa
    tau_r_R = 1 + (gamma_c - 1)/2*M0_R**2
    Tt4_R = 914.95 # K
    if M0_R <= 1:
        eta_r = 1
    else:
        eta_r = 1 - 0.075*(M0_R - 1)**1.35

    pi_d_R = pi_d_max*eta_r
    m0_R = 0.2 # kg/s (?)
    Tt2_R = T0_R*tau_r_R
    pi_c_R = 2.25
    tau_c_R = pi_c_R**((gamma_c-1)/gamma_c)
    pi_r_R = 0.98
    Pt9_P9_R = 1.1


#  Equations
R_c = (gamma_c - 1)/gamma_c*cp_c # J/(kg.K)
R_t = (gamma_t - 1)/gamma_t*cp_t # J/(kg.K)
a0 = (gamma_c*R_c*T0)**(1/2) # m/s
V0 = a0*M0
tau_r = 1 + (gamma_c - 1)/2*M0**2
pi_r = tau_r**(gamma_c/(gamma_c - 1))
if M0 <= 1:
    eta_r = 1
else:
    eta_r = 1 - 0.075*(M0 - 1)**1.35

pi_d = pi_d_max*eta_r
Tt2 = T0*tau_r
tau_c = 1 + (tau_c_R - 1)*Tt4/Tt2/(Tt4_R/Tt2_R)
pi_c = (1 + eta_c*(tau_c - 1))**(gamma_c/(gamma_c - 1))
tau_lambda = cp_t*Tt4/(cp_c*T0)
f = (tau_lambda - tau_r*tau_c)/(hpr*eta_b/(cp_c*T0) - tau_lambda) # kgFuel/kgAir
m0 = m0_R*P0*pi_r*pi_d*pi_c/(P0_R*pi_r_R*pi_d_R*pi_c_R)*(Tt4_R/Tt4)**(1/2) # kg/s
Pt9_P9 = P0_P9*pi_r*pi_d*pi_c*pi_b*pi_t*pi_n
M9 = (2/(gamma_t - 1)*(Pt9_P9**((gamma_t - 1)/gamma_t) - 1))**(1/2)
T9_T0 = tau_lambda*tau_t/(Pt9_P9**((gamma_t - 1)/gamma_t))*cp_c/cp_t
V9_a0 = M9*(gamma_t*R_t/(gamma_c*R_c)*T9_T0)**(1/2)
F_m0 = a0*((1 + f)*V9_a0 - M0 + (1 + f)*R_t*T9_T0/(R_c*V9_a0)*(1 - P0_P9)/gamma_c) # N/(kg/s)
F = F_m0*m0 # N
S = f/F_m0 # (kgFuel/s)/N
eta_T = a0**2*((1 + f)*V9_a0**2 - M0**2)/(2*f*hpr)
eta_P = 2*V0*F_m0/(a0**2*((1 + f)*V9_a0**2 - M0**2))
eta_Total = eta_P*eta_T
N_NR = (T0*tau_r/(T0_R*tau_r_R)*(pi_c**((gamma_c - 1)/gamma_c) - 1)/(pi_c_R**((gamma_c - 1)/gamma_c) - 1))**(1/2)
A9_A9R = (Pt9_P9/Pt9_P9_R)**((gamma_t + 1)/(2*gamma_t))*((Pt9_P9_R**((gamma_t - 1)/gamma_t) - 1)/(Pt9_P9**((gamma_t - 1)/gamma_t) - 1))**(1/2)
mc2_mc2_R = pi_c/pi_c_R*((Tt4_R/Tt2_R)/(Tt4/Tt2))**(1/2) # vazão mássica corrigida no compressor

#  Outputs extras
V9 = V9_a0*a0 # m/s
AF = 1/f  # kgAir/kgFuel
Pt4 = P0*pi_r*pi_d*pi_c*pi_b # Pa
Pt9 = P0*pi_r*pi_d*pi_c*pi_b*pi_t # Pa
T9 = T0*T9_T0 # K



