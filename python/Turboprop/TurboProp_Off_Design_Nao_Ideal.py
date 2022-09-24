from ambiance import Atmosphere

# Inputs do PT6A-68C
# Escolhas
M0 = 0.03
H = 1000 #m
# [T0, a0, P0] = atmosisa(H)
atmosphere = Atmosphere(H)
T0 = atmosphere.temperature
P0 = atmosphere.pressure
a0 = atmosphere.speed_of_sound

T0 = 288.15
P0 = 101325
Tt4 = 1350 #k
# Constantes
gamma_c = 1.4
cp_c = 1004 #J/(kg.K)
gamma_t = 1.3
cp_t = 1108 #J/(kg.K)
hpr = 42*10**6 #J/kg
pi_d_max = 0.99
pi_b = 1
pi_tH = 0.499
pi_n = 1
tau_tH = 0.823
eta_c = 0.95
eta_b = 0.98
eta_tL = 0.98
eta_mL = 0.98
eta_g = 0.9409
eta_prop = 0.24
# Referências
M0_R = 0.03
T0_R = 288.15 #K
P0_R = 101325 #Pa
tau_r_R = 1.00018
tau_lambda_R = 5.17
pi_r_R = 1.00063014176276
Tt4_R = 1350 #k
pi_d_R = 0.98
pi_c_R = 8.70
pi_tL_R = 0.421
tau_c_R = 1.92
tau_tL_R = 0.785
M9_R = 0.18
m0_R = 4.78 #kg/s

# Equações
R_c = (gamma_c - 1)/gamma_c*cp_c #J/(kg.K)
R_t = (gamma_t - 1)/gamma_t*cp_t #J/(kg.K)
a0 = (gamma_c*R_c*T0)**(1/2) #m/s
V0 = a0*M0 #m/s
tau_r = 1 + (gamma_c - 1)/2*M0**2
pi_r = tau_r**(gamma_c/(gamma_c - 1))
eta_r = 1
pi_d = pi_d_max*eta_r
tau_c = 1 + Tt4/T0/(Tt4_R/T0_R)*(tau_r_R)/(tau_r)*(tau_c_R - 1)
pi_c = (1 + eta_c*(tau_c - 1))**(gamma_c/(gamma_c - 1))
tau_lambda = cp_t*Tt4/(cp_c*T0)
f = (tau_lambda - tau_r*tau_c)/(hpr*eta_b/(cp_c*T0) - tau_lambda) #kgFuel/kgAir
m0 = m0_R*P0*pi_r*pi_d*pi_c/(P0_R*pi_r_R*pi_d_R*pi_c_R)*(Tt4_R/Tt4)**(1/2) #kg/s

# teste = 10
# while teste > 0.0001:
#     pi_tL = pi_tL_R
#     tau_tL = 1 - eta_tL*(1 - pi_tL**((gamma_t - 1)/gamma_t))
#     Pt9_P0 = pi_r*pi_d*pi_c*pi_b*pi_tH*pi_tL*pi_n

#     # import pdb; pdb.set_trace()
#     if Pt9_P0 >= ((gamma_t + 1)/2)**(gamma_t/(gamma_t - 1)):
#         M9 = 1
#         Pt9_P9 = ((gamma_t + 1)/2)**(gamma_t/(gamma_t - 1))
#         P0_P9 = Pt9_P9/Pt9_P0
#     else:
#         P0_P9 = 1
#         Pt9_P9 = Pt9_P0
#         M9 = (2/(gamma_t - 1)*(Pt9_P0**((gamma_t - 1)/gamma_t) - 1))**(1/2)

#     M9
#     MFP_M9 = M9*(gamma_t/R_t)**(1/2)*(1 + (gamma_t - 1)/2*M9**2)**((gamma_t + 1)/(2*(1 - gamma_t)))
#     MFP_M9_R = M9_R*(gamma_t/R_t)**(1/2)*(1 + (gamma_t - 1)/2*M9_R**2)**((gamma_t + 1)/(2*(1 - gamma_t)))
#     pi_tL = pi_tL_R*(tau_tL/tau_tL_R)**(1/2)*MFP_M9_R/MFP_M9
#     teste = abs(pi_tL - pi_tL_R)
#     pi_tL_R = pi_tL

pi_tL_R = 0.421


T9_T0 = tau_lambda*tau_tH*tau_tL/((Pt9_P9)**((gamma_t - 1)/gamma_t))
V9_a0 = M9*(gamma_t*R_t/(gamma_c*R_c)*T9_T0)**(1/2)
C_c = (gamma_c - 1)*M0*((1 + f)*V9_a0 - M0 + (1 + f)*R_t/R_c*T9_T0/V9_a0*(1 - P0_P9)/gamma_c)
C_prop = eta_prop*eta_g*eta_mL*(1 + f)*tau_lambda*tau_tH*(1 - tau_tL)
C_Total = C_prop + C_c
F_m0 = C_Total*cp_c*T0/V0
S = f/F_m0 #kg/s/N
W_m0 = C_Total*cp_c*T0 #W/(kg/s)
S_P = f/(C_Total*cp_c*T0) #kgFuel/s/W vazão mássica por potência
F = m0*F_m0 #N
W = m0*W_m0 #W
eta_P = C_Total/(C_prop/eta_prop + ((gamma_c - 1)/2)*((1 + f)*V9_a0**2 - M0**2))
eta_T = C_Total*cp_c*T0/(f*hpr)
eta_Total = eta_P*eta_T
N_NR_core = (T0*tau_r/(T0_R*tau_r_R)*(tau_c - 1)/(tau_c_R - 1))**(1/2)
N_NR_power = (Tt4/(Tt4_R)*(1 - tau_tL)/(1 - tau_tL_R))**(1/2)


