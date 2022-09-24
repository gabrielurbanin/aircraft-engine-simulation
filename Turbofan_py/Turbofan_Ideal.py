# Inputs

M0 = 0.1
T0 = 288.7 # Kelvin (K)
gamma = 1.4
cp = 1020; # J/(kg.K)
hpr = 120 * 10**6; # J/kg
Tt4 = 1700; # K
pi_c = 15
pi_f = 1.7
alfa = 5

# Equations

R = (gamma - 1)/gamma * cp
a0 = (gamma * R * T0)**(1/2)
tau_r = 1 + (gamma - 1)/2 * M0**2
tau_lambda = Tt4/T0
tau_c = pi_c**((gamma - 1) / gamma)
tau_f = pi_f**((gamma - 1) / gamma)
tau_t = 1 - tau_r / tau_lambda * (tau_c - 1)
V9 = (2/(gamma - 1) * (tau_lambda - tau_r * (tau_c - 1 + alfa * (tau_f - 1)) - tau_lambda/(tau_r * tau_c)))**(1/2) * a0
V19 = (2/(gamma - 1) * (tau_r * tau_f - 1))**(1/2) * a0
F_m0 = a0 * 1/(1 + alfa) * (V9/a0 - M0 + alfa * (V19/a0 - M0))
f = cp * T0/hpr * (tau_lambda - tau_r * tau_c)
S = f/((1 + alfa) * F_m0)
eta_T = 1 - 1/(tau_r * tau_c)
eta_P = 2 * M0 * (V9/a0 - M0 +alfa * (V19/a0 - M0))/(V9**2/(a0**2) - M0**2 + alfa * (V19**2/(a0**2) - M0**2))
eta_Total = eta_P * eta_T
FR = (V9/a0 - M0)/(V19/a0 - M0)

# F_m0/1000 %N/(kg/s)
# f
# S %(kg/s)/N
# eta_T
# eta_P
# eta_Total
# FR

