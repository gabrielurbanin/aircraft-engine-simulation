% Inputs do GE CF34-10
% Escolhas
M0 = 0.78;
H = 38000/3.281; %m %1 m = 3.28084 ft
[T0, a0, P0] = atmosisa(H);
% T0 = 288.15; %K
% P0 = 101325; %Pa
Tt4 = 925; %K
% Constantes
gamma_c = 1.4;
cp_c = 1004; %J/(kg.K)
gamma_t = 1.3;
cp_t = 1100; %J/(kg.K)
hpr = 42*10^6; %J/kg
pi_d_max = 0.97;
pi_b = 0.98;
pi_tH = 0.496748866560457;
pi_n = 0.92;
pi_fn = 0.99;
tau_tH = 0.857795681568010;
eta_f = 0.978446467818899;
eta_cL = 0.944819603365960;
eta_cH = 0.933486026234520;
eta_b = 0.7;
eta_mL = 0.98;
eta_mH = 0.98;
eta_tL = 0.953641682437914;
% Referência
M0_R = 0;
T0_R = 216.65; %K
P0_R = 2.064801539810828e+04; %Pa
tau_r_R = 1;
tau_lambda_R = 9.355607801887306;
pi_r_R = 1;
Tt4_R = 1850; %K
pi_d_R = 0.97;
pi_f_R = 1.70;
pi_cH_R = 8.29;
pi_cL_R = 2.02;
pi_tL_R = 0.506344666186594;
tau_f_R = 1.167310183961620;
tau_cH_R = pi_cH_R^((gamma_c - 1)/(gamma_c));
tau_cL_R = pi_cL_R^((gamma_c - 1)/(gamma_c));
tau_tL_R = 0.861401297789846;
alfa_R = 5;
M9_R = 1.634209811710537;
M19_R = 1.204746603708882;
m0_R = 242; %kg/s

% Equações
R_c = (gamma_c - 1)/gamma_c*cp_c; %J/(kg.K)
R_t = (gamma_t - 1)/gamma_t*cp_t; %J/(kg.K)
a0 = (gamma_c*R_c*T0)^(1/2); %m/s
V0 = a0*M0; %m/s
tau_r = 1 + (gamma_c - 1)/2*M0^2;
pi_r = tau_r^(gamma_c/(gamma_c - 1));
if M0 <= 1
    eta_r = 1;
else
    eta_r = 1 - 0.075*(M0 - 1)^1.35;
end
pi_d = pi_d_max*eta_r;
tau_lambda = cp_t*Tt4/(cp_c*T0);

teste = 10;
while teste > 0.0001
    tau_tL = tau_tL_R;
    tau_f = tau_f_R;
    tau_cL = tau_cL_R;
    pi_tL = pi_tL_R;
    pi_cL = pi_cL_R;
    tau_cH = 1 + Tt4/T0/(Tt4_R/T0_R)*(tau_f_R*tau_cL_R)/(tau_r*tau_cL)*(tau_cH_R - 1);
    pi_cH = (1 + eta_cH*(tau_cH - 1))^(gamma_c/(gamma_c - 1));
    pi_f = (1 + (tau_f - 1)*eta_f)^(gamma_c/(gamma_c - 1));
    Pt19_P0 = pi_r*pi_d*pi_f*pi_fn;
    if Pt19_P0 < ((gamma_c + 1)/2)^(gamma_c/(gamma_c - 1))
        Pt19_P19 = Pt19_P0;
    else
        Pt19_P19 = ((gamma_c + 1)/2)^(gamma_c/(gamma_c - 1));
    end
    M19 = (2/(gamma_c - 1)*(Pt19_P19^((gamma_c - 1)/gamma_c) - 1))^(1/2);
    Pt9_P0 = pi_r*pi_d*pi_f*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n;
    if Pt9_P0 < ((gamma_t + 1)/2)^(gamma_t/(gamma_t - 1))
        Pt9_P9 = Pt9_P0;
    else
        Pt9_P9 = ((gamma_t + 1)/2)^(gamma_t/(gamma_t - 1));
    end
    M9 = (2/(gamma_t - 1)*(Pt9_P9^((gamma_t - 1)/gamma_t) - 1))^(1/2);
    MFP_M19 = M19*(gamma_c/R_c)^(1/2)*(1 + (gamma_c - 1)/2*M19^2)^((gamma_c + 1)/(2*(1 - gamma_c)));
    MFP_M19_R = M19_R*(gamma_c/R_c)^(1/2)*(1 + (gamma_c - 1)/2*M19_R^2)^((gamma_c + 1)/(2*(1 - gamma_c)));
    alfa = alfa_R*pi_cL_R*pi_cH_R/pi_f_R/(pi_cL*pi_cH/pi_f)*((tau_lambda)/(tau_r*tau_f)/(tau_lambda_R/(tau_r_R*tau_f_R)))^(1/2)*MFP_M19/MFP_M19_R;
    tau_f = 1 + (tau_f_R - 1)*((1 - tau_tL)/(1 - tau_tL_R)*(tau_lambda/tau_r)/(tau_lambda_R/tau_r_R)*(tau_cL_R - 1 + alfa_R*(tau_f_R - 1))/(tau_cL_R - 1 + alfa*(tau_f_R - 1)));
    tau_cL = 1 + (tau_f - 1)*(tau_cL_R - 1)/(tau_f_R - 1);
    pi_cL = (1 + eta_cL*(tau_cL - 1))^(gamma_c/(gamma_c - 1));
    tau_tL = 1 - eta_tL*(1 - pi_tL^((gamma_t - 1)/gamma_t));
    MFP_M9 = M9*(gamma_t/R_t)^(1/2)*(1 + (gamma_t - 1)/2*M9^2)^((gamma_t + 1)/(2*(1 - gamma_t)));
    MFP_M9_R = M9_R*(gamma_t/R_t)^(1/2)*(1 + (gamma_t - 1)/2*M9_R^2)^((gamma_t + 1)/(2*(1 - gamma_t)));
    pi_tL = pi_tL_R*(tau_tL/tau_tL_R)^(1/2)*MFP_M9_R/MFP_M9;
    teste = abs(tau_tL - tau_tL_R);
    tau_tL_R = tau_tL;
    tau_f_R = tau_f;
    tau_cL_R = tau_cL;
    pi_tL_R = pi_tL;
end
tau_tL_R = 0.861401297789846;
tau_f_R = 1.167310183961620;
tau_cL_R = pi_cL_R^((gamma_c - 1)/(gamma_c));
pi_tL_R = 0.506344666186594;

m0 = m0_R*(1 + alfa)/(1 + alfa_R)*P0*pi_r*pi_d*pi_cL*pi_cH/(P0_R*pi_r_R*pi_d_R*pi_cL_R*pi_cH_R)*(Tt4_R/Tt4)^(1/2); %kg/s
f = (tau_lambda - tau_r*tau_cL*tau_cH)/(hpr*eta_b/(cp_c*T0) - tau_lambda); %kgFuel/kgAir
T9_T0 = tau_lambda*tau_tH*tau_tL/(Pt9_P9^((gamma_t - 1)/gamma_t))*cp_c/cp_t;
V9_a0 = M9*(gamma_t*R_t/(gamma_c*R_c)*T9_T0)^(1/2);
T19_T0 = tau_r*tau_f/(Pt19_P19^((gamma_c - 1)/gamma_c));
V19_a0 = M19*(T19_T0)^(1/2);
P19_P0 = Pt19_P0/(1 + (gamma_t - 1)/2*M19^2)^(gamma_t/(gamma_t - 1));
P9_P0 = Pt9_P0/(1 + (gamma_c - 1)/2*M9^2)^(gamma_c/(gamma_c - 1));
FF_m0 = alfa/(1 + alfa)*a0*(V19_a0 - M0 + T19_T0/V19_a0*(1 - 1/P19_P0)/gamma_c); %N/(kg/s)
FC_m0 = 1/(1 + alfa)*a0*((1 + f)*V9_a0 - M0 + (1 + f)*R_t*T9_T0/(R_c*V9_a0)*(1 - 1/P9_P0)/gamma_c); %N/(kg/s)
F_m0 = FF_m0  + FC_m0; %N/(kg/s)
S = f/((1 + alfa)*F_m0); %(kgFuel/s)/N
F = m0*F_m0; %N
N_NR_fan = (T0*tau_r/(T0_R*tau_r_R)*(pi_f^((gamma_c - 1)/gamma_c) - 1)/(pi_f_R^((gamma_c - 1)/gamma_c) - 1))^(1/2);
N_NR_H = (T0*tau_r*tau_cL/(T0_R*tau_r_R*tau_cL_R)*(pi_cH^((gamma_c - 1)/gamma_c) - 1)/(pi_cH_R^((gamma_c - 1)/gamma_c) - 1))^(1/2);
eta_T = a0^2*((1 + f)*V9_a0^2 + alfa*(V19_a0^2)- (1 + alfa)*M0^2)/(2*f*hpr);
eta_P = 2*V0*(1 + alfa)*F_m0/(a0^2*((1 + f)*V9_a0^2 + alfa*V19_a0^2 - (1 + alfa)*M0^2));
eta_Total = eta_P*eta_T;
mf = S*F; %kgFuel/s vazao de combustivel
AF = 1/f; %kgAir/kgFuel
mC = m0*1/(1 + alfa); %kgAir/s vazao de ar pelo Core
mF = m0*alfa/(1 + alfa); %kgAir/s vazao de ar pelo Fan
