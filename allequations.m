
C_o = Co;
C_f = 15 *10^-3;
Q_f = 20;
b_fm = bfm;
L_E = 10 *10^3;
R = 2 *10^-3/365/24/60/60;
T_T = 12 *60*60;
H = 1.4/2;
b_f = x(1);
d_f = x(2);
d_m = x(3);

H_w = .2413*(tanh((.493*(g*((d_f+max(0,d_f-2*H))/2)/v_w^2)^.75))*tanh((3.13*10^-3*(g*(2*b_f)/v_w^2)^.57)/tanh((.493*(g*((d_f+max(0,d_f-2*H))/2)/v_w^2)^.75))))^.87*v_w^2/g;
T_w = 7.518*(tanh((.331*(g*((d_f+max(0,d_f-2*H))/2)/v_w^2)^1.01))*tanh((5.215*10^-4*(g*(2*b_f)/v_w^2)^.73)/tanh((.331*(g*((d_f+max(0,d_f-2*H))/2)/v_w^2)^1.01))))^.37*v_w/g;
k_w = 2*pi/T_w/sqrt(g*((d_f+max(0,d_f-2*H))/2));
f_w = 0.4*(H_w/k_0/sinh(k_w*((d_f+max(0,d_f-2*H))/2)))^-.75;
u_w = pi*H_w/T_w/sinh(k_w*((d_f+max(0,d_f-2*H))/2));
tau = rho_w*f_w*u_w^2/2;

bed_erosion = max(0,(1/2-1/pi*asin((H-d_f)/H))*E_0*(tau-tau_c)/tau_c*b_f*L_E);
ocean_in = (max(((d_f*b_f+d_m*(b_fm-b_f))*L_E)/T_T-Q_f,0))*C_o;
river_in = Q_f*C_f;
TF_deposition = min((1/2-1/pi*asin((H-d_f)/H))*((R-(k_B*(B_max*(1-(0.5*(H-d_m)/H)/(-0.5*(H-d_m)/H+1)))))*T_T*rho_s/d_m)*b_f*omega_s*L_E,((R-(k_B*(B_max*(1-(0.5*(H-d_m)/H)/(-0.5*(H-d_m)/H+1)))))*T_T*rho_s/d_m)*b_f*d_f*L_E/T_T);
M_deposition = ((R-(k_B*(B_max*(1-(0.5*(H-d_m)/H)/(-0.5*(H-d_m)/H+1)))))*T_T*rho_s/d_m)*(b_fm-b_f)*d_m*L_E/T_T;
export = ((R-(k_B*(B_max*(1-(0.5*(H-d_m)/H)/(-0.5*(H-d_m)/H+1)))))*T_T*rho_s/d_m)*b_f*min(d_f,2*H)*L_E/T_T;

G(1) = max(0,(1/2-1/pi*asin((H-d_f)/H))*E_0/rho_s*(tau-tau_c)/tau_c)-min((1/2-1/pi*asin((H-d_f)/H))*((R-(k_B*(B_max*(1-(0.5*(H-d_m)/H)/(-0.5*(H-d_m)/H+1)))))*T_T*rho_s/d_m)*omega_s/rho_s,((R-(k_B*(B_max*(1-(0.5*(H-d_m)/H)/(-0.5*(H-d_m)/H+1)))))*T_T*rho_s/d_m)*d_f/T_T/rho_s)+R;
G(2) = (bed_erosion + ocean_in + river_in - TF_deposition - M_deposition - export)/rho_s/L_E/b_fm;
F = (k_e*(gamma*(pi/k_w/T_w*(1+2*k_w*((d_f+max(0,d_f-2*H))/2)/sinh(2*k_w*((d_f+max(0,d_f-2*H))/2))))*H_w^2/16))-(k_a*omega_s*((R-(k_B*(B_max*(1-(0.5*(H-d_m)/H)/(-0.5*(H-d_m)/H+1)))))*T_T*rho_s/d_m)/rho_s);
