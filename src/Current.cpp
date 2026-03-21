struct MVParameters
{
  // --- Valori limite del voltaggio ---
  double u_o        = 0.0;    // voltaggio a riposo (adimensionale)
  double u_u        = 1.61;   // voltaggio massimo upstroke

  // --- Soglie per le funzioni di Heaviside ---
  double theta_v    = 0.3;    // soglia per J_fi e gate v
  double theta_w    = 0.13;   // soglia per J_so, J_si, gate w
  double theta_v_m  = 0.1;    // soglia per tau_v- (quale ramo)
  double theta_o    = 0.005;  // soglia per tau_o e w_inf
 
  // --- Parametri per tau_v- (costante di tempo v in chiusura) ---
  double tau_v1m    = 80.0;   // tau_v - quando u < theta_vm
  double tau_v2m    = 1.4506; // tau_v- quando u > theta_vm
 
  // --- Parametri per tau_v+ (costante di tempo v in apertura) ---
  double tau_v_p    = 1.4506;
 
  // --- Parametri per tau_w- (costante di tempo w in chiusura) ---
  double tau_w1_m   = 70.0;   // valore minimo di tau_w-
  double tau_w2_m   = 8.0;    // valore massimo di tau_w-
  double k_w_m      = 200.0;  // slope della sigmoide per tau_w-
  double u_w_m      = 0.016;  // punto di mezzo della sigmoide
 
  // --- Parametri per tau_w+ (costante di tempo w in apertura) ---
  double tau_w_p    = 280.0;
 
  // --- Parametri per J_fi (fast inward - sodio) ---
  double tau_fi     = 0.078;
 
  // --- Parametri per J_so (slow outward - potassio) ---
  double tau_o1     = 410.0;  // tau_o quando u < theta_o
  double tau_o2     = 7.0;    // tau_o quando u > theta_o
  double tau_so1    = 91.0;   // tau_so minimo
  double tau_so2    = 0.8;    // tau_so massimo
  double k_so       = 2.1;    // slope della sigmoide per tau_so
  double u_so       = 0.6;    // punto di mezzo della sigmoide
 
  // --- Parametri per s (quarta variabile, morfologia AP) ---
  double tau_s1     = 2.7342; // tau_s quando u < theta_w
  double tau_s2     = 4.0;    // tau_s quando u > theta_w
  double k_s        = 2.0994; // slope della tanh per s_inf
  double u_s        = 0.9087; // punto di mezzo della tanh
 
  // --- Parametri per J_si (slow inward - calcio) ---
  double tau_si     = 3.3849;
 
  // --- Parametri per w_inf ---
  double tau_w_inf  = 0.01;   // usato nel calcolo di w_inf
  double w_inf_star = 0.5;    // valore di w_inf quando u > theta_o
 
  // --- Diffusione ---
  // D = 1.171 cm^2/s dal paper (Appendice A)
  // Qui in unità adimensionali del modello
  double D = 1.171e-4;        // adattato alle unità del problema
};