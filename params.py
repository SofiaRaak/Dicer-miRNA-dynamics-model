import numpy as np

#experimental data from figure 1, tsutsumi et al.
WT_data = np.array([0, 0.11144276160503169, 0.16566679779700877, 0.23905143587726366, 0.2954956726986665, 0.2946793863099961])
short_data = np.array([0, 0.0033684107002276975, 0.007599822974028003, 0.010019177812737812, 0.009603658536577298, 0.01242378048779691])
time = np.array([0, 5, 10, 20, 40, 60])
minutes = 60

Kd_wt = 25.4 #nM, experimental values
Kd_short = 147.7 #nM, experimental values

#K_d = k_off / k_on

WT_init = 1 #nM, from experimental setup
short_init = 1 #nM
dicer_init = 5 #nM
mirna_init = 0
WT_dicer_init = 0
short_dicer_init = 0

#optimised from CMA
ka1	= np.exp(-1.50569425)#np.exp(-4.91241215e+00)
ka2	= np.exp(1.8658517)#np.exp(-3.15282165e+00)
kc1	= np.exp(-4.39587544)#np.exp(-4.43307349e+00)
kc2	= np.exp(-3.9140815)#np.exp(-6.44167537e+03)
kd1	= np.exp(1.39394916)#np.exp(1.49475764e+00)
kd2 = np.exp(-3.3699480)#np.exp(-2.86279214e+00)

theta = [ka1, ka2, kc1, kc2, kd1, kd2]

init_pMiR1 = 1
init_pMiR2 = 1
init_pMiR1_dcr = 0
init_pMiR2_dcr = 0
init_pMiR1_dcr_star = 0
init_pMiR2_dcr_star = 0
init_dcr1 = 5
init_dcr2 = 5
init_MiR1 = 0
init_MiR2 = 0

init_values = [init_pMiR1, init_pMiR2, init_pMiR1_dcr, init_pMiR2_dcr, init_pMiR1_dcr_star,\
               init_pMiR2_dcr_star, init_dcr1, init_dcr2, init_MiR1, init_MiR2]

#Plotting colour scheme
colors = [('#332288', '#88CCEE'), ('#44AA99', '#117733'), ('#999933', '#DDCC77'), ('#CC6677', '#882255')]