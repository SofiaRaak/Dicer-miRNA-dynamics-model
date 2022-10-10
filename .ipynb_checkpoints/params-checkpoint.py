import numpy as np

Kd_wt = 25.4 #nM, experimental values

ka1 = 1.0285752
ka_1 = Kd_wt * ka1
ka2 = 0.0180311 * 10
kb1 = 1.0285752
kb_1 = Kd_wt * kb1
kb2 = 0.0180311 * 10
kc1 = 1.0285752
kc_1 = Kd_wt * kc1
kc2 = 0.0180311 * 10
kd1 = 1.0285752
kd_1 = Kd_wt * kd1
kd2 = 0.0180311 * 10
ke1 = 1.0285752
ke_1 = Kd_wt * ke1
ke2 = 0.0180311 * 10

theta = np.array([ka1, ka_1, ka2, kb1, kb_1, kb2, kc1, kc_1, kc2, kd1, kd_1, kd2, ke1, ke_1, ke2])

#initialised at initial values for test tube model
init_dicer = 5 * 10
init_premirna1 = 1
init_premirna2 = 1
init_premirna3 = 1
init_premirna4 = 1
init_premirna5 = 1
init_premirna1_dicer = 0
init_premirna2_dicer = 0
init_premirna3_dicer = 0
init_premirna4_dicer = 0
init_premirna5_dicer = 0
init_mirna1 = 0
init_mirna2 = 0
init_mirna3 = 0
init_mirna4 = 0
init_mirna5 = 0

init_values = np.array([init_dicer, init_premirna1, init_premirna2, init_premirna3, init_premirna4,
                        init_premirna5, init_premirna1_dicer, init_premirna2_dicer, init_premirna3_dicer,
                        init_premirna4_dicer, init_premirna5_dicer, init_mirna1, init_mirna2, init_mirna3,
                        init_mirna4, init_mirna5])

minutes = 60
dt = 0.1