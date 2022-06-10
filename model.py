import numpy as np
from scipy.integrate import solve_ivp
import tqdm
import params

def toy_model(t, init_values): #not bothering with passing reaction rates since not optimising
    """
    Toy model consisting of 5 miRNAs competing for the same Dicer pool
    """
    
    init_dicer, init_premirna1, init_premirna2, init_premirna3, init_premirna4, init_premirna5, \
    init_premirna1_dicer, init_premirna2_dicer, init_premirna3_dicer, init_premirna4_dicer, \
    init_premirna5_dicer, init_mirna1, init_mirna2, init_mirna3, init_mirna4, init_mirna5 = init_values #this needs to be automated if using >20 miRNAs
    
    #Proper model
    premirna1 = init_premirna1_dicer*ka_1 - init_premirna1*init_dicer*ka1
    premirna2 = init_premirna2_dicer*kb_1 - init_premirna2*init_dicer*kb1
    premirna3 = init_premirna3_dicer*kc_1 - init_premirna3*init_dicer*kc1
    premirna4 = init_premirna4_dicer*kd_1 - init_premirna4*init_dicer*kd1
    premirna5 = init_premirna5_dicer*ke_1 - init_premirna5*init_dicer*ke1
    
    dicer = (init_premirna1_dicer*(ka_1 + ka2) + init_premirna2_dicer*(kb_1 + kb2) + 
             init_premirna3_dicer*(kc_1 + kc2) + init_premirna4_dicer*(kd_1 + kd2) +
             init_premirna5_dicer*(ke_1 + ke2)) - \
            init_dicer*(init_premirna1*ka1 + init_premirna2*kb1 + init_premirna3+kc1 +
                        init_premirna4*kd1 + init_premirna5*kd1)
    
    mirna1 = init_premirna1_dicer*ka2
    mirna2 = init_premirna2_dicer*kb2
    mirna3 = init_premirna3_dicer*kc2
    mirna4 = init_premirna4_dicer*kd2
    mirna5 = init_premirna5_dicer*ke2
    
    return dicer, premirna1, premirna2, premirna3, premirna4, premirna5, premirna1_dicer, \
           premirna2_dicer, premirna3_dicer, premirna4_dicer, premirna5_dicer, mirna1, mirna2, \
           mirna3, mirna4, mirna5
                                                              