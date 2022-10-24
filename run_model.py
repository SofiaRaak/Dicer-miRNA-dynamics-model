import numpy as np
import matplotlib.pyplot as plt
import params
import utils
from scipy.integrate import solve_ivp
import pandas as pd
from tqdm import tqdm

#setup all same concentration, two populations with different concentrations
init_premirna = [
    [params.init_premirna1]*params.n_mirna,
    [params.init_premirna1]*int(params.n_mirna/2) + [params.init_premirna1*0.1]*int(params.n_mirna/2)
]

init_premirna_dicer = [
    [params.init_premirna1_dicer]*params.n_mirna,
    [params.init_premirna1_dicer]*int(params.n_mirna/2) + [params.init_premirna1_dicer*0.1]*int(params.n_mirna/2)
]

init_mirna = [
    [params.init_mirna1]*params.n_mirna, 
    [params.init_mirna1]*int(params.n_mirna/2) + [params.init_mirna1*0.1]*int(params.n_mirna/2)
]

conc_names = ["conc_hom", "conc_het"]
    
#setup all same reaction rate, two populations with different reaction rates
reaction_rates = [
    [[params.ka1]*params.n_mirna, [params.ka_1]*params.n_mirna, [params.ka2]*params.n_mirna],
    [[params.ka1]*int(params.n_mirna/2) + [params.ka1*params.k_dec]*int(params.n_mirna/2),
     [params.ka_1]*int(params.n_mirna/2) + [params.ka_1*params.k_dec]*int(params.n_mirna/2),
     [params.ka2]*int(params.n_mirna/2) + [params.ka2*params.k_dec]*int(params.n_mirna/2)]
]

rate_names = ["rate_hom", "rate_het"]

reaction_conditions = [[1]*3, 
                       [params.k_inc]+[1]*2, [params.k_dec]+[1]*2,
                       [1] + [params.k_inc] + [1], [1] + [params.k_dec] + [1],
                       [1]*2 + [params.k_inc], [1]*2 + [params.k_dec]]

reaction_names = ["basal", "high_ka", "low_ka", "high_kb", "low_kb", "high_kc", "low_kc"]

for i in tqdm(range(len(reaction_rates))):
    
    for j in tqdm(range(len(init_premirna))):
        inits, inits_dict = utils.make_inits(params.n_mirna, init_premirna[j], \
                                     init_premirna_dicer[j], \
                                     init_mirna[j], \
                                     params.init_dicer)
   
        reactions = reaction_rates[i]

        ks = np.zeros((len(reactions), params.n_mirna))
        for k in tqdm(range(len(reaction_conditions))):
            for l in range(len(reactions)):
                ks[l] = np.multiply(reactions[l], reaction_conditions[k][l])
                
            ODEs = utils.makeODEs(inits, ks[0], ks[1], ks[2])
            
            model = utils.makeModel(ODEs)
            
            exec(model)
            
            res = solve_ivp(runModel, (0, int(params.minutes)), inits, args = (ks[0], ks[1], ks[2]))
            
            mirna_length = int((len(res.y) - 1)/3)
            
            index = []
            
            for l in range(mirna_length):
                index = index + ["premirna" + str(l), "pmiR_dicer" + str(l), "mirna" + str(l)]
            
            index = index + ["dicer"]
                
            
            df = pd.DataFrame(res.y, index = index)
            
            df = df.transpose()
            
            df.insert(loc = 0, column = "time_min", value = res.t)
            
            name = "1pDicer_" + rate_names[i] + "_" + conc_names[j] + "_" + reaction_names[k] + ".csv"
            
            df.to_csv(r'./output/'+name, index = False)