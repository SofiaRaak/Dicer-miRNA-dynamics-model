import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp
import params
import utils
from tqdm import tqdm

kd_wt = params.Kd_wt

theta = [params.ka1, params.ka2]

solvers = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG', \
           'L-BFGS-B', 'TNC', 'COBYLA', 'SLSQP', 'trust-constr', \
           'dogleg', 'trust-ncg', 'trust-exec', 'trust-krylov']
bounds = ((0, np.inf), (0, np.inf))
solver_array = []
data_array = np.zeros(len(solvers))
extr_array = np.zeros(len(solvers))

data_dict = {}

#custom error function
def ErrorODE(theta, time, data_values = None):
    """
    Error function for ODE model based on error function described in utils module
    """
    
    ka, kc = theta
    
    kb = kd_wt * ka
    
    inits = [params.init_premirna1, params.init_premirna1_dicer, \
             params.init_dicer, params.init_mirna1]
    
    def model(t, inits):
        
        
        premiR = inits[1]*ka - inits[0]*inits[2]*kb
        premiR_dicer = inits[0]*inits[2]*kb - inits[1]*(ka + kc)
        dicer = inits[1]*(ka + kc) - inits[0]*inits[2]*kb
        mirna = inits[1]*kc
        
        return premiR, premiR_dicer, dicer, mirna
    
    sol = solve_ivp(model, (0, int(params.minutes)), inits)
    
    premiR = sol.y[0]
    
    diced = np.zeros(len(premiR))
    
    for i in range(len(diced)):
        diced[i] = (premiR[0] - premiR[i]) / premiR[0]
    
    if data_values is None:
        return utils.curveError(diced, time, sol.t)
    else:
        return utils.Error(data_values, diced, time, sol.t)

for i in tqdm(range(len(solvers))):
    
    utils.run_test()
    data_res = minimize(ErrorODE, theta, args = (params.time, params.WT_data),\
                        method = solvers[i], bounds=bounds)
    
    extr_res = minimize(ErrorODE, theta, args = (params.time), method = solvers[i], \
                       bounds=bounds)
    
    data_dict['data_'+solvers[i]] = data_res.x
    data_dict['extr_'+solvers[i]] = extr_res.x
    
    solver_array.append(solvers[i])
    data_array[i] = ErrorODE(data_res.x, params.time, params.WT_data)
    extr_array[i] = ErrorODE(extr_res.x, params.time)
    
df_error = pd.DataFrame(data_dict)
df_arrays = pd.DataFrame({"solver": solver_array, "data_points": data_array, "extrapolated": extr_array})

df_error.to_csv(r'optimisation/errors.csv', index=False)
df_arrays.to_csv(r'optimisation/optim_params.csv', index=False)