"""
ODE based model of miRNA and Dicer interactions

Background and data obtained from Tsutsumi et al. 2011, Nat Struct Mol Biol, 10.1038/nsmb.2125

"""

import numpy as np

def test():
    print('model_utils imported OK!')
    
    
def make_inits(num_preMiR, pMiR, dicer_pMiR, dicer_pMiR_star, miR, dicer):
    """
    Function to generate vector of initial concentrations for all species.
    
    Args:
    num_preMiR:      int, number of pre-miRNA species
    pMiR:            list of floats, initial concentration of pre-miRNA
    dicer_pMiR:      list of floats, initial concentration Dicer and pre-miRNA complexes
    dicer_pMiR_star: list of floats, initial concentration Dicer and pre-miRNA complexes that \
                     cannot dissociate
    miR:             list of floats, initial concentration of mature miRNA
    dicer:           float, initial concentration of dicer
    
    Returns:
    list of initial concentrations of (pre-miR, pre-miR x dicer, miR) * n and dicer
    dictionary of initial concentrations of each species
    """
    
    assert num_preMiR == len(pMiR), "The number of concentrations given do \
    not correspond with the number of species"
    assert num_preMiR == len(dicer_pMiR), "The number of concentrations given do not \
    correspond with the number of species"
    
    init_concs = [pMiR[0], dicer_pMiR[0], dicer_pMiR_star[0], miR[0]]
    init_species = {"pMiR0": pMiR[0], "dicer_pMiR0": dicer_pMiR[0], \
                    "dicer_pMiR0_star": dicer_pMiR0_star[0], "miR": miR[0]}
    
    for i in range(1, num_preMiR):
        init_concs.append(pMiR[i])
        init_concs.append(dicer_pMiR[i])
        init_concs.append(dicer_pMiR_star[0])
        init_concs.append(miR[i])
        
        init_species["preMiR"+str(i)] = pMiR[i]
        init_species["dicer_preMiR"+str(i)] = dicer_pMiR[i]
        init_species["miR"+str(i)] = miR[i]
    
    init_concs.append(dicer)
    init_species["dicer"] = dicer
        
    return init_concs, init_species

def makeODEs(init_values, ka, kb, kc, kd):
    """
    Function to generate string of functions for passing to ODE solver
    
    init_values:    List of initial concentration of free pre-miRNA concentration for
                    each pre-mirna, pre-mirna x dicer concentration, pre-mirna x dicer* 
                    concentration, mirna concentration, and dicer concentration. 
                    Alternating pre-mirna (a), pre-mirna dicer (b), pre-mirna dicer* (c)
                    and mirna(d), with final dicer concentration (e), i.e. 
                    [a, b, c, d, a, b, c, d, e]. Can be automatically generated with fun
                    make_inits().
    ka:             Array of reaction rates for formation of pre-miR - dicer commplex
    kb:             Array of reaction rates for dissociation of pre-miR - dicer complex
    kc:             Array of reaction rates for formation of pre-miR - dicer complex 
    kd:             Array of reaction rates for catalysis of pre-miR to miR by dicer
    
    Outputs a string with all functions separated by new lines.
    """
    
    assert (len(init_values)-1) % 4 == 0, "You have supplied an uneven number of concentrations"
    assert int((len(init_values) - 1)/4) == len(ka), "The number of reactions do not \
    agree with the number of reaction rates given."
    assert int((len(init_values) - 1)/4) == len(kb), "The number of reactions do not \
    agree with the number of reaction rates given."
    assert int((len(init_values) - 1)/4) == len(kc), "The number of reactions do not \
    agree with the number of reaction rates given."
    assert int((len(init_values) - 1)/4) == len(kd), "The number of reactions do not \
    agree with the number of reaction rates given."
    
    #calculate number of miRNAs
    n_mirna = int((len(init_values) - 1)/4)
    output = ""
    names = ""
    concs = ""
    dicer_strings = "        dicer = "
    
    #create system
    j = 0
    for i in range(n_mirna):
        
        #init_values[0] = pMiR
        #init_values[1] = pMiR_dicer
        #init_values[2] = pMiR_dicer_star >> never in ODE, needed for solver
        #init_values[3] = miR >> never in ODE, needed for solver
        #init_values[-1] = dicer

        output = output + '\n\n        premirna' + str(i) + ' = init_values[' + str(j+1) + '] * kb[' + str(i) + '] - init_values[' + str(j) + '] * init_values[-1] * (ka[' + str(i) + '] + kc[' + str(i) + '])' + '\n        pmiR_dicer' + str(i) + ' = init_values[' + str(j) + '] * init_values[-1] * ka[' + str(i) + '] - init_values[' + str(j+1) + '] * (kb[' + str(i) + '] + kd[' + str(i) + '])' + '\n        pmiR_dicer_star' + str(i) + ' = init_values[' + str(j) + '] * init_values[-1] * kc[' + str(i) + ']' + '\n        mirna' + str(i) + ' = init_values[' + str(j+1) + '] * kd[' + str(i) + ']'
        
        if i == 0:
            dicer_strings = dicer_strings + 'init_values[' + str(j+1) + '] * (kb[' + str(i) + '] + kd[' + str(i) + ']) - init_values[-1] * init_values[' + str(j) + '] * (ka[' + str(i) + '] + kc[' + str(i) + '])'
            names = names + 'premirna' + str(i) + ', pmiR_dicer' + str(i) + ', pmiR_dicer_star' + str(i) + ', mirna' + str(i)
        
        else:
            dicer_strings = dicer_strings + '+ init_values[' + str(j+1) + '] * (kb[' + str(i) + '] + kd[' + str(i) + ']) - init_values[-1] * init_values[' + str(j) + '] * (ka[' + str(i) + '] + kc[' + str(i) + '])'
            names = names + ', premirna' + str(i) + ', pmiR_dicer' + str(i) + ', pmiR_dicer_star' + str(i) + ', mirna' + str(i)
        
        j += 4
    
    output = output + '\n\n' + dicer_strings
    names = names + ', dicer'
    
    return output, names

def makeModel(ODEs):
    """
    Function to generate string of functions for passing to ODE solver
    
    Args:
    ODEs:    List of strings with all ODEs in system, output variable names.
    
    Returns string version of final model to be passed to exec function.
    """
    
    model = '''def runModel(t, init_values, ka, kb, kc, kd):
        """
        Function to run model in ODE solver
        
        Function to generate string of functions for passing to ODE solver
    
        init_values:    List of initial concentration of free pre-miRNA concentration for
                        each pre-mirna, pre-mirna x dicer concentration, pre-mirna x dicer* 
                        concentration, mirna concentration, and dicer concentration. 
                        Alternating pre-mirna (a), pre-mirna dicer (b), pre-mirna dicer* (c)
                        and mirna(d), with final dicer concentration (e), i.e. 
                        [a, b, c, d, a, b, c, d, e]. Can be automatically generated with fun
                        make_inits().
        ka:             Array of reaction rates for formation of pre-miR - dicer commplex
        kb:             Array of reaction rates for dissociation of pre-miR - dicer complex
        kc:             Array of reaction rates for formation of pre-miR - dicer complex 
        kd:             Array of reaction rates for catalysis of pre-miR to miR by dicer
        """
        
        
        
        ''' + ODEs[0] + '''
        
        return ''' + ODEs[1] 
    
    return model


def createGaussian(n, mu, sigma, noise = 0):
    """
    Function to create normal distribution of input array with noise
    
    Args:
    n:        size of generated data
    mu:       mean of distribution
    sigma:    standard deviation of generated data
    noise:    introduction of noise into the system. Default 0
    
    Returns:
    array of normally distributed values of size n
    """
    
    gauss = np.random.normal(loc = mu, scale = sigma, size = n)
    
    return gauss + noise * np.random.normal(size = len(gauss))

def createGamma(n, mu, sigma, noise = 0):
    """
    Function to create Gamma distribution of input array with noise
    
    Args:
    n:        size of generated data
    mu:       mean of distribution
    sigma:    standard deviation of generated data
    noise:    introduction of noise into the system. Default 0
    
    Returns:
    array of gamma distributed values of size n
    """
    
    shape = (mu**2)/sigma
    scale = sigma/mu
    
    gamma = np.random.gamma(shape = shape, scale = scale, size = n)
    
    return gamma + noise * np.random.gamma(shape = shape, scale = scale, size = n)

def createLogNormal(n, mu, sigma, noise = 0):
    """
    Function to create log-normal distribution of input data with noise
    
    Args:
    n:        size of generated data
    mu:       mean of distribution
    sigma:    standard deviation of generated data
    noise:    introduction of noise into the system. Default 0
    
    Returns:
    array of log-normal distributed values of size n
    """
    
    lognorm = np.random.lognormal(mean = mu, sigma = sigma, size = n)
    
    return lognorm + noise * np.random.lognormal(mean = mu, sigma = sigma, size = n)

