def make_inits(num_preMiR, pMiR, dicer_pMiR, miR, dicer):
    """
    Function to generate vector of initial concentrations for all species.
    
    Args:
    num_preMiR:    int, number of pre-miRNA species
    pMiR:          list of floats, initial concentration of pre-miRNA
    dicer_pMiR:    list of floats, initial concentration Dicer and pre-miRNA complexes
    miR:           list of floats, initial concentration of mature miRNA
    dicer:         float, initial concentration of dicer
    
    Returns:
    list of initial concentrations of (pre-miR, pre-miR x dicer, miR) * n and dicer
    dictionary of initial concentrations of each species
    """
    
    assert num_preMiR == len(pMiR), "The number of concentrations given do \
    not correspond with the number of species"
    assert num_preMiR == len(dicer_pMiR), "The number of concentrations given do not \
    correspond with the number of species"
    
    init_concs = [pMiR[0], dicer_pMiR[0], miR[0]]
    init_species = {"pMiR0": pMiR[0], "dicer_pMiR0": dicer_pMiR[0], "miR": miR[0]}
    
    for i in range(1, num_preMiR):
        init_concs.append(pMiR[i])
        init_concs.append(dicer_pMiR[i])
        init_concs.append(miR[i])
        
        init_species["preMiR"+str(i)] = pMiR[i]
        init_species["dicer_preMiR"+str(i)] = dicer_pMiR[i]
        init_species["miR"+str(i)] = miR[i]
    
    init_concs.append(dicer)
    init_species["dicer"] = dicer
        
    return init_concs, init_species

def makeODEs(init_values, ka, kb, kc):
    """
    Function to generate string of functions for passing to ODE solver
    
    init_values:    List of initial concentration of free pre-miRNA concentration for
                    each pre-mirna, pre-mirna x dicer concentration, mirna concentration
                    and dicer concentration. Alternating pre-mirna (a), pre-mirna dicer (b),
                    and mirna(c), with final dicer concentration (d), i.e. [a, b, c, a, b, c, d]
    ka:             Array of reaction rates for formation of pre-miR - dicer commplex
    kb:             Array of reaction rates for dissociation of pre-miR - dicer complex
    kc:             Array of reaction rates for catalysis of pre-miR to miR by dicer
    
    Outputs a string with all functions separated by new lines.
    """
    
    assert (len(init_values)-1) % 3 == 0, "You have supplied an uneven number of concentrations"
    assert int((len(init_values) - 1)/3) == len(ka), "The number of reactions do not \
    agree with the number of reaction rates given."
    assert int((len(init_values) - 1)/3) == len(kb), "The number of reactions do not \
    agree with the number of reaction rates given."
    assert int((len(init_values) - 1)/3) == len(kc), "The number of reactions do not \
    agree with the number of reaction rates given."
    
    #calculate number of miRNAs
    n_mirna = int((len(init_values) - 1)/3)
    output = ""
    names = ""
    concs = ""
    dicer_strings = "        dicer = "
    
    #create system
    j = 0
    for i in range(n_mirna):

        output = output + '\n        premirna' + str(i) + ' = init_values[' + str(j+1) + '] * kb[' + str(i) + '] - init_values[' + str(j) + '] * init_values[-1] * ka[' + str(i) + ']' + '\n        pmiR_dicer' + str(i) + ' = init_values[' + str(j) + '] * init_values[-1] * ka[' + str(i) + '] - init_values[' + str(j+1) + '] * (kb[' + str(i) + '] + kc[' + str(i) + '])' + '\n        mirna' + str(i) + ' = init_values[' + str(j+1) + '] * kc[' + str(i) + ']'
        
        if i == 0:
            dicer_strings = dicer_strings + 'init_values[' + str(j+1) + '] * (kb[' + str(i) + '] + kc[' + str(i) + ']) - init_values[-1] * init_values[' + str(j) + '] * ka[' + str(i) + ']'
            names = names + 'premirna' + str(i) + ', pmiR_dicer' + str(i) + ', mirna' + str(i)
        
        else:
            dicer_strings = dicer_strings + '+ init_values[' + str(j+1) + '] * (kb[' + str(i) + '] + kc[' + str(i) + ']) - init_values[-1] * init_values[' + str(j) + '] * ka[' + str(i) + ']'
            names = names + ', premirna' + str(i) + ', pmiR_dicer' + str(i) + ', mirna' + str(i)
        
        j += 3
    
    output = output + '\n' + dicer_strings
    names = names + ', dicer'
    
    return output, names

def makeModel(ODEs):
    """
    Function to generate function runModel(t, init_values, ka, kb, kc)
    
    Args:
    ODEs:    List of strings with all ODEs in system, output variable names.
    
    Returns string version of final model to be passed to exec function.
    """
    
    model = '''def runModel(t, init_values, ka, kb, kc):
        """
        Function to run model in ODE solver
        
        init_values:    List of initial concentration of free pre-miRNA concentration for
                        each pre-mirna, pre-mirna x dicer concentration, mirna concentration
                        and dicer concentration. Alternating pre-mirna (a), pre-mirna dicer (b),
                        and mirna(c), with final dicer concentration (d), i.e. [a, b, c, a, b, c, d]
        ka:             Array of reaction rates for formation of pre-miR - dicer commplex
        kb:             Array of reaction rates for dissociation of pre-miR - dicer complex
        kc:             Array of reaction rates for catalysis of pre-miR to miR by dicer
        """
        
        
        
        ''' + ODEs[0] + '''
        
        return ''' + ODEs[1] 
    
    return model