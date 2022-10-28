# Dynamic miR-Dicer model

Computational model describing interactions between a single pool of Dicer and n pre-miRNAs and their maturation.

## Model background

Parameters have been obtained from [Tsutsumi _et. al_ 2011](https://www.nature.com/articles/nsmb.2125) and through optimisation

## Model structure

The model is designed as a system of ordinary differential equations where each species occupies a specific index:

```math
\frac{dpMiR_i}{dt} = Dicer\_pMiR_i*k_{b_i} - Dicer*pMiR_i*k_{a_i}

\frac{dDicer\_pMiR_i}{dt} = Dicer*pMiR_i*k_{a_i} - Dicer\_pMiR_i*(k_{b_i}+k_{c_i})

\frac{dMiR_i}{dt} = Dicer\_pMiR_i*k_{c_i}

\frac{dDicer}{dt} = \sum^{n}_{i=0} Dicer\_pMiR_i*(k_{b_i}+k_{c_i}) - \sum^{n}_{i=0} Dicer*pMiR_i*k_{a_i}
```

where

* $`pMiR_i`$ is the concentration of pre-miRNA species at index $`i`$ ($`nM`$)
* $`Dicer\_pMiR_i`$ is the concentration of pre-miRNA Dicer complexes at index &`i`& ($`nM`$)
* $`MiR_i`$ is the concentration of mature miRNA species at index $`i`$ ($`nM`$
* $`Dicer`$ is the concentraton of free Dicer in the system ($`nM`$)
* $`n`$ is the number of miRNA species in the system
* $`k_{a_i}`$ is the rate of association of pre-miRNA at index $`i`$ and Dicer ($`s^{-1}`$)
* $`k_{b_i}`$ is the rate of dissociation of pre-miRNA at index $`i`$ and Dicer ($`s^{-1}`$)
* $`k_{c_i}`$ is the maturation rate of miRNA at index $`i`$ ($`s^{-1}`$
