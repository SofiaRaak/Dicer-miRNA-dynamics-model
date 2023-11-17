

# Dynamic miR-Dicer model

Computational model describing interactions between a single pool of Dicer and n pre-miRNAs and their maturation.

## Model background

Parameters have been obtained from [Tsutsumi _et. al_ 2011](https://www.nature.com/articles/nsmb.2125) and through optimisation

## Model structure

The model is designed as a system of ordinary differential equations where each species occupies a specific index:

```math

\frac{dpMiR_i}{dt} = pMiR_i\_Dicer * k_{b_i} - pMiR_i * Dicer * k_{a_i}
```

```math
\frac{dpMiR_i\_Dicer}{dt} = pMiR_i * Dicer * k_{a_i} - pMiR_i\_Dicer * (k_{b_i} + k_{c_i} + k_{d_i})
```

```math
\frac{dDicer}{dt} = pMiR_i\_Dicer * (k_{b_i} + k_{d_i}) - pMiR_i * Dicer * k_{a_i}
```

```math
\frac{dpMiR_i\_Dicer^*}{dt} = pMiR_i\_Dicer * k_{c_i}
```

```math
\frac{dmiR_i}{dt} = pMIR_i\_Dicer * k_{d_i}
```

where

* $pMiR_i$ is the concentration of pre-miRNA species at index $i$ $(nM)$
* $Dicer\_pMiR_i$ is the concentration of pre-miRNA Dicer complexes at index $i$ $(nM)$
* $MiR_i$ is the concentration of mature miRNA species at index $i$ $(nM)$
* $Dicer$ is the concentraton of free Dicer in the system $(nM)$
* $n$ is the number of miRNA species in the system
* $k_{a_i}$ is the rate of association of pre-miRNA at index $i$ and Dicer $(s^{-1})$
* $k_{b_i}$ is the rate of dissociation of pre-miRNA at index $i$ and Dicer $(s^{-1})$
* $k_{c_i}$ is the stalling rate of pre-miRNA at index $i$ $(s^{-1})$
* $k_{d_i}$ is the maturation rate of pre-miRNA at index $i$ $(s^{-1})$
