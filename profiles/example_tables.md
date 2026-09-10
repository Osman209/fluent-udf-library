# Example input files for the profiles pack

Copy these into your Fluent working directory and edit.

## mu_T.txt - water viscosity against temperature

    # T [K]    mu [Pa s]
    273.15     0.001792
    283.15     0.001307
    293.15     0.001002
    303.15     0.000798
    313.15     0.000653
    323.15     0.000547
    333.15     0.000467
    343.15     0.000404
    353.15     0.000355
    363.15     0.000315
    373.15     0.000282

## transient_1.txt - a pressure that ramps then holds

    # t [s]    p [Pa gauge]
    0.0        0.0
    0.5        0.0
    2.0    50000.0
    10.0   50000.0

## rho_TP.txt - density on a temperature-pressure grid

    3 2
    280.0 300.0 320.0
    100000.0 200000.0
    999.9 1000.0
    996.5 996.6
    989.4 989.5
