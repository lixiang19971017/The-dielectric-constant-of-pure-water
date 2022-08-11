import numpy as np


def PureWater(T, freq):
    """
The dielectric constant and the complex refractive index of pure water.
Background information: Meissner, Thomas, and Frank J. Wentz.
                        "The complex dielectric constant of pure and sea water from microwave satellite observations."
                        IEEE Transactions on Geoscience and remote Sensing 42.9 (2004): 1836-1849.
Epsilon_water, Kw2 = PureWater(T, freq)
    T:          temperature [Centigrade]
    freq:       frequency [GHz]
"""
    Epsilon_s = (3.70886 * pow(10, 4) - 8.2168 * pow(10, 1) * T) \
                / (4.21854 * pow(10, 2) + T)
    Epsilon_1 = 5.7230 + (2.2379 * pow(10, -2) * T) + (-7.1237 * pow(10, -4) * T ** 2)
    Epsilon_inf = 3.6143 + (2.8841 * pow(10, -2) * T)
    Nu_1 = (45 + T) / (5.0478 + (-7.0315 * pow(10, -2) * T) + (6.0059 * pow(10, -4) * T ** 2))
    Nu_2 = (45 + T) / (1.3652 * pow(10, -1) + (1.4825 * pow(10, -3) * T) + (2.4166 * pow(10, -4) * T ** 2))

    Epsilon_water = (Epsilon_s - Epsilon_1) / (1 + 1j * freq / Nu_1) + \
                    (Epsilon_1 - Epsilon_inf) / (1 + 1j * freq / Nu_2) + Epsilon_inf

    Kw2 = np.abs((Epsilon_water - 1) / (Epsilon_water + 2)) ** 2
    return Epsilon_water, Kw2
