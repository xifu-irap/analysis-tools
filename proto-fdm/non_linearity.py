import numpy as np
import matplotlib.pyplot as plt


def calculate_INL_DNL(measured_values):

    n_codes = len(measured_values)
    digital_codes = np.arange(n_codes)

    # Calcul de l'écart entre deux valeurs de code consécutives
    ideal_step = (measured_values[-1] - measured_values[0]) / (n_codes - 1)

    # Calcul des valeurs différentielles pour les DNL
    dnl_values = (measured_values[1:] - measured_values[:-1]) / ideal_step - 1

    # Calcul de l'INL en intégrant les valeurs DNL
    inl_values = np.cumsum(dnl_values)

    # Tracé des courbes avec des styles personnalisés
    plt.figure(figsize=(10, 6))

    # Courbe INL
    plt.subplot(2, 1, 1)
    plt.plot(digital_codes[1:], inl_values, marker='o', label='INL', color='b')
    plt.xlabel('Digital code')
    plt.ylabel('INL (LSB)')
    plt.title('Integral Nonlinearity (INL)')
    plt.grid()

    # Courbe DNL
    plt.subplot(2, 1, 2)
    plt.plot(digital_codes[1:], dnl_values, marker='s', label='DNL', color='r')
    plt.xlabel('Digital code')
    plt.ylabel('DNL (LSB)')
    plt.title('Differential Nonlinearity (DNL)')
    plt.grid()

    plt.tight_layout()
    plt.show()

    return inl_values, dnl_values
