"""
Este string estima la constante de Hubble a traves de un ajuste lineal,
y entrega un intervalo de confianza al 95% utilizando una una simulacion
Bootstrap. Utiliza datos de data/SNIa.dat (Freedman et al. 2000)
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
np.random.seed(100)


def bootstrap(v, D, alpha):
    """Retorna una estimacion bootstrap para un intervalo de confianza
    al 100(1-alpha)%"""
    N = len(v)
    Nboot = int(10000)
    H_values = np.zeros(Nboot)
    for i in range(Nboot):
        indice = np.random.randint(low=0, high=N, size=N)
        v_muestra = v[indice]
        D_muestra = D[indice]
        H_values[i] = ajuste(D_muestra, v_muestra)
    H_values = np.sort(H_values)
    limite_inferior = H_values[int(Nboot*(alpha/2))]
    limite_superior = H_values[int(Nboot*(1 - alpha/2))]
    fig2 = plt.figure(2)
    plt.hist(H_values, bins=30, normed=True, label='histograma para H0')
    plt.axvline(np.mean(H_values), color='r', label='valor optimo')
    plt.legend(loc='upper left')
    plt.xlabel('H0 [km/s/Mpc]')
    plt.title('Histograma de H0 con datos de Freedman (et al.)')
    return limite_inferior, limite_superior


def ajuste(D, v):
    H0_1 = np.sum(v*D)/np.sum(D**2)
    H0_2 = np.sum(v**2)/np.sum(v*D)
    H0_promedio = (H0_1 + H0_2)/2
    return H0_promedio

# setup

datos = np.loadtxt('data/SNIa.dat', usecols=(1, 2))
D = datos[:, 1]  # distancia
v = datos[:, 0]  # velocidad

# Constante de Hubble obtenida de metodo minimos cuadrados
H0 = ajuste(D, v)
recta = D*H0

# Intervalo de confianza al 95% de H0

alpha = 0.05
limite_inferior, limite_superior = bootstrap(v, D, alpha)
print 'H0 = ' + str(H0)
print 'intervalo de confianza = ' + '[ ' + str(limite_inferior) + ' , ' \
                                + str(limite_superior) + ']'

# Grafico

fig = plt.figure(1)
plt.plot(D, v, '*', label='datos Freedman (et al.)')
plt.plot(D, recta, label='ajuste lineal')
plt.xlabel('distancia [Mpc]')
plt.ylabel('velocidad de recesion [km/s]')
plt.title('Ley de Hubble datos Freedman (et al)')
plt.legend(loc='upper left')
plt.show()
