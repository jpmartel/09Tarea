"""
Este script realiza un modelo lineal para la relacion entre el flujo en la
banda i y la banda z de cuasares con datos de Data Release 9 (data/DR9Q.dat)
del survey SDSS.
Incluye intervalo de confianza al 95% utilizando una simulacion de Monte Carlo
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
np.random.seed(21)


def Monte_Carlo(Bi, Bz, Bi_error, Bz_error):
    N_mc = 10000
    pendiente = np.zeros(N_mc)
    constante = np.zeros(N_mc)
    N = len(Bi)
    for i in range(N_mc):
        r = np.random.normal(0, 1, size=N)
        muestra_Bi = Bi + Bi_error*r
        muestra_Bz = Bz + Bz_error*r
        pendiente[i], constante[i] = np.polyfit(muestra_Bi, muestra_Bz, 1)
    pendiente = np.sort(pendiente)
    constante = np.sort(constante)
    lim_inferior_pendiente = pendiente[int(N_mc*(alpha/2))]
    lim_superior_pendiente = pendiente[int(N_mc*(1-alpha/2))]
    lim_inferior_constante = constante[int(N_mc*(alpha/2))]
    lim_superior_constante = constante[int(N_mc*(1-alpha/2))]
    return [[lim_inferior_pendiente, lim_superior_pendiente],
            [lim_inferior_constante, lim_superior_constante]]

# Setup

datos = np.loadtxt('data/DR9Q.dat', usecols=(80, 81, 82, 83))
Bi = datos[:, 0]*3.631
Bi_error = datos[:, 1]*3.631
Bz = datos[:, 2]*3.631
Bz_error = datos[:, 3]*3.631

# Ajuste lineal
a, b = np.polyfit(Bi, Bz, 1)
ajuste = a*Bi+b
print 'constantes ajuste lineal: '
print 'pendiente = ' + str(a)
print 'constante = ' + str(b)
print 'ajuste: Bz = ' + str(a) + '*Bi + ' + str(b)

# Intervalo de confianza al 95%

alpha = 0.05
Intervalo_confianza = Monte_Carlo(Bi, Bz, Bi_error, Bz_error)
[lim_inferior_pendiente, lim_superior_pendiente] = Intervalo_confianza[0]
[lim_inferior_constante, lim_superior_constante] = Intervalo_confianza[1]
print 'intervalo de confianza pendiente = ' \
                                    + '[ ' + str(lim_inferior_pendiente) \
                                    + ' , ' + str(lim_superior_pendiente) + ']'
print 'intervalo de confianza constante = ' \
                                    + '[ ' + str(lim_inferior_constante) \
                                    + ' , ' + str(lim_superior_constante) + ']'

# Grafico

plt.plot(Bi, Bz, '*', label='datos Data Release 9')
plt.plot(Bi, ajuste, label='ajuste lineal')
plt.xlabel('flujo banda i [ 1e-6 Jy ]')
plt.ylabel('flujo banda z [ 1e-6 Jy ]')
plt.title('Relacion flujo banda i y banda z')
plt.legend(loc='upper left')
plt.show()
