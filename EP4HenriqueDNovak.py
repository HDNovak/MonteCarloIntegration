import numpy as np
import math
import random
from matplotlib import pyplot as plt
from IPython.display import clear_output

PI = 3.1415926
e = 2.71828
C = 1.398147061
R = 1.10349828

from math import *
from RandomArray import *
from matplotlib.pylab import *

def sdnorm(s):
    """
    Normalização da PDF, tendo em vista sua definição exclusivamente positiva.
    Núcleo normal esperado.
    """
    return exp(-z*z/2.)/sqrt(2*pi)

n = 10000
alpha = 2
x = 0
k = 0
vec = []
bec = []
bec.append(k)
vec.append(x)
innov= np.random.uniform(0, alpha, n) #Parte aleatória para a proporção
for i in xrange(1,n):
    can = x + innov[i] #candidato
    aprob = min([1.,sdnorm(can)/sdnorm(x)]) #Metrópolis
    u = uniform(0,1)
    if u < aprob:
        x = can
        vec.append(x)
    vec = np.array(vec)
    return vec

for i in krange(1,n):
    can2 = k + innov[i] #candidate
    bprob = [sdnorm(can)/(sdnorm(k) + sdnorm(can))] #Barker
    u = uniform(0,1)
    if u < bprob:
        k = can2
        vec.append(k)
    bec = np.array(bec)
    return bec
"""
#Algoritimo para plotar a curva

x = arange(-3,3,.1)
y = sdnorm(x)
subplot(211)
title('Metropolis-Hastings')
plot(vec)
subplot(212)

hist(vec, bins=30,normed=1)
plot(x,y,'ro')
ylabel('Frequency')
xlabel('x')
legend(('PDF','Samples'))
show()

"""
# Função integranda f(x)

def f_of_x(x):
    """
    Essa é a função do qual queremos integrar
    """
    
    return x*([C*(e**(-Cx))]*abs(math.cos(Rx)))

def numero_aleatorioMetropolis(min_value, max_value):
    """
    Esse função encontra um valor aleatório entre a distribuição de ambos inputs
    - min_value (float)
    - max_value (float)
    Return:
    - Número aleatório entre esse intervalo (float)
    """
    range = max_value - min_value
    choice = vec
    return min_value + range*choice

def numero_aleatorioBarker(min_value, max_value):
    """
    Esse função encontra um valor aleatório entre a distribuição de ambos inputs
    - min_value (float)
    - max_value (float)
    Return:
    - Número aleatório entre esse intervalo (float)
    """
    range = max_value - min_value
    choice = bec
    return min_value + range*choice

def crude_monte_carlo_metropolis(num_samples=250):
    """
    Essa função faz o nosso Crude Monte Carlo
    Basicamente gera um loop que através do somátorio dos valores aleatórios
    contabiliza o nosso valor estimado
    
    """
    lower_bound = 1
    upper_bound = 2
    
    sum_of_samples = 0
    for i in range(num_samples):
        x = numero_aleatorioMetropolis(lower_bound, upper_bound)
        sum_of_samples += f_of_x(x)
    
    return (upper_bound - lower_bound) * float(sum_of_samples/num_samples)

def crude_monte_carlo_barker(num_samples=250):
    """
    Essa função faz o nosso Crude Monte Carlo
    Basicamente gera um loop que através do somátorio dos valores aleatórios
    contabiliza o nosso valor estimado
    
    """
    lower_bound = 1
    upper_bound = 2
    
    sum_of_samples = 0
    for i in range(num_samples):
        x = numero_aleatorioBarker(lower_bound, upper_bound)
        sum_of_samples += f_of_x(x)
    
    return (upper_bound - lower_bound) * float(sum_of_samples/num_samples)

def get_crude_MC_varience_metropolis(num_samples):
    """
    Essa função retorna a variância do nosso monte_carlo
    Simplesmente gerando o valor esperado
    """
    int_max = 1 
    
    # média dos quadrados através do loop
    running_total = 0
    for i in range(num_samples):
        x = numero_aleatorioMetropolis(0, int_max)
        running_total += f_of_x(x)**2
    sum_of_sqs = running_total*int_max / num_samples
    
    # quadrado da média com o loop
    running_total = 0
    for i in range(num_samples):
        x = numero_aleatorioMetropolis(0, int_max)
        running_total = f_of_x(x)
    sq_ave = (int_max*running_total/num_samples)**2
    
    return sum_of_sqs - sq_ave

def get_crude_MC_varience_barker(num_samples):
    """
    Essa função retorna a variância do nosso monte_carlo
    Simplesmente gerando o valor esperado
    """
    int_max = 1 
    
    # média dos quadrados através do loop
    running_total = 0
    for i in range(num_samples):
        x = numero_aleatorioBarker(0, int_max)
        running_total += f_of_x(x)**2
    sum_of_sqs = running_total*int_max / num_samples
    
    # quadrado da média com o loop
    running_total = 0
    for i in range(num_samples):
        x = numero_aleatorioBarker(0, int_max)
        running_total = f_of_x(x)
    sq_ave = (int_max*running_total/num_samples)**2
    
    return sum_of_sqs - sq_ave


# A simulação foi feita diversas vezes com valores <300, sendo 300 o que gerou o melhor erro

MC_samples = 250
var_samples = 250 # Valores usados para a variância
crude_barker = crude_monte_carlo_barker(MC_samples)
crude_metropolis = crude_monte_carlo_metropolis(MC_samples)
variance_barker = get_crude_MC_variance_barker(var_samples)
variance_metropolis = get_crude_MC_variance_metropolis(var_samples)
error_barker = math.sqrt(variance_barker/MC_samples)
error_metropolis = math.sqrt(variance_metropolis/MC_samples)


# Resultados
print(f"Crude Monte Carlo: {crude_barker}")
print(f"Crude Monte Carlo: {crude_metropolis}")
print(f"Variância: {variance_metropolis}")
print(f"Variância: {variance_barker}")
print(f"Erro estimado: {error_barker}")
print(f"Erro estimado: {error_metropolis}")

