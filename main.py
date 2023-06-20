# -*- coding: utf-8 -*-
"""
@author: Guilhermee Peinador Gomes
"""
import numpy as np
#import sympy as smp
import matplotlib.pyplot as plt

''''Motor escolhido P&W PW4098'''

'''CÁLCULO DO TURBOFAN'''

K = 273.15                      #Temperatura kelvin
Mo = np.linspace(0.1, 1, 10)     #Mach de voo variando de 0.1 a 1
pif = 1.7                       #Taxa de compressão LPC Motor escolhido = 1.7 - 1.8
pic = 34.2                      #Taxa de compressão HPC Motor Escolhido = 34.2 - 42.8
T04 = 1944                      #Temperatura máxima câmara combustão
PCI = 42.8E6                    #Poder calorifico
cp = 1004
gamma = 1.4
Po = 22700                      #Pressão na entrada **** Colocar no relatório
To = -30 + K                    #Temperatura em K na entrada
B = 5.8                         #Bypass Motor Escolhido = 5.8 - 6.4
R = 287                         #Constante R
PSL= 101325                     #Pressão em nível do mar


'''Intake'''

P01 = Po*(1+((gamma-1)/2)*Mo**2)**(gamma/(gamma-1))
T01 = To*(1+((gamma-1)/2)*Mo**2)

'''Fan LPC'''

P02f = P01
T02f = T01

P03f = pif*P02f
tauf = pif**((gamma-1)/gamma)
T03f = tauf*T02f
#print(tauf)

'''Compressor HPC'''

P03 = P03f*pic
tauc = pic**((gamma-1)/gamma)
T03 = T03f*tauc
#print(tauc)

'''Câmara de combustão'''

P04 = P03
f = ((T04/T03)-1)/((PCI/(cp*T03))-(T04/T03))
#print('f', f)

'''Conectar HPT-HPC'''

T05 = T04-((T03-T03f)/(1+f))
P05 = P04*((T05/T04)**((gamma)/(gamma-1)))

'''Conectar LPT-LPC'''

T06 = T05-((1+B)*(T03f-T02f))/(1+f)
P06 = P05*(T06/T05)**((gamma)/(gamma-1))

'''Bocal 3' - 7' '''

P07f = P03f

NPRf = P07f/Po
Pef = P07f/1.892

'''Bocal 6-7'''

P07 = P06

NPR = P07/Po
Pe  = P07/1.892

Tef = T03f/(1+((gamma-1)/2))
Te  = T06/(1+((gamma-1)/2))

Uef = (gamma*R*(Tef))**(1/2)
Ue  = (gamma*R*(Te))**(1/2)
Uo  = Mo*((gamma*R*To)**(1/2))

rhof = Pef/(R*T03f)
rho  = Pe/(R*T06)

'''Cálculo das vazões mássicas e área de saída dos bocais'''

a = np.array([[5, -1], [1, 1]])
b = np.array([0, 115])
x = np.linalg.solve(a, b)

mac = x[1]
mah = x[0]

Aef = mac/(rhof*Uef)
Ae  = mah/(rho*Ue)
#print('Ae',Ae)
#print('Aef',Aef)


'''Cálculo do impulso'''

fnh = (mah*((1+f)*(Ue-Uo)))+((Pe-Po)*Ae)
fnc = (B*mac*(Uef-Uo))+((Pef-Po)*Aef)
FN  = 2*(fnc+fnh)

#print('fnh',fnh)
#print('fnc',fnc)
#print('FN',FN)

'''
plt.figure(1)
plt.xlabel('Mo')
plt.ylabel('FN')
plt.plot(Mo, FN, 'r')
plt.show()
'''

'''ACOPLAMENTO MOTOR-AERONAVE'''


Wcf = 2000                                      #Peso de combustível final desejado **
delta1 = Po/101325                              #Pressão adimensional
Wempty = 155500                                 #Peso vazio sem motor

Clmax = 1.8                                     #CLMAX Admitido
Aw = 427.8                                      #Área das asas

rhocomb = 0.81
qtdcomb = 25000                                 #Litros
Wmotor = 7264                                   #Peso motor P&W PW4098
Wpas = 368*70                                   #peso total de passageiros com peso médio de 70kg 3 classes - config padrão
Wci = rhocomb*qtdcomb                           #Peso combustível possível
                                                #Densidade combustível (JET A)
#qtdcomb = Wci/(rhocomb*9.81)                   #Quantidade de combustível máxima em litros

print('Wci',Wci)
print('Quantidade de combustível', qtdcomb)

W   = (155500+2*(Wmotor)+Wci+Wpas)*9.81         #Peso total inicial
print('W', W)

Cl = (2*W) / (delta1*101325*(Mo**2)*Aw)         #coeficiente de sustentação
Cd = 0.026 + 0.0431*((Cl)**2)                   #Polar de arrasto


FNdelta1 = FN/delta1

WW = np.linspace(Wci, Wcf, 10)                      #Variação do peso inicial até o peso final step by step
newcl = (2*WW) / (delta1*gamma*101325*(Mo**2)*Aw)
cd2 = 0.026 + 0.0431*((newcl)**2)                   #Polar de arrasto
FNdelta2 = (gamma/2) * 101325 * cd2 * (Mo**2) * Aw
#print('Variação do peso inicial até o final', WW)
'''
plt.figure("Impulso Disponível / Requirido")
plt.plot(Mo, FNdelta1, 'r', label='Impulso disponível')
plt.plot(Mo, FNdelta2, 'b', label='Impulso requirido')
plt.legend()
plt.xlabel('Mo')
plt.ylabel('FN')
plt.show()
'''

rhoarSL = (101325/(R*298))                      #Densidade do ar Sea level
Vstall = (2*W/(Clmax*rhoarSL*Aw))**(1/2)        #Velocidade de stall
Vlo = (1.2)*Vstall                              #Velocidade de liftoff simplificada
FNDEC = 0.707*Vlo                               #FN Decolagem simplificada
amed = ((FNDEC-0.02*W) / W) * 9.81              #Aceleração média
machfnedc = FNDEC/347
SG = (Vlo**2)/(2*amed)
print('vstall', Vstall)
print('FNDEC', FNDEC)
print('amed', amed)
print('SG', SG)

'''Consumo e Tempo para decolagem'''
#ma - Calcular
ma = 115
mf = f[1] * ma
TSFC = mf/FN[1]
tdec = Vlo / amed                               #Tempo de descida
consumodec = TSFC*FN[1]*tdec                    #Consumo na descida
print('TSFC', TSFC)
print('Tempo de descida', tdec)
print('Consumo descida',consumodec)

'''Atualizando dados de peso / Subida e Descida'''

W2 = W-consumodec
D = Cd * 0.5*rhoarSL*Vlo*Aw
dhdt = ((FN[1]-D[1])*Vlo)/W
tsubida = (11000/dhdt)                          #Tempo de subida
consumosubida = TSFC*FN[1]*tsubida              #Consumo na subida

tdescida = tsubida                              #Porque a velocidade é constante
#consumodescida = consumosubida

print('dhdt',dhdt)
print('Consumo Subida', consumosubida)
print('Tempo de subida',tsubida)

'''Voo em cruzeiro'''

W3 = Wci - (consumodec+consumosubida)
rhocruz = Po/(R*To)
Vcruz = Mo*(gamma*R*To)**(1/2)
S = ((2/(9.81*TSFC)) * ((2/(Aw*rhocruz))**(1/2)) * (Cl[1] / Cd[1]**2) * ((W2**(1/2)) - (W3**(1/2))))
tcruzeiro = S / Vcruz
consumocruz = TSFC*FN[1] * tcruzeiro

