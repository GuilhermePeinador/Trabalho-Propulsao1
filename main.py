# -*- coding: utf-8 -*-
"""
@author: Guilhermee Peinador Gomes
"""
import numpy as np
#import sympy as smp
import matplotlib.pyplot as plt

''''Motor escolhido P&W PW4098'''

'''DADOS PARA O CÁLCULO DO TURBOFAN'''

K = 273.15                      #Temperatura kelvin
Mo = 0.84                       #Mach de voo variando de 0.1 a 1 0.84 -> Cruzeiro np.linspace(0.1,1.2,12)
pif = 1.7                       #Taxa de compressão LPC Motor escolhido = 1.7 - 1.8
pic = 34.2                      #Taxa de compressão HPC Motor Escolhido = 34.2 - 42.8
T04 = 1944                      #Temperatura máxima câmara combustão
PCI = 42.8E6                    #Poder calorifico do combustível
cp = 1004
gamma = 1.4
Po = 22700                      #Pressão na entrada **** Colocar no relatório
To = -56.5 + K                  #Temperatura em K na entrada
B = 5.8                         #Bypass Motor Escolhido = 5.8 - 6.4
R = 287                         #Constante R
PSL= 101325                     #Pressão em nível do mar

print('Mach',Mo)
'''Intake'''

P01 = Po*(1+((gamma-1)/2)*Mo**2)**(gamma/(gamma-1))
T01 = To*(1+((gamma-1)/2)*Mo**2)
print("--------- Intake ---------")
print("P01", P01)
print("T01", T01)

'''Fan LPC'''

P02f = P01
T02f = T01
print("--------- Fan LPC---------")
print("P02f", P02f)
print("T02f", T02f)

P03f = pif*P02f
tauf = pif**((gamma-1)/gamma)
T03f = tauf*T02f
print("P03f", P03f)
print("tauf",tauf)
print("T03f",T03f)

'''Compressor HPC'''

P03 = P03f*pic
tauc = pic**((gamma-1)/gamma)
T03 = T03f*tauc
print("--------- Compressor HPC ---------")
print("P03", P03)
print("tauc",tauc)
print("T03",T03)

'''Câmara de combustão'''

P04 = P03
f = ((T04/T03)-1)/((PCI/(cp*T03))-(T04/T03))
print("--------- Câmara de Combustão ---------")
print("P04", P04)
print('f', f)

'''Conectar HPT-HPC'''

T05 = T04-((T03-T03f)/(1+f))
P05 = P04*((T05/T04)**((gamma)/(gamma-1)))
print("--------- Conexão HPT-HPC ---------")
print("P05", P05)
print("T05",T05)

'''Conectar LPT-LPC'''

T06 = T05-((1+B)*(T03f-T02f))/(1+f)
P06 = P05*(T06/T05)**((gamma)/(gamma-1))
print("--------- Conexão LPT-LPC ---------")
print("P06", P06)
print("T06",T06)

'''Bocal 3' - 7' '''

P07f = P03f
NPRf = P07f/Po
Pef = P07f/1.892
print("--------- Bocal 3'-7' ---------")
print("P07f", P07f)
print("NPRf",NPRf)
print("Pef", Pef)

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
print("--------- Bocal 6-7' ---------")
print("P07", P07)
print("NPR",NPR)
print("Pe", Pe)
print("Tef",Tef)
print("Te", Te)
print("Uef",Uef)
print("Ue", Ue)
print("Uo",Uo)
print("rhof", rhof)
print("rho",rho)

'''Cálculo das vazões mássicas e área de saída dos bocais'''

ma = rho*Uo*np.pi*1.4**2         #Fluxo de masssa de ar que entra no motor
mac = (B*ma)/(B+1)                  #Fluxo 'frio'
mah = mac/B                         #Fluxo'quente'

Aef = mac/(rhof*Uef)                #Área do bocal secundário
Ae  = mah/(rho*Ue)                  #Área do bocal primário
print("--------- Vazões mássicas e áreas de saída ---------")
print("ma", ma)
print("mac", mac)
print("mah", mah)
print('Ae',Ae)
print('Aef',Aef)


'''Cálculo do impulso'''

fnh = mah*((1+f)*Ue-Uo)+((Pe-Po)*Ae)    #Impulso 'quente'
fnc = B*mac*(Uef-Uo)+(Pef-Po)*Aef       #Impulso 'frio'
FN  = 2*(fnc + fnh)                     #impulso total
print("--------- Cálculo dos Impulsos ---------")
print('fnh',fnh)
print('fnc',fnc)
print('FN total (2 motores)',FN)

'''ACOPLAMENTO MOTOR-AERONAVE'''

delta1 = Po/101325                              #Pressão adimensional
Wempty = 155500                                 #Peso vazio sem motor
Wmotor = 44452.052                              #Peso motor PW4098

Clmax  = 1.8                                    #CLMAX Admitido
Aw     = 427.8                                  #Área das asas

rhocomb = 0.81                                  #Densidade combustível (JET A)
qtdcomb = 50000                                 #Litros

Wpas = 386*70                                   #Peso total de passageiros com peso médio de 70kg 3 classes - config padrão + bagagem max
Wci = rhocomb*qtdcomb                           #Peso combustível inicial
Wcf = 5000                                      #Peso combustível final
print("--------- ACOPLAMENTO MOTOR-AERONAVE ---------")
print('Quantidade de combustível em litros', qtdcomb)
print('Wci',Wci)
print('Peso de passageiros', Wpas)

W = 155500+2*Wmotor+Wci+Wpas                    #Peso total inicial
print('Peso total da aeronave (W) na decolagem', W)     #Limite = 2935816.8105N

Cl = (2*W)/(delta1*101325*(Mo**2)*Aw)                   #Coeficiente de sustentação
Cd = 0.026 + 0.0431*((Cl)**2)                           #Polar de arrasto

FNdelta1 = FN/delta1
print('FNdelta1', FNdelta1)

'''DECOLAGEM Mo[1]'''
rhoarSL = 1.2                                   #Densidade do ar Sea level
Vstall = (2*W/(Clmax*rhoarSL*Aw))**(1/2)        #Velocidade de stall
Vlo = 1.2*Vstall                                #Velocidade de liftoff simplificada
Pista = 3700                                    #Tamanho da pista em metros  #altera pista pra FNdec aproximar do FN
amed = Vlo**2/(2*Pista)
FNdec = amed*W+0.02*W
print("--------- DECOLAGEM ---------")
print('vstall', Vstall)
print('Vlo', Vlo)
print('FNdec', FNdec)
print('amed', amed)
print('Tamanho da pista', Pista)

WW = np.linspace(Wci, Wcf, 12)                      #Variação do peso inicial até o peso final step by step
newcl = (2*WW) / (delta1*gamma*101325*(Mo**2)*Aw)
cd2 = 0.026 + 0.0431*((newcl)**2)                   #Polar de arrasto
FNdelta2 = (gamma/2) * 101325 * cd2 * (Mo**2) * Aw
#print('Variação do peso inicial até o final', WW)




'''PLOT IMPULSOS'''
plt.figure("Impulso Disponível / Requirido")
plt.plot(Mo, FNdelta1, 'r', label='Impulso disponível')
plt.plot(Mo, FNdelta2, 'b', label='Impulso requirido')
plt.legend()
plt.xlabel('Mo')
plt.ylabel('FN')
#plt.show()


'''CONSUMO E TEMPO DE DESCIDA'''

mf = f * ma
TSFC = mf/FN
tdec = Vlo/amed                                     #Tempo de descida
consumodec = TSFC*FN*tdec                           #Consumo na descida
print('mf',mf)
print("--------- TSFC ---------")
print('TSFC', TSFC)
print("--------- DESCIDA ---------")
print('Tempo de descida', tdec)
print('Consumo descida',consumodec)

'''ATUALIZANDO DADOS DE PESO'''

W2 = W - consumodec
D = Cd * 0.5*rhoarSL*Vlo*Aw
dhdt = ((FN-D)*Vlo)/W
tsubida = (11000/dhdt)                          #Tempo de subida
consumosubida = TSFC*FN*tsubida                 #Consumo na subida

tdescida = tsubida                              #Porque a velocidade é constante

consumodescida = consumosubida                  #Aproximação apra o trabalho

print("--------- ATUALIZAÇÃO DOS DADOS DE PESO ---------")
print('W2',W2)
print('dhdt',dhdt)
print('Consumo Subida', consumosubida)
print('Tempo de subida',tsubida)
print('Tempo de descida',tdescida)

'''VOO DE CRUZEIRO Mo[8]'''

W3 = W2 - (consumodec + consumosubida)
S = 8000000                                     #Distância em metros da missão
Vcruz = 288.12                                  #Velocidade em metros por segundo
tcruzeiro = S / Vcruz
consumocruz = TSFC*FN * tcruzeiro
print("--------- VOO DE CRUZEIRO ---------")
print('W3',W3)
print('Velocidade de cruzeiro',Mo)
print('S',S)
print('Tempo em cruzeiro',tcruzeiro/3600)
print('Consumo em voo de cruzeiro',consumocruz)