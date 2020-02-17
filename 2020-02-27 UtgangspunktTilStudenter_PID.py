# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 17:34:29 2019

@author: cakir
"""

#IMPORTERER BIBLIOTEK
import numpy as np                             #Trenger dette biblioteket til sinus, pi, opprette array og matriseregning
import matplotlib.pyplot as mat                #Trenger dette biblioteket til plotting

#SETTER VERDIER
A_v= 30                                         #Areal innervegg
A_f=5                                           #Areal vindu
A_y=10                                          #Areal yttervegg
tykkelse_vegg = 0.01                            #Tykkelse vegg
V= A_v * tykkelse_vegg                          #Regner ut volum vegg
U_y = 0.3                                       #U-verdi yttervegg
V_til=100                                       #luftmengde tilluft m3/h   
inf=10*2.7*0.2                                  #Luftmengde infiltrasjon høyde 27m n=0.2
C_p= 960                                        #Spesifikk varmekapasitet betong
rho=2200                                        #tetthet betong
lamda_delta = 0.2/0.1                           #konduktiv ledningsevne inn til veggen
alpha_k=3.5                                     #konveksjon
alpha_s=4.5                                     #stråling
U_f=0.8                                         #U-verdi vindu
u_f=(1)/((1/U_f)-(1/(alpha_k+alpha_s)))         #u-verdi innerglass
sol_trans= 0.6                                  #Transmisjon sol
lys_trans= 0.5                                  #Transmisjon lys
sol_abs=0.4                                     #sol absorbsjon
t_s=3600                                        #tidsskritt 1t
T_til=18                                        #Temperatur tilluft
P_person=0                                      #Effekt fra personer
P_data=0                                        #Effekt fra data
P_lys=0                                         #Effekt fra lys
P_varme=0                                       #Effekt fra varme - initiell verdi, skal regnes ut senere
P_kjol=0                                        #Effekt fra kjøling initiell verdi, skal regnes ut senere
    
masse_init=20                                   #Massetemperatur ved start
rom_init=20
T_onsket=21                                     #Ønsket temperatur i rommet

T_ute=20                                        #Utetemperatur (konstant variant) - lages kurve for senere
P_sol=0                                         #Effekt fra sol (konstant variant) - lages kurve for senere

n=120                    #Antall tidsskritt i en time. 3600/120 = 30s per tidsskritt


#OPPRETTER MATRISER
A=np.zeros(shape= (4,4), dtype='float')                     #Oppretter 4x4 matrise som en "array" med variabeltype "float" = desimaltall - A matrise som skal fylles inn
B=np.zeros(shape= (4,9), dtype='float')                     #Oppretter 4x9 matrise som en "array" med variabeltype "float" = desimaltall - B matrise som skal fylles inn

#FYLLER INN VERDIER I A OG B MATRISEN
#Temperatur ROM - ledd
A[0,1]=A_v*alpha_k                                          #Trom-Tvegg: konveksjon ganger areal innervegg
A[0,2]=A_f*alpha_k                                          #Trom-Tvindu: konveksjon ganger areal vindu
B[0,0]=V_til*0.335                                          #Trom-Ttilluft: Luftmengde tilluft ganger 0.335 (W per luftmengde K)
B[0,1]=inf*0.335 + U_y*A_y                                  #Trom-Tute: U-verdi yttervegg ganger areal yttervegg + infiltrasjon ganger 0.335
B[0,3]=1                                                    #Trom-Ppers: 1. All varme til rommet. 
B[0,4]=0.5                                                  #Trom-Plys: konveksjon 0.5. (Lys 0,5 konveksjon og 0,5 transmisjon)
B[0,5]=1                                                    #Trom-Pdata: 1. All varme til rommet. 
B[0,6]=1                                                    #Trom-Pvarme: 1. All varme til rommet
B[0,7]=1                                                    #Trom-Pkjøling: 1. All kjøling til rommet
A[0,0]=-A[0,1]-A[0,2]-A[0,3]-B[0,0]-B[0,1]-B[0,2]           #Trom-Trom: summen av alle temperaturledd. (A:0 til 3, B 0 til 2)

#Temperatur VEGG - ledd
A[1,0]=A[0,1]                                               #Tvegg-Trom = Trom-Tvegg (konveksjon ganger areal innervegg)
A[1,2]=(A_v)*(A_f/(A_v+A_f))*alpha_s                        #Tvegg-Tvindu: formfaktor (hvor mye veggen "ser" av vinduet, arealforhold) ganger stråling
A[1,3]=A_v*lamda_delta                                      #Tvegg-Tmasse: Areal av innervegg ganger lambda_delta - konduktiv ledningsevne inn til massesenter
B[1,4]=0.5                                                  #Tvegg-Plys: transmisjon 0.5
B[1,8]=sol_trans                                            #Tvegg-Psol: sol transmisjon
A[1,1]=-A[1,0]-A[1,2]-A[1,3]-B[1,0]-B[1,1]-B[1,2]           #Tvegg-Tvegg: summen av alle temperaturledd

#Temperatur Vindu - ledd
A[2,0]=A[0,2]                                               #Tvindu-Trom = Trom-Tvindu (konveksjon ganger areal vindu)
A[2,1]=A[1,2]                                               #Tvindu-Tvegg = Tvegg-Tvindu (Areal vindu ganger formfaktor vindu-vegg ganger stråling)
B[2,1]=u_f*A_f                                              #Tvindu-Tute: Areal vindu ganger u-verdi innervindu
B[2,8]=sol_abs                                              #T-vindu-Psol: Solabsorbsjon
A[2,2]=-A[2,0]-A[2,1]-A[2,3]-B[2,0]-B[2,1]-B[2,2]           #Tvindu-Tvindu: summen av alle temperaturledd

#Temperatur masse - leddfile:///J:/30_Tekniske_systemer/310_Industriell_Elektro_Aut/BIM/AoE/%23C/Python/OsloMet/2020/2020-01-23 Fra Forelesning_Eksaktregulering.py
A[3,1]=A[1,3]                                               #Tmasse-Tvegg = Tvegg-Tmasse (Areal vindu ganger lambda-ro, inn til massesenter)                
B[3,2]=V*C_p*rho/(t_s/n)                                                    #Tmasse-Tmasse_forrigesteg: Volum vegg ganger varmekapasitet betong ganger tetthet betong delt på tidsskritt
A[3,3]=-A[3,0]-A[3,1]-A[3,2]-B[3,0]-B[3,1]-B[3,2]           #Tmasse-Tmasse: summen av alle temperaturledd                                    

#FINNER R MATRISEN
#Regner ut R=A-1 * (-1)B
A_inv=np.linalg.inv(A)              #Inverterer A ved hjelp av np-funksjon
B=(-1)*B                            #Ganger B med (-1). Siden B er en array kan en utføre dette direkte, -1 ganges da med alle tall i B
R=np.dot(A_inv,B)                   #R er dot-produktet mellom invertert A og negativ B


#Funksjon for å lage utetemperaturer etter en sinuskurve
def utetemperatur(T_min, T_maks, Klokkeslett):
    T_middel=(T_maks+T_min)/2
    Amplitude=(T_maks-T_min)/2
    dklokke=4
    utetemp=(np.sin(((Klokkeslett-dklokke)/24)*2*np.pi))*Amplitude + T_middel
    return(utetemp)

#Lager utetemperaturkurve
T_ute=[0]*24  
for i in range(24):      
    T_ute[i]=utetemperatur(0,20,i)      #For alle klokkeslett fra 0 til 23 kalles utetemperatur-funksjonen slik at vi får 24 temperaturer
 

#Funksjon for å lage solstråling etter sinuskurve
def solstraaling(S_min, S_maks, Klokkeslett):
    S_middel=(S_maks+S_min)/2
    Amplitude=(S_maks-S_min)/2
    dklokke=4
    solstr=(np.sin(((Klokkeslett-dklokke)/24)*2*np.pi))*Amplitude + S_middel
    return(solstr)

#Lager kurve over solstråling    
P_sol=[0]*24 
for i in range(24):                     #For alle klokkeslett fra 0 til 23 kalles solstråling-funksjonen slik at vi får 24 effekter
    P_sol[i]=solstraaling(0,900,i)
        
P_sol_faktor=2                          #Gir mulighet for å gange opp solstrålingen med en fakor
for i in range(len(P_sol)):
    P_sol[i]*=P_sol_faktor

#Funksjon for plotting av temperaturer
def rommodell_plot_temp(T_ute, T_rom, T_vegg, T_masse):     #Lager en funksjon som plotter Tute, Trom og Tmasse, må da ta disse tabellene inn som input
    t=np.linspace(1,24,24)                          #Oppretter en tabell med tidene, fra 1 til 24
    mat.plot(t, T_ute)                              #T_ute kan da plottes som funksjon av t
    mat.plot(t, T_vegg)
    mat.plot(t, T_rom)                              #T_rom kan da plottes som funksjon av t
    mat.plot(t, T_masse)                            #T_masse plottes som funksjon av t
    mat.xlabel("Tid [h] ")
    mat.ylabel("Temperatur [C]")
    mat.legend(["Utetemperatur", "Temperatur vegg", "Romtemperatur", "Massetemperatur"])
    mat.title("Temperaturkurver")
    mat.show()                                         #Kaller plotte-temperaturfunksjonen for resultatene med regulering

#Funksjon for å plotte effekter
def rommodell_plot_effekt(P_sol, P_varme, P_kjol):                                              #Lager funksjon for å plotte effekter
    t=np.linspace(1,24,24)
    mat.plot(t, P_sol)
    mat.plot(t, P_varme)
    mat.plot(t, P_kjol)
    mat.xlabel("Tid [h] ")
    mat.ylabel("Effekt [W]")
    mat.legend(["Solstråling", "Varme", "Kjøling"])
    mat.title("Effektkurver")
    mat.show()
    
    ##START PID###
    
#P-vektor
P=[0]*9
P[0]=T_til
P[3]=P_person
P[4]=P_lys
P[5]=P_data

Ts=t_s/n

def PID_Innregulering(Kp, Ti, Td):
    T_rom=[0]*24
    T_rom[0]=rom_init
    T_masse=[0]*24
    T_masse[0]=masse_init
    
    deltaT=[0]*24
    I=[0]
    
    P[1]=-10
    P[8]=0
    
    for i in range(1,24):
        P[2]=T_masse[i-1]
        deltaT[i]=T_onsket-T_rom[i-1]
        u=Kp*deltaT[i] + I[i-1] + (Kp*Ts/Ti)*deltaT[i] + (Kp*Td/Ts)*(deltaT[i]-deltaT[i-1])
        I.append(I[i-1] + (Kp*Ts/Ti)*deltaT[i])
        if(u>0):
            P[6]=u
            P[7]=0
        else:
            P[6]=0
            P[7]=u
        T=np.dot(R,P)
        T_rom[i]=T[0]
        T_masse[i]=T[3]
    t=np.linspace(1,24,24)
    mat.plot(t, T_rom)
    mat.show()

Kp=80
Ti=10000
Td=0
PID_Innregulering(Kp, Ti, Td)

Kpu=80

def P_kontroller(Kpu):
    Kp=0.5*Kpu
    Ti=100000
    Td=0
    return Kp, Ti, Td

def PI_kontroller(Kpu):
    Kp=0.45*Kpu
    Ti=2*Ts/1.2
    Td=0
    return Kp, Ti, Td

def PID_kontroller(Kpu):
    Kp=0.6*Kpu
    Ti=2*Ts/2
    Td=Ti/4
    return Kp, Ti, Td

P_Kp, P_Ti, P_Td=P_kontroller(Kpu)
PID_Innregulering(P_Kp, P_Ti, P_Td)

PI_Kp, PI_Ti, PI_Td=PI_kontroller(Kpu)
PID_Innregulering(PI_Kp, PI_Ti, PI_Td)

PID_Kp, PID_Ti, PID_Td=PID_kontroller(Kpu)
PID_Innregulering(PID_Kp, PID_Ti, PID_Td)

