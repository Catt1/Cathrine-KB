# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 08:32:38 2019

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
T_til=20                                        #Temperatur tilluft
P_person=0                                      #Effekt fra personer
P_data=0                                        #Effekt fra data
P_lys=0                                         #Effekt fra lys
P_varme=0                                       #Effekt fra varme - initiell verdi, skal regnes ut senere
P_kjol=0                                        #Effekt fra kjøling initiell verdi, skal regnes ut senere
    
masse_init=20                                   #Massetemperatur ved start
T_onsket=25                                     #Ønsket temperatur i rommet

T_ute=20                                        #Utetemperatur (konstant variant) - lages kurve for senere
P_sol=0                                         #Effekt fra sol (konstant variant) - lages kurve for senere

"""
print("u_f: ")
print(u_f)                                     #tester koden fram til nå, sjekker at u_f blir 0.8888
print("\n")
"""

#OPPRETTER MATRISER
A=np.zeros(shape= (4,4), dtype='float')                     #Oppretter 4x4 matrise som en "array" med variabeltype "float" = desimaltall - A matrise som skal fylles inn
B=np.zeros(shape= (4,9), dtype='float')                     #Oppretter 4x9 matrise som en "array" med variabeltype "float" = desimaltall - B matrise som skal fylles inn

"""                                                         #Skriver ut A og B for å se 0-matrisene
print("A før verdier er fylt inn: ")
print(A)
print("\n")
print("B før verdier er fylt inn: ")
print(B)
print("\n")
"""
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

#Temperatur masse - ledd
A[3,1]=A[1,3]                                               #Tmasse-Tvegg = Tvegg-Tmasse (Areal vindu ganger lambda-ro, inn til massesenter)                
B[3,2]=V*C_p*rho/t_s                                        #Tmasse-Tmasse_forrigesteg: Volum vegg ganger varmekapasitet betong ganger tetthet betong delt på tidsskritt
A[3,3]=-A[3,0]-A[3,1]-A[3,2]-B[3,0]-B[3,1]-B[3,2]           #Tmasse-Tmasse: summen av alle temperaturledd                                    

"""
print("A: ")                                                #Skriver ut A og B matrisene for kontroll
print(A)
print("\n")
print("B:")
print(B)
print("\n")
"""

"""
print("Tester at likningene er i balanse. Sum av temp-ledd i A og B = 0 :")        #Summerer alle temperaturleddene for hver rad. Sjekker at de er lik 0. 
kontroll_A_B=[]
for i in range(4):                  #Dette skal gjøres for alle 4 rader
    sum=0
    for j in range(4):              #Temperaturledd i A er alle kolonnene: 0, 1, 2, 3. Range(4) gir 0, 1, 2, 3 
        sum+=A[i,j]    
    for j in range(3):              #Temperaturledd i B er de tre første kolonnene: 0, 1, 2. Range (3) gir 0, 1, 2
        sum+=B[i,j]  
    kontroll_A_B.append(sum)        #Legger til summen i en tabell som vi kan skrive ut, for å få resultatet av kontroll på alle 4 rader
print(kontroll_A_B)                 #Gir svar 0 eller e-15 tilsvarende lik 0 - OK
print("\n")
"""

#FINNER R MATRISEN
#Regner ut R=A-1 * (-1)B
A_inv=np.linalg.inv(A)              #Inverterer A ved hjelp av np-funksjon
B=(-1)*B                            #Ganger B med (-1). Siden B er en array kan en utføre dette direkte, -1 ganges da med alle tall i B
R=np.dot(A_inv,B)                   #R er dot-produktet mellom invertert A og negativ B

"""
print("R:")                         #Skriver ut R for kontroll
print(R)
print("\n")
"""

"""
print("Tester at R er riktig. Sum av alle templedd = 1 : ")         #Summerer alle temperaturleddene for de tre første kolonnene = temperaturledd, sjekker at sum er 1. 
kontroll=[]
for i in range(4):
    sum=0
    for j in range(3):
        sum+=R[i,j]
    kontroll.append(sum)
print(kontroll)
print("\n")
"""

#Oppretter P-vektor. Pådrag
P=[0]*9                             #Matrisen skal bestå av 9 elementer. Oppretter først en tabell av 9 0-er. Fyller deretter inn verdiene. 
P[0]=T_til                          #P[0]: Temperatur tilluft - konstant
P[1]=T_ute                          #P[1]: Temperatur ute - settes først lik konstant. Vil senere variere. 
P[2]=masse_init                     #P[2]: Massetemperatur ved start - må bruke en initiell verdi
P[3]=P_person                       #P[3]: Effekt fra personbelastning - konstant
P[4]=P_lys                          #P[4]: Effekt fra lys - konstant
P[5]=P_data                         #P[5]: Effekt fra data - konstant
P[6]=P_varme                        #P[6]: Effekt fra varmekilde - skal beregnes ved regulering. Settes først til 0 for å regne temperatur uten
P[7]=P_kjol                         #P[7]: Effekt fra kjøling - skal beregnes ved regulering. Settes først til 0 for å regne temperatur uten
P[8]=P_sol                          #P[8]: Effekt fra solen - settes først lik konstant. Vil senere variere. 

"""
print("P: ")                        #Skriver ut P for kontroll
print(P)
print("\n")
"""

"""
print("T: ")                        #Sjekker at T blir 20 grader for alle ledd ved initielle verdier
T=np.dot(R,P)
print(T)
print("\n")
"""

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
 
"""
print("Temperaturkurve: ")
print(T_ute)
print("\n")
"""

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

"""    
print("Solstrålingkurve: ")
print(P_sol)
print("\n")
"""
        
P_sol_faktor=2                          #Gir mulighet for å gange opp solstrålingen med en fakor
for i in range(len(P_sol)):
    P_sol[i]*=P_sol_faktor

"""
print("Solstrålingskurve etter ganging med solfaktor: ")
print(P_sol)
print("\n")
"""

#Funksjon for å beregne utetemperatur uten regulering
def rommodell_utenReg(T_ute, P_sol):                #Funksjonen tar inn en kurve av utetemperaturer og en kurve av solstråling
    T_rom=[0]*24                                    #Vi skal samle på temperatur i rom som et resultat, oppretter en tabell for disse verdiene
    T_vegg=[0]*24
    T_masse=[0]*24                                  #Vi skal samle på temperatur masse som et resultat, oppretter en tabell for disse verdiene
    P_varme=[0]*24                                  #Vi skal samle på varmeeffekt som et resultat, oppretter en tabell for disse verdiene
    P_kjol=[0]*24                                   #Vi skal samle på kjøleeffekt som et resultat, oppretter en tabell for disse verdiene
    T_masse_prev=[0]*25                             #Vi må vite hva massetemp i forrige ledd var. Denne vil fylles inn for neste gang. Trenger derfor 25 verdier i denne. 
    
    T_masse_prev[0]=masse_init                      #Setter masse fra forrige ledd lik initiell verdi for første kjøring
    
    for i in range(24):                             #Itererer over de 24 timene
        P[1]=T_ute[i]                               #For hver time hentes riktig verdi fra tabellen T_ute ut og legges inn i P vektor på plass P[1] 
        P[2]=T_masse_prev[i]                        #For hver time hentes riktig verdi fra tabellen T_masse_prev ut og legges inn i P vektor på plass P[2] 
        P[6]=0                                      #Først skal temperatur om ingen kjøling eller varme regnes ut, derfor settes varme P[6] til 0
        P[7]=0                                      #Først skal temperatur om ingen kjøling eller varme regnes ut, derfor settes kjøling P[7] til 0
        P[8]=P_sol[i]                               #For hver time hentes riktig verdi fra tabellen P_sol ut og legges inn i P vektor på plass P[8]  
        T=np.dot(R,P)                               #T-vektor regnes ut som dot produkt av R-matrise og P-vektor
        T_rom[i]=T[0]                               #T[0] er temperatur i rommet, denne skrives til T_rom-tabellen der vi samler opp resultater
        T_vegg[i]=T[1]
        T_masse[i]=T[3]                             #T[3] er temperatur masse, denne skrives til T_masse-tabellen der vi samler opp resultater
        T_masse_prev[i+1]=T[3]                      #T_masse_prev[i+1] er massetemper_forrige neste steg, denne settes lik massetemp nå
    return T_rom, T_vegg, P_varme, P_kjol, T_masse          #Funksjonen retunerer en tabell over romtemp, varme, kjøling og Tmasse

T_rom_ureg, T_vegg_ureg, P_varme_ureg, P_kjol_ureg, T_masse_ureg = rommodell_utenReg(T_ute, P_sol)       #Kall til funksjonen med T_ute og P_sol som vi har laget

"""
print("T rom uten regulering: ")                    #Skriver ut romtemperaturen for kontroll
print(T_rom_ureg)
print("\n")
"""

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
    mat.show()

print("Temperaturkurve uten regulering: ")
rommodell_plot_temp(T_ute, T_rom_ureg, T_vegg_ureg, T_masse_ureg)        #Kaller plottefuksjonen med tabell for utetemp som ble generert, og trom og tmasse beregnet uten regulering

#Funksjon for å beregne utetemperatur med regulering
def rommodell_medReg(T_ute, P_sol):         
    T_rom=[0]*24
    T_vegg=[0]*24
    T_masse=[0]*24
    P_varme=[0]*24
    P_kjol=[0]*24
    T_masse_prev=[0]*25
    
    T_masse_prev[0]=masse_init
    
    for i in range(24):
        P[1]=T_ute[i]
        P[2]=T_masse_prev[i]
        P[6]=0
        P[7]=0
        P[8]=P_sol[i]
        T=np.dot(R,P)
        delta_T=T_onsket-T[0]                               #Etter T er beregnet kan vi finne ut avviket: T_onsket-T i rommet for dette steget
        if delta_T>=0:                                      #Om avviket er positivt trenger vi varme
            varmebehov=delta_T/R[0,6]                       #Varmepådraget er lik delta T delt på R_rom_Varme = [R0,6]
            P_varme[i]=varmebehov                           #Fyller inn varmebehovet i resultattabellen
            P[6]=varmebehov                                 #Setter pådraget likt varmebehovet for neste utregning
        else:                                               #Om avviket er negativt trenger vi kjøling
            kjolebehov=delta_T/R[0,7]                       #Kjølebehovet er lik delta T delt på R_rom_Kjøling = [R0,7]
            P_kjol[i]=kjolebehov                            #Fyller inn kjølebehovet i resultattabellen
            P[7]=kjolebehov                                 #Setter pådraget likt kjølebehovet for neste utregning
        T=np.dot(R,P)                                       #Nå som P[6]-varme eller P[7]-kjøling er endret regner vi ut T på nytt
        T_rom[i]=T[0]                                       #Temperaturen i rommet T[0] blir nå den endelige verdien for steget, fylles inn i resultattabellen
        T_vegg[i]=T[1]
        T_masse[i]=T[3]                                         
        T_masse_prev[i+1]=T[3]
    return T_rom, T_vegg ,P_varme, P_kjol, T_masse

T_rom_mreg, T_vegg_mreg, P_varme_mreg, P_kjol_mreg, T_masse_mreg = rommodell_medReg(T_ute, P_sol)            #Kaller rommodellen med regulering for T_ute og Psol

"""
print("T rom med regulering: ")
print(T_rom_mreg)
print("\n")
"""

print("Temperaturkurve med regulering")
rommodell_plot_temp(T_ute, T_rom_mreg, T_vegg_mreg, T_masse_mreg)                                            #Kaller plotte-temperaturfunksjonen for resultatene med regulering

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

print("Effektkurve med regulering: ")    
rommodell_plot_effekt(P_sol, P_varme_mreg, P_kjol_mreg)                                     #Kaller plot-effekt-funksjonen for resultatene med regulering
        
    
    
    
    


