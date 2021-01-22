# Moskvina V, Harold D, Russo G, Vedernikov A, Sharma M, Saad M, Holmans P, Bras JM, Bettella F, Keller MF, Nicolaou N, Simón-Sánchez J, Gibbs JR, Schulte C, Durr A, Guerreiro R, Hernandez D, Brice A, Stefansson H, Majamaa K, Gasser T, Heutink P, Wood N, Martinez M, Singleton AB, Nalls MA, Hardy J, Owen MJ, O'Donovan MC, Williams J, Morris HR, Williams NM, IPDGC and GERAD Investigators. Analysis of genome-wide association studies of Alzheimer disease and of Parkinson disease to determine if these 2 diseases share a common genetic risk. JAMA Neurol. 2013 Oct;70(10):1268–76. 

# study 1(AD) case/ control samples
S1.case = 3177
S1.cont = 7277

# study 2(PD) case/ control samples
S2.case = 5333
S2.cont = 12298

# study 1-2 shared samples (case/ control) 
shared.case = 0
shared.cont = 5571

# N_S1 N_S2
S1.tot = S1.case + S1.cont
S2.tot = S2.case + S2.cont

# estimate correlation due to overlapping samples
r_ij = ( shared.cont * sqrt( (S1.case * S2.case)  / (S1.cont * S2.cont) ) + shared.case * sqrt( (S1.cont*S2.cont) / (S1.case * S2.case) ) ) / sqrt(S1.tot * S2.tot)

## print out the result
print(r_ij)
