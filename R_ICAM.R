#15/03/2021
#analizaremos las proteinas ICAM-1 que estan en NCBI, incluyen las que son solo predicciónes. 
#Previamente se ha realizado un formato fasta con los aminoacidos correspondientes de 151 secuencias
library(Biostrings)
library(BiocGenerics)
library(parallel)
library(msa)
library(seqinr)
library(ape)
###Primero realizaremos un árbol filogenetico####
ICAM_cw <- readAAStringSet ("6to semestre/Genomica funcional/ICAM1_VIRO/ICAM_1.fasta")
ICAM_cw

ICAM_CW_R<- msa(ICAM_cw, "ClustalW")
ICAM_CW_R
ICAM_CO_R<- msa(ICAM_cw, "ClustalOmega")
ICAM_CO_R
ICAM_M_R<- msa(ICAM_cw, "Muscle")
ICAM_M_R

ICAM_Aln_CW <- msaConvert(ICAM_CW_R, type="seqinr::alignment")
ICAM_Aln_CO <- msaConvert(ICAM_CO_R, type="seqinr::alignment")
ICAM_Aln_M <- msaConvert(ICAM_M_R, type="seqinr::alignment")


d1 <- dist.alignment(ICAM_Aln_CW, "identity")
Matr_homo <-as.matrix(d1)[2:150, "sp|P05362.2|Homo sapiens|ICAM1_HUMAN RecName: Full=Intercellular adhesion molecule 1; Short=ICAM-1; AltName: Full=Major group rhinovirus receptor; AltName: CD_antigen=CD54; Flags: Precursor", drop=FALSE]
View (Matr_homo)
Matr_rino <-as.matrix(d1)[1:149, "pdb|1D3L|rinovirus|A Chain A, PROTEIN (INTERCELLULAR ADHESION MOLECULE-1)", drop=FALSE]
View (Matr_rino)

d2 <- dist.alignment(ICAM_Aln_CW, "identity")
Matr_homo2 <-as.matrix(d2)[2:150, "sp|P05362.2|Homo sapiens|ICAM1_HUMAN RecName: Full=Intercellular adhesion molecule 1; Short=ICAM-1; AltName: Full=Major group rhinovirus receptor; AltName: CD_antigen=CD54; Flags: Precursor", drop=FALSE]
View(Matr_homo2)
Matr_rino2 <-as.matrix(d2)[1:149, "pdb|1D3L|rinovirus|A Chain A, PROTEIN (INTERCELLULAR ADHESION MOLECULE-1)", drop=FALSE]
View(Matr_rino2)

d3 <- dist.alignment(ICAM_Aln_CW, "identity")
Matr_homo3 <-as.matrix(d3)[2:150, "sp|P05362.2|Homo sapiens|ICAM1_HUMAN RecName: Full=Intercellular adhesion molecule 1; Short=ICAM-1; AltName: Full=Major group rhinovirus receptor; AltName: CD_antigen=CD54; Flags: Precursor", drop=FALSE]
View(Matr_homo3)
Matr_rino3 <-as.matrix(d3)[1:149, "pdb|1D3L|rinovirus|A Chain A, PROTEIN (INTERCELLULAR ADHESION MOLECULE-1)", drop=FALSE]
View(Matr_rino3)

#árbol filogenetico
CW_Tree <- njs(d1)
plot(CW_Tree, main="Arbol Filogenetico de ICAM-1, Clusta W")

CW_Tree <- njs(d2)
plot(CW_Tree, main="Arbol Filogenetico de ICAM-1, Clusta Omega")

CW_Tree <- njs(d3)
plot(CW_Tree, main="Arbol Filogenetico de ICAM-1, Muscle")


#Ciclo para elegír a las temperaturas pertinentes

#1) Inicializar los nombres y valores de los poderes :
poderes <- c("Volar", "Dinero", "Invisibilidad", "Telepatia",
             "Telara?a", "Hipnosis", "Vision-laser", "Super-fuerza",
             "Teletransportacion", "Sarcasmo")
grado <- c(23, 45, 66, 19, 80, 3, 24, 90, 1, 78)

#2) Inicializar los nombres de los superheroes :
heroe <- c("Falcon", "Batman", "Mujer-invisible",
           "Jean", "Kaliman", "Spiderman", "Ciclope",
           "Superman", "Nightcrawler", "X" )

#3) Asociar los valores con los nombres :
names(grado) <- poderes

#4) Inicializar una variable umbral para la frase del Tio Ben :


umbral <- .5 #arbitrario pero menor al valor de telara?a

#5) Indicar el loop :
for (i in 1:length(Matr_rino)) {
  if(grado[i] < umbral) {
    print(paste(Matr_rino [i], #para ver el nombre y no solo el index en el vector "grado"
                ", un gran poder conlleva una gran responsabilidad... Atte. Ben Parker"))
  }
}
