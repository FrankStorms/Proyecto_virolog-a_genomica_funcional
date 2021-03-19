#15/03/2021
#analizaremos las proteinas ICAM-1 que estan en NCBI, incluyen las que son solo predicciónes. 
#Previamente se ha realizado un formato fasta con los aminoacidos correspondientes de 148 secuencias, 
#(por errores de dedo ya no obtuvimos las 152 secuencias que estan en NCBI al día de hoy)
library(Biostrings)
library(BiocGenerics)
library(parallel)
library(msa)
library(seqinr)
library(ape)
###Primero realizar alineación####
ICAM_cw <- readAAStringSet ("6to semestre/Genomica funcional/ICAM1_VIRO/ICAM_1.fasta")
ICAM_cw

#Usamos los tres metodos con los que cuneta la librerya, y que de manera general, son los que comunmente se usan
ICAM_CW_R<- msa(ICAM_cw, "ClustalW")
ICAM_CW_R
ICAM_CO_R<- msa(ICAM_cw, "ClustalOmega")
ICAM_CO_R
ICAM_M_R<- msa(ICAM_cw, "Muscle")
ICAM_M_R

#Pasamos los objetos a un formato en el que se puedan leer para las funciones de otra librerya
ICAM_Aln_CW <- msaConvert(ICAM_CW_R, type="seqinr::alignment")
ICAM_Aln_CO <- msaConvert(ICAM_CO_R, type="seqinr::alignment")
ICAM_Aln_M <- msaConvert(ICAM_M_R, type="seqinr::alignment")


#De esta manera generamos las matrices de distancia
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

###Árbol filogenetico###
CW_Tree <- njs(d1)
plot(CW_Tree, main="Arbol Filogenetico de ICAM-1, Clusta W")

CW_Tree <- njs(d2)
plot(CW_Tree, main="Arbol Filogenetico de ICAM-1, Clusta Omega")

CW_Tree <- njs(d3)
plot(CW_Tree, main="Arbol Filogenetico de ICAM-1, Muscle")

# Al final generamos los tres árboles, sin embargo, al visualizarlos se veian los tres practicamente igual, y decidi optar por el de Clustal W. 
