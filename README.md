# CHP2

4/12/2019
Code : découpage du domaine
Problème : sans recouvrement, on obtient pas deux solution (ce que l'on devrait avoir)
Commande de compilation:
mpif90 -fbounds-check fonctions.f90 gc.f90 ma.f90 parametres.f90 pg1.f90 -o pg1
mpirun -n 2 pg1
