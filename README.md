# Projet-Science-Num-rique
## Conduction instationnaire en 2D avec conditions de Dirichlet

Description du projet

Ce projet implémente une résolution numérique de l'équation de la chaleur en régime instationnaire dans un domaine rectangulaire 2D avec conditions aux limites de type Dirichlet. La résolution numérique utilise un schéma implicite et est comparée à la solution analytique afin d'évaluer la précision.

Contexte physique

On considère un domaine rectangulaire :
-Dimensions : $L_x \times L_y$
-Conditions aux limites :
-Température imposée $T_1$ sur la face $y=0$.
-Température imposée $T_2$ sur les trois autres faces ($x=0$, $x=L_x$, $y=L_y$).

Méthodes utilisées

L'équation instationnaire de la chaleur est une parabolique.

1. Solution analytique


2. Solution numérique

Le programme utilise une méthode implicite pure avec la technique des directions alternées (ADI), pour résoudre l'équation de la chaleur où $\alpha_x = \frac{a \Delta t}{\Delta x^2}$ et $\alpha_y = \frac{a \Delta t}{\Delta y^2}$.

3. Comparaison et visualisation

Solution numérique obtenue par itération temporelle.
Solution analytique calculée directement pour comparaison.
Différence entre les deux solutions visualisée sous forme de carte thermique.



