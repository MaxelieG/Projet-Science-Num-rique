# Projet-Science-Num-rique
## Conduction instationnaire en 2D avec conditions de Dirichlet

Description du projet

Ce projet implémente une résolution numérique de l'équation de la chaleur en régime instationnaire dans un domaine rectangulaire 2D avec conditions aux limites de type Dirichlet. La résolution numérique utilise un schéma implicite et est comparée à la solution analytique afin d'évaluer la précision.

Le code est subdivisé en plusieurs fichiers. Le principal, \textit{main.py} permet d'obtenir la comparaison entre les solutions analytique et numérique ainsi que l'affichage et l'enregistrement de l'animation représentant le régime transitoire de la solution numérique. Le second code qu'il est possible d'exécuter est \textit{NumericalError.py}, il réutilise les mêmes fonctions que \textit{main.py} mais permet de traiter la dernier question du problème, c'est à dire l'approche de la norme $L_2$ de l'erreur par une formule du type $Ch^\alpha$. Les autres fichiers pythons sont utilisés par les deux principaux et leurs fonctionnement seront expliqués au fur et à mesure du rapport. Le fichier \textit{evolution\_temperature.mp4} est la dernière animation du régime transitoire de la solution numérique sauvegardée.