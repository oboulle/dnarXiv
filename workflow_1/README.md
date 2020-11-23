# Workflow_1

workflow utilisant les différents simulateurs

fonctionne sans interface graphique, l'utilisateur doit juste donner les paramètres demandés

le cycle fait : séquences -> ajout des primers -> synthèse -> séquençage -> basecalling -> demultiplexing -> consensus

necessite d'avoir les projets du Github en local

1) Utilisation interactive avec python

mettre les chemins des projets locaux dans le fichier project_paths.txt

commande : 

*python3 workflow_1.py*

2) Utilisation automatique avec bash

modifier tous les paramètres nécessaires dans le fichier workflow_1.sh

commande : 

*./workflow_1.sh*
