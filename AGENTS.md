# AGENTS.md

Guide de conventions pour l'assistant et ses agents travaillant sur le projet
**analysis-tools** (ATHENA X-IFU DRE data analysis tools, IRAP Toulouse).

## Vue d'ensemble du projet

- Tous les scripts Python se trouvent dans `dmx-dm/` (analyse de données DEMUX :
  bruit, non-linéarité, délais, scans, HK/température, etc.).
- `dmx-dm/constants.py` : constantes partagées (dimensions DEMUX, fréquences,
  répertoires, modèles).
- `dmx-dm/general_tools.py` : utilitaires généraux.
- `dmx-dm/readData.py` : lecture des données (science, dump, scan, HK).
- `dmx-dm/noise_models/` : modèles de bruit.
- Des notebooks Jupyter (`.ipynb`) sont aussi présents dans `dmx-dm/`.

## Conventions de codage Python

Appliquer les recommandations de codage Python standard (PEP 8 et bonnes
pratiques) sur **tout code nouveau ou modifié** :

- **PEP 8** : lignes ≤ 79 caractères, indentations à 4 espaces, deux lignes
  vides entre les fonctions de module.
- **Type hints** : annoter les paramètres et valeurs de retour de toutes les
  fonctions.
- **Docstrings** : format NumPy, rédigées en anglais, avec les sections
  `Parameters`, `Returns` (et `Raises` si pertinent).
- **Nommage** : `snake_case` pour les fonctions et variables. Ne pas renommer
  l'API publique existante (les signatures doivent rester compatibles avec les
  appelants).
- **Imports** : en tête de fichier, groupés (stdlib, tiers, locaux), sans
  import dans le corps des fonctions sauf nécessité.
- **Commentaires** : en anglais de préférence.
- Ne pas ajouter de commentaires redondants au code.

Note : une partie du code existant est encore en style ancien (pas de type
hints, continuations par backslash). Les modifications doivent moderniser ce
qu'elles touchent, sans casser l'API.

## Licence et en-tête de fichier

Tous les fichiers Python commencent par l'en-tête de licence GPL-3.0 existant
(bloc `# ------...`). Conserver cet en-tête intact lors des modifications et le
reproduire pour tout nouveau fichier `.py`.

## Vérifications

- Syntaxe : `python -m py_compile <fichier.py>`
- Style : `python -m flake8 <fichier.py>` (si `flake8` est disponible)
- Import : `python -c "import <module>"` dans `dmx-dm/`

Exécuter ces vérifications après toute modification d'un fichier Python.

## Gestion des dépendances

Aucun gestionnaire de paquets ni environnement virtuel n'est configuré dans le
dépôt. Ne pas en introduire sans demande explicite.

## Conventions de commit

- Ne committer/pusher que sur demande explicite de l'utilisateur.
- Message de commit concis, en anglais, décrivant le changement.
