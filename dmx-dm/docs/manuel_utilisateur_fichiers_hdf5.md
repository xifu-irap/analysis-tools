# Manuel utilisateur des fichiers HDF5 X-IFU

Lecture Python, extraction de configuration et exemples pratiques.

**Portée** : fichiers produits par DispatcherStudio et BackupConverter.
**API utilisée dans les exemples** : `h5py` + NumPy.

---

## 1. Vue d'ensemble

Les fichiers HDF5 du projet se répartissent en quatre familles. Les exemples de
ce manuel lisent les données sans charger inutilement un fichier complet en
mémoire.

| Type | Jeux de données principaux | Configuration |
| --- | --- | --- |
| Science continue | `ctrl [N]`, `pixels [N, 34]` | `/Configuration/DMXA` et `DMXB` |
| Scan | `ctrl [N]`, `x [N]`, `pixels [N, 34]` | `/Configuration/DMXA` et `DMXB` (v2) |
| Pulses | `/Pulses/FrameNum`, `col`, `pixel`, `pulse` | `/Configuration/DMX` |
| Dump | `Col0..Col3 [N, 1360]`, `Errors [N, 1360]` | Absente |

Prérequis : `pip install h5py numpy matplotlib`

Les tableaux sont indexés à partir de zéro en Python : les colonnes vont de 0 à
3 et les pixels de 0 à 33.

### Inspection rapide d'un fichier

```python
import h5py

def show_hdf5_tree(filename):
    with h5py.File(filename, "r") as h5:
        print("Attributs racine:", dict(h5.attrs))
        def visitor(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"/{name}: shape={obj.shape}, dtype={obj.dtype}")
            else:
                print(f"/{name}/")
        h5.visititems(visitor)

show_hdf5_tree("mon_fichier.h5")
```

---

## 2. Extraction de la configuration

Les valeurs stockées dans le groupe `/Configuration` sont déjà décodées. Les
paramètres fractionnaires sont exposés en nombres flottants ; les autres
conservent un type entier.

```python
def read_instrument_configuration(hdf5_filename):
    """Retourne /Configuration sous forme de dictionnaire Python."""
    import h5py
    import numpy as np

    def to_python(value):
        if isinstance(value, bytes):
            return value.decode("utf-8")
        if isinstance(value, np.ndarray):
            return [to_python(item) for item in value.tolist()]
        if isinstance(value, np.generic):
            return value.item()
        return value

    def read_group(group):
        result = {}
        if group.attrs:
            result["_attributes"] = {
                name: to_python(value)
                for name, value in group.attrs.items()
            }
        for name, item in group.items():
            if isinstance(item, h5py.Group):
                result[name] = read_group(item)
            else:
                result[name] = to_python(item[()])
        return result

    with h5py.File(hdf5_filename, "r") as h5:
        if "Configuration" not in h5:
            raise KeyError(
                f"{hdf5_filename!r} ne contient pas /Configuration"
            )
        return read_group(h5["Configuration"])
```

Un scan ancien au format v1 et un dump ne possèdent pas de groupe
`/Configuration`. La fonction l'indique explicitement avec une exception
`KeyError`.

> **Note d'implémentation** : dans `readData.py`, la fonction
> `read_instrument_configuration(hdf5_filename, required=False)` retourne
> `None` par défaut quand `/Configuration` est absent, pour simplifier le
> traitement des anciens fichiers. Passez `required=True` pour obtenir le
> comportement `KeyError` décrit ci-dessus.

---

## 3. Utiliser les paramètres de configuration

```python
configuration = read_instrument_configuration("scan_test_C0.h5")

# Paramètre scalaire de DMXA
mode = configuration["DMXA"]["DATA_ACQ_MODE"]
print("DATA_ACQ_MODE =", mode)

# Registre compact : 4 bits par colonne
offset_mode_register = configuration["DMXA"]["AMP_SQ_OFFSET_MODE"]
print("AMP_SQ_OFFSET_MODE =", offset_mode_register)

# Décoder la valeur de la colonne 2
column = 2
offset_mode_c2 = (offset_mode_register >> (4 * column)) & 0xF

# 34 pixels de la colonne 0
mux_fb0_c0 = configuration["DMXA"]["MUX_SQ_FB0"][0]

# Même paramètre pour les deux DMX
for dmx in ("DMXA", "DMXB"):
    print(dmx, configuration[dmx]["DATA_ACQ_MODE"])

# Sous-ensemble de paramètres
names = ["DATA_ACQ_MODE", "AMP_SQ_OFFSET_MODE", "MUX_SQ_FB_ON_OFF"]
selected = {name: configuration["DMXA"][name] for name in names}
print(selected)
```

### Métadonnées de configuration

```python
meta = configuration["_attributes"]
print("Version du format:", meta["FORMAT_VERSION"])
print("Date de capture (ms Unix):", meta["CAPTURE_DATE_MS"])
```

### Cas des fichiers de pulses

Un fichier de pulses ne contient que le DMX qui a produit l'acquisition. La clé
est donc `DMX`, et non `DMXA`/`DMXB`.

```python
configuration = read_instrument_configuration("pulses.h5")
mode = configuration["DMX"]["DATA_ACQ_MODE"]
gain_c1 = configuration["DMX"]["MUX_SQ_INPUT_GAIN"][1]
print(mode, gain_c1)
```

---

## 4. Fichiers de science continue

Chaque fichier correspond à une colonne DMX sélectionnée. Le suffixe `_C0.h5`
à `_C3.h5` identifie cette colonne.

| Chemin | Type / forme | Contenu |
| --- | --- | --- |
| `/ctrl` | `uint8 [N]` | Mot de contrôle de chaque trame |
| `/pixels` | `int16 [N, 34]` | 34 échantillons par trame |
| `/Configuration/DMXA` | groupe | Configuration DMX A complète |
| `/Configuration/DMXB` | groupe | Configuration DMX B complète |

### Lire une plage de trames

```python
import h5py
import numpy as np

with h5py.File("science_C0.h5", "r") as h5:
    # Lecture partielle : trames 1000 à 1999
    ctrl = h5["ctrl"][1000:2000]
    pixels = h5["pixels"][1000:2000, :]

    # Série temporelle du pixel 12
    pixel_12 = pixels[:, 12]

    # Moyenne de chaque pixel sur la plage
    mean_per_pixel = pixels.mean(axis=0)

print("Forme:", pixels.shape)
print("Moyenne pixel 12:", float(np.mean(pixel_12)))
```

### Sélectionner les trames par mot de contrôle

```python
with h5py.File("science_C0.h5", "r") as h5:
    ctrl = h5["ctrl"][:]
    wanted = (ctrl == 0xE0)
    selected_pixels = h5["pixels"][:][wanted]

print("Nombre de trames sélectionnées:", len(selected_pixels))
```

---

## 5. Fichiers de scan

Chaque fichier correspond à une colonne. Les lignes des trois datasets sont
alignées : l'indice `i` désigne la même trame dans `ctrl`, `x` et `pixels`.

| Chemin | Type / forme | Contenu |
| --- | --- | --- |
| `/ctrl` | `uint8 [N]` | Mot de contrôle |
| `/x` | `int32 [N]` | Consigne Offset ou Feedback |
| `/pixels` | `float32 [N, 34]` | Mesure des 34 pixels |
| `/Configuration/...` | groupes | Configuration DMXA et DMXB (scan v2) |

Attributs utiles : `COLUMN`, `COLUMN_MASK`, `SCNSRC`, `X_LABEL`,
`EXPECTED_FRAME_COUNT`, `ACTUAL_FRAME_COUNT`, `AMP_SQ_OFFSET_MODE` et
`TP_REG0` à `TP_REG4`.

### Lire et tracer un pixel

```python
import h5py
import matplotlib.pyplot as plt

pixel = 7
with h5py.File("scan_test_C0.h5", "r") as h5:
    x = h5["x"][:]
    y = h5["pixels"][:, pixel]
    x_label = h5.attrs["X_LABEL"]
    if isinstance(x_label, bytes):
        x_label = x_label.decode("utf-8")

plt.plot(x, y, marker=".")
plt.xlabel(x_label)
plt.ylabel(f"Pixel {pixel} (ADU)")
plt.grid(True)
plt.show()
```

---

## 6. Exploitation avancée d'un scan

### Vérifier que le scan est complet

```python
import h5py

with h5py.File("scan_test_C0.h5", "r") as h5:
    expected = int(h5.attrs["EXPECTED_FRAME_COUNT"])
    actual = int(h5.attrs["ACTUAL_FRAME_COUNT"])
    dataset_rows = h5["pixels"].shape[0]

complete = actual == expected == dataset_rows
print("Scan complet:", complete)
print(f"{actual} / {expected} trames")
```

### Moyenner les répétitions d'une même consigne

```python
import h5py
import numpy as np

with h5py.File("scan_test_C0.h5", "r") as h5:
    x = h5["x"][:]
    y = h5["pixels"][:, 7]

levels = np.unique(x)
mean_y = np.array([y[x == level].mean() for level in levels])
std_y = np.array([y[x == level].std() for level in levels])
for level, mean, std in zip(levels, mean_y, std_y):
    print(f"x={level}: {mean:.3f} +/- {std:.3f}")
```

### Lire les régions du test pattern

```python
def read_test_pattern(h5):
    regions = []
    for index in range(5):
        value = h5.attrs[f"TP_REG{index}"]
        if isinstance(value, bytes):
            value = value.decode("utf-8")
        regions.append([int(v) for v in value.split(",")])
    return regions

with h5py.File("scan_test_C0.h5", "r") as h5:
    test_pattern = read_test_pattern(h5)
    offset_mode = int(h5.attrs["AMP_SQ_OFFSET_MODE"])

print(test_pattern)
print("AMP_SQ_OFFSET_MODE =", offset_mode)
```

---

## 7. Fichiers de pulses

| Chemin | Type / forme | Contenu |
| --- | --- | --- |
| `/Pulses/FrameNum` | `uint64 [N]` | Numéro de trame |
| `/Pulses/col` | `uint8 [N]` | Colonne 0 à 3 |
| `/Pulses/pixel` | `uint8 [N]` | Pixel 0 à 33 |
| `/Pulses/pulse` | `int16 [N, S]` | S échantillons par pulse |
| `/Configuration/DMX` | groupe | Configuration du DMX source |

Les attributs racine `DATE`, `DATE_MS` et `DMX_ID` identifient l'acquisition. Le
groupe `/Pulses` possède `PULSE_SIZE` et `NROWS`.

### Lire les pulses d'un pixel

```python
import h5py
import numpy as np

with h5py.File("pulses.h5", "r") as h5:
    columns = h5["Pulses/col"][:]
    pixels = h5["Pulses/pixel"][:]
    selected = (columns == 1) & (pixels == 12)
    # h5py ne permet pas tous les masques multidimensionnels directement :
    indices = np.flatnonzero(selected)
    pulses = h5["Pulses/pulse"][indices, :]
    frame_numbers = h5["Pulses/FrameNum"][indices]

mean_pulse = pulses.mean(axis=0)
peak_per_pulse = pulses.max(axis=1)
print("Pulses trouvés:", len(indices))
print("Trames:", frame_numbers[:10])
print("Pic moyen:", float(peak_per_pulse.mean()))
```

---

## 8. Fichiers de dump

Un dump contient les quatre colonnes dans un seul fichier. Chaque ligne contient
1360 valeurs. Aucun groupe `/Configuration` n'est actuellement enregistré.

| Chemin | Type / forme | Contenu |
| --- | --- | --- |
| `/Col0 ... /Col3` | `int16 [N, 1360]` | Données des quatre colonnes |
| `/Errors` | `uint8 [N, 1360]` | Drapeaux d'erreur associés |
| Attribut `DATE` | texte | Date lisible |

---

## 9. Fonctions de lecture dans `readData.py`

Le module `readData.py` fournit des fonctions prêtes à l'emploi pour lire chaque
type de fichier. Les anciens fichiers (sans `/Configuration`) restent pris en
charge.

### Détection du type de fichier

```python
detect_hdf5_type(filename)
# -> "science" | "scan" | "pulses" | "dump" | "inconnu"
```

### Science continue

| Fonction | Retour |
| --- | --- |
| `get_science_from_hdf5(full_file_name)` | `(data, ctrl)` : `data` de forme `[34, N]` en unités s(16,2) (divisé par 4), `ctrl` mots de contrôle |
| `read_science_from_file(full_file_name, flatten=False, remove_dc=True, verbose=True)` | `col_data` : données de science d'une colonne |
| `read_col_science_from_dir(data_path, col_id, flatten=False, remove_dc=True, verbose=True)` | `(col_data, file_exists)` : lecture depuis un répertoire de fichiers `_C{col_id}.h5` |

### Scan

| Fonction | Retour |
| --- | --- |
| `read_scan(hdf5_file)` | `(x_name, ctrl, x_values, pixels_data)` : nom de l'axe X, mots de contrôle, consignes X, pixels `[34, N]` |
| `read_scan_type(hdf5_file)` | `x_name` : le nom de l'axe X |

### Dump

| Fonction | Retour |
| --- | --- |
| `read_dump_from_hdf5(hdf5_file)` | `(dump, adc_error)` : les quatre colonnes `[4, 1360]` et les erreurs ADC |

### Pulses

| Fonction | Retour |
| --- | --- |
| `read_pulses_from_hdf5(hdf5_file, column=None, pixel=None)` | `(frame_num, columns, pixels, pulses)` : numéros de trame, colonnes, pixels et échantillons `[M, S]` des pulses sélectionnés. Si `column`/`pixel` sont fournis, seuls les pulses correspondants sont retournés |

### Configuration instrument

| Fonction | Retour |
| --- | --- |
| `read_instrument_configuration(hdf5_filename, required=False)` | Le contenu de `/Configuration` sous forme de dict Python, ou `None` si absent |

---

## 10. Bonnes pratiques

1. Utiliser la lecture partielle `h5[dataset][start:stop]` pour les gros
   fichiers de science.
2. Ne pas supposer qu'un scan est complet : comparer `EXPECTED_FRAME_COUNT` et
   `ACTUAL_FRAME_COUNT`.
3. Tester l'existence de `/Configuration` avant de l'extraire. Les scans v1 et
   les dumps n'en ont pas.
4. Pour les scans, garder les datasets `x`, `ctrl` et `pixels` alignés avec le
   même découpage de lignes.
5. Pour les pulses, sélectionner les lignes à partir de `col` et `pixel` avant
   de lire `pulse`.
6. Les paramètres fractionnaires de `/Configuration` sont déjà convertis en
   valeurs physiques flottantes selon leur encodage fixe.

### Fonction utilitaire : type de fichier

```python
def detect_hdf5_type(filename):
    import h5py
    with h5py.File(filename, "r") as h5:
        if "Pulses" in h5:
            return "pulses"
        if "x" in h5 and "pixels" in h5:
            return "scan"
        if "ctrl" in h5 and "pixels" in h5:
            return "science"
        if all(f"Col{i}" in h5 for i in range(4)) and "Errors" in h5:
            return "dump"
        return "inconnu"
```

### Résumé des conventions

| Convention | Valeur |
| --- | --- |
| Colonnes | 0 à 3 |
| Pixels | 0 à 33 |
| Suffixes science/scan | `_C0.h5` à `_C3.h5` |
| Configuration science/scan | `/Configuration/DMXA` et `/DMXB` |
| Configuration pulse | `/Configuration/DMX` |

Document établi à partir des structures effectivement écrites par
DispatcherStudio et BackupConverter.
