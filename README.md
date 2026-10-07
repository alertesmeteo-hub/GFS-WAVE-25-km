# GFS-Wave (NOAA) — vagues 0,25° pour la carte marine

Grilles de valeurs des vagues du modèle **GFS-Wave** (WaveWatch III, NOAA/NCEP), lues par la carte « Météo marine »
de https://www.alertes-meteo.com/meteo-cotier.

- Source : bucket public AWS `noaa-gfs-bdp-pds` (`gfs.AAAAMMJJ/HH/wave/gridded/`), runs 00/06/12/18 UTC, +120 h par pas de 3 h.
  Seuls les messages utiles sont téléchargés (plages d'octets lues dans l'index `.idx`).
- Sortie : branche `data` — `index.json` et `maps/index.json` (grilles de valeurs HKV1 : hauteur, période de pic,
  direction, mer du vent, houle principale).
- Workflow : `.github/workflows/update-gfswave.yml` (horaire, relancé aussi par le déclencheur du site).

```bash
pip install -r requirements.txt
cd scripts && python update_gfswave.py --output-dir ../build/national --force
```
