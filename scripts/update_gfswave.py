#!/usr/bin/env python3
"""Vagues GFS-Wave (NOAA / NCEP, WaveWatch III global 0,25°) : grilles de valeurs pour la carte marine du site.

Source : bucket public AWS `noaa-gfs-bdp-pds` (gfs.AAAAMMJJ/HH/wave/gridded/gfswave.tHHz.global.0p25.fFFF.grib2),
runs 00, 06, 12, 18 UTC. Chaque fichier est accompagné d'un index .idx : on ne télécharge que les messages utiles
(plages d'octets). Les champs sont rééchantillonnés sur une grille Web Mercator et publiés au format « HKV1 »
(voir wavegrid.py), avec un manifeste maps/index.json de même forme que celui de MFWAM.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import requests
from eccodes import codes_get, codes_get_double_array, codes_new_from_message, codes_release

from wavegrid import MercatorResampler, RegularGrid, write_hkv

LOGGER = logging.getLogger("gfswave")
BUCKET = "https://noaa-gfs-bdp-pds.s3.amazonaws.com/"
USER_AGENT = "alertes-meteo.com/gfswave-noaa/1.0"
PIPELINE_VERSION = "1.0.0"
DEFAULT_CURRENT_METADATA_URL = "https://raw.githubusercontent.com/alertesmeteo-hub/GFS-WAVE-25-km/data/index.json"

BOUNDS = {"south": 30.0, "west": -10.0, "north": 50.0, "east": 20.0}
PROBE_WIDTH = 120  # maille native 0,25°

# Clé de la grille publiée -> (libellé idx NCEP, direction ?). Mêmes clés que MFWAM / EWAM.
# SWELL / SWPER / SWDIR « 1 in sequence » = première partition de houle (la plus énergétique).
PROBES = {
    "hs": ("HTSGW:surface", False),
    "tp_pic": ("PERPW:surface", False),
    "dir": ("DIRPW:surface", True),
    "wind_h": ("WVHGT:surface", False),
    "wind_tp": ("WVPER:surface", False),
    "wind_dir": ("WVDIR:surface", True),
    "swell_h": ("SWELL:1 in sequence", False),
    "swell_tp": ("SWPER:1 in sequence", False),
    "swell_dir": ("SWDIR:1 in sequence", True),
}
# Échéances : toutes les 3 h jusqu'à +120 h.
LEADS = list(range(0, 121, 3))
MISSING = 1e10


def iso(value: datetime) -> str:
    return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def file_url(run: datetime, lead: int) -> str:
    return f"{BUCKET}gfs.{run:%Y%m%d}/{run:%H}/wave/gridded/gfswave.t{run:%H}z.global.0p25.f{lead:03d}.grib2"


def latest_complete_run(session: requests.Session) -> datetime | None:
    """Run le plus récent dont la dernière échéance publiée (+120 h) est disponible."""
    now = datetime.now(timezone.utc)
    candidates = []
    for days_back in (0, 1):
        day = (now - timedelta(days=days_back)).replace(minute=0, second=0, microsecond=0)
        for hour in (18, 12, 6, 0):
            run = day.replace(hour=hour)
            if run <= now:
                candidates.append(run)
    for run in sorted(candidates, reverse=True):
        response = session.head(file_url(run, LEADS[-1]) + ".idx", timeout=30)
        if response.status_code == 200:
            return run
    return None


def already_published(url: str, run: datetime) -> bool:
    try:
        response = requests.get(url, timeout=30, headers={"Cache-Control": "no-cache"})
        if response.status_code != 200:
            return False
        return response.json().get("model", {}).get("run_time") == iso(run)
    except (requests.RequestException, ValueError):
        return False


def get_bytes(session: requests.Session, url: str, byte_range: str | None = None) -> bytes:
    headers = {"Range": f"bytes={byte_range}"} if byte_range else {}
    last: Exception | None = None
    for attempt in range(4):
        try:
            response = session.get(url, headers=headers, timeout=(15, 180))
            if response.status_code in (200, 206):
                return response.content
            last = RuntimeError(f"HTTP {response.status_code}")
        except requests.RequestException as exc:
            last = exc
        if attempt < 3:
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"Téléchargement impossible : {url} ({last})")


def parse_idx(text: str) -> list[tuple[str, int]]:
    """[(« NOM:niveau », octet de début)] dans l'ordre du fichier."""
    rows = []
    for line in text.strip().splitlines():
        parts = line.split(":")
        if len(parts) >= 5:
            rows.append((f"{parts[3]}:{parts[4]}", int(parts[1])))
    return rows


def read_field(data: bytes) -> tuple[np.ndarray, RegularGrid]:
    gid = codes_new_from_message(data)
    try:
        ni = int(codes_get(gid, "Ni"))
        nj = int(codes_get(gid, "Nj"))
        lat_first = float(codes_get(gid, "latitudeOfFirstGridPointInDegrees"))
        lat_last = float(codes_get(gid, "latitudeOfLastGridPointInDegrees"))
        lon_first = float(codes_get(gid, "longitudeOfFirstGridPointInDegrees"))
        di = float(codes_get(gid, "iDirectionIncrementInDegrees"))
        dj = float(codes_get(gid, "jDirectionIncrementInDegrees"))
        values = codes_get_double_array(gid, "values").reshape(nj, ni).astype(np.float64)
    finally:
        codes_release(gid)
    values[values >= MISSING] = np.nan
    grid = RegularGrid(lat_first=lat_first, lon_first=lon_first, lat_step=-dj if lat_first > lat_last else dj, lon_step=di, nj=nj, ni=ni)
    return values, grid


def fetch_step(session: requests.Session, run: datetime, lead: int) -> dict[str, tuple[np.ndarray, RegularGrid]]:
    url = file_url(run, lead)
    idx = parse_idx(get_bytes(session, url + ".idx").decode("utf-8"))
    starts = [start for _, start in idx]
    wanted = {}
    for key, (label, _) in PROBES.items():
        position = next((i for i, (name, _) in enumerate(idx) if name.startswith(label)), None)
        if position is None:
            raise RuntimeError(f"Champ {label} absent de l'index à +{lead} h")
        end = starts[position + 1] - 1 if position + 1 < len(starts) else ""
        wanted[key] = f"{starts[position]}-{end}"
    with ThreadPoolExecutor(max_workers=6) as pool:
        payloads = dict(zip(wanted, pool.map(lambda r: get_bytes(session, url, r), wanted.values())))
    return {key: read_field(data) for key, data in payloads.items()}


def build(run: datetime, workdir: Path) -> Path:
    result = workdir / "result"
    maps = result / "maps"
    maps.mkdir(parents=True)
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})

    resampler: MercatorResampler | None = None
    steps = []
    for position, lead in enumerate(LEADS, start=1):
        fields = fetch_step(session, run, lead)
        if resampler is None:
            resampler = MercatorResampler(BOUNDS, PROBE_WIDTH, fields["hs"][1], lon_wrap=True)
            LOGGER.info("Grille Mercator %sx%s", resampler.width, resampler.height)
        written: dict[str, str] = {}
        for key, (_, is_direction) in PROBES.items():
            sampled = resampler.sample(fields[key][0], nearest=is_direction)
            relative = f"maps/values/{key}/{lead:03d}.hkv.gz"
            if write_hkv(result / relative, sampled):
                written[key] = relative
        if "hs" not in written:
            raise RuntimeError(f"Hauteur significative vide à +{lead} h")
        steps.append({"lead_hour": lead, "valid_time": iso(run + timedelta(hours=lead)), "files": {}, "probes": written})
        LOGGER.info("Échéance +%03d h (%s/%s)", lead, position, len(LEADS))

    assert resampler is not None
    generated_at = iso(datetime.now(timezone.utc))
    manifest = {
        "schema_version": 1,
        "status": "ok",
        "module_version": PIPELINE_VERSION,
        "generated_at": generated_at,
        "run_time": iso(run),
        "projection": "EPSG:3857",
        "bounds": BOUNDS,
        "probe_grid": {"width": resampler.width, "height": resampler.height},
        "layers": {},
        "steps": steps,
    }
    with (maps / "index.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, separators=(",", ":"))
        handle.write("\n")
    index = {
        "schema_version": 1,
        "status": "ok",
        "generated_at": generated_at,
        "model": {
            "name": "GFS-Wave",
            "provider": "NOAA / NCEP",
            "dataset": "GFS Wave (WaveWatch III) global 0,25°",
            "resolution_km": 28,
            "run_time": iso(run),
            "pipeline_version": PIPELINE_VERSION,
            "source_url": BUCKET,
            "license": "Domaine public (NOAA)",
        },
        "coverage": {"label": "Méditerranée occidentale et centrale, golfe de Gascogne, Manche (30N-50N, 10W-20E)"},
        "maps": {"status": "ok", "manifest": "maps/index.json", "steps": len(steps)},
    }
    with (result / "index.json").open("w", encoding="utf-8") as handle:
        json.dump(index, handle, ensure_ascii=False, separators=(",", ":"))
        handle.write("\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="build/national")
    parser.add_argument("--current-metadata-url", default=DEFAULT_CURRENT_METADATA_URL)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=os.getenv("LOG_LEVEL", "INFO"), format="%(asctime)s | %(levelname)s | %(message)s")

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    run = latest_complete_run(session)
    if run is None:
        LOGGER.info("Aucun run GFS-Wave complet pour le moment.")
        return 0
    LOGGER.info("Run GFS-Wave complet le plus récent : %s", iso(run))
    if not args.force and already_published(args.current_metadata_url, run):
        LOGGER.info("Ce run est déjà publié, rien à faire.")
        return 0
    with tempfile.TemporaryDirectory(prefix="gfswave-build-") as tmp:
        result = build(run, Path(tmp))
        destination = Path(args.output_dir)
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(result, destination)
    LOGGER.info("Fichiers prêts dans %s", args.output_dir)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        LOGGER.exception("Échec de la mise à jour GFS-Wave")
        raise SystemExit(1)
