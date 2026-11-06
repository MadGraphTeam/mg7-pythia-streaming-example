# MadGraph7 gridpack → Pythia → HepMC3, without parton-level files

`shower_to_hepmc3.py` streams unweighted `p p > j j j` events from a MadGraph7
gridpack into Pythia in memory and writes only the showered events to a HepMC3 file.

```
gridpack EventStream  →  StreamLHAup (LHAup subclass)  →  Pythia  →  HepMC3 file
        (in memory, no LHE file)
```

Pythia's Python interface only batches its *output* (`nextBatch`); it has no bulk
input, so Pythia pulls the parton-level events one at a time through `LHAup.setEvent`.

## Setup

Requires Python ≥ 3.12 and a C++ compiler.

**1. Python environment**

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt   # pythia8mc (>= 8.318), pyhepmc
```

**2. MadGraph7** (main branch; tested at `580dd8a07a`)

```sh
git clone git@github.com:MadGraphTeam/MadGraph7.git
git -C MadGraph7 checkout main
```

(To keep a separate checkout next to an existing clone: `git worktree add ../MadGraph7-main origin/main`.)

**3. Gridpack**

Event streaming needs `run_mode = fix_max_weight`. `make_gridpack.mg5` contains

```
generate p p > j j j
output mg7 PROC_ppjjj
launch PROC_ppjjj
set save_gridpack true
set run_mode fix_max_weight
done
```

Run it (the first run builds `madspace`, and the whole thing took about 7 min on 8 cores):

```sh
MadGraph7/bin/madgraph make_gridpack.mg5
```

The gridpack is written to `PROC_ppjjj/Events/run_01/gridpack` (≈ 50 MB).

## Run

```sh
.venv/bin/python shower_to_hepmc3.py \
    --gridpack PROC_ppjjj/Events/run_01/gridpack \
    -n 2000 -o ppjjj.hepmc3.gz --seed 7
```

- `-n`: number of showered events; `--seed`: seeds the gridpack stream and Pythia.
- `-o`: output file; a `.gz`/`.bz2`/`.xz` suffix compresses it (uncompressed ≈ 250 MB per 1000 events).
- The output holds the beams, hard process and final-state particles. The full
  shower history is not written: the HepMC3 ASCII writer of `pyhepmc` produces files
  its own reader rejects for full hadron-level records. Adjust `keep()` to change this.
- Read it back with `pyhepmc.open("ppjjj.hepmc3.gz")`. The event weights are 1 (larger for overweight events); normalise to the printed cross section.

## Notes

- This is an inclusive 3-jet sample showered without matching or merging, a demo rather than a physics-ready setup.
- The script unsets `PYTHIA8DATA` (e.g. set by a conda `pythia8`), which otherwise conflicts with the XML data bundled in the `pythia8mc` wheel.
- The conda `pythia8` package (8.312) has no `LHAup` binding, so the pip `pythia8mc` is required.
