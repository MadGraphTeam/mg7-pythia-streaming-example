#!/usr/bin/env python
"""Showered p p > j j j events from a MadGraph7 gridpack, written to a HepMC3 file.

    MadGraph7 EventStream -> StreamLHAup -> Pythia -> HepMC3 file

The parton-level events are handed to Pythia in memory by ``StreamLHAup``, a Python
subclass of Pythia's ``LHAup`` (``Beams:frameType = 5``): no LHE file or pipe. Pythia
pulls them one at a time, since its Python interface has no bulk input. Only the
showered/hadronised events (beams, hard process and final state) are written, as
HepMC3 ASCII.

Requirements: pythia8mc (the conda ``pythia8`` has no ``LHAup`` binding), pyhepmc, and a
gridpack made with ``run_mode = fix_max_weight``; see README.md.

    python shower_to_hepmc3.py --gridpack PROC_ppjjj/Events/run_01/gridpack \
        -n 2000 -o ppjjj.hepmc3
"""

import argparse
import os
import sys

# the conda PYTHIA8DATA points at a different Pythia version than the pip wheel
os.environ.pop("PYTHIA8DATA", None)

import pyhepmc
import pythia8mc


class StreamLHAup(pythia8mc.LHAup):
    """Serves the parton-level events of a gridpack ``EventStream`` to Pythia as
    external Les Houches events. Pythia pulls events one at a time through
    ``setEvent`` (the Python interface has no bulk input), so the events are taken
    from the stream in chunks and passed on as they are, without conversion."""

    def __init__(self, stream, chunk=1000):
        super().__init__()
        self.stream, self.chunk = stream, chunk
        self.events, self.pos = [], 0

    def setInit(self):
        meta = self.stream.lhe_meta()
        self.setBeamA(meta.beam1_pdg_id, meta.beam1_energy, 0, meta.beam1_pdf_id)
        self.setBeamB(meta.beam2_pdg_id, meta.beam2_energy, 0, meta.beam2_pdf_id)
        self.setStrategy(meta.weight_mode)  # 3: unweighted events, weight +-1
        for proc in meta.processes:
            self.addProcess(proc.process_id, proc.cross_section,
                            proc.cross_section_error, proc.max_weight)
        return True

    def setEvent(self, idProcIn=0):
        if self.pos >= len(self.events):
            self.events = self.stream.next_events(self.chunk)  # list of ms.LHEEvent
            self.pos = 0
        ev = self.events[self.pos]
        self.pos += 1
        self.setProcess(ev.process_id, ev.weight, ev.scale, ev.alpha_qed, ev.alpha_qcd)
        for p in ev.particles:
            self.addParticle(p.pdg_id, p.status_code, p.mother1, p.mother2,
                             p.color, p.anti_color, p.px, p.py, p.pz,
                             p.energy, p.mass, p.lifetime, p.spin, -1.0)
        return True


def keep(event, i):
    """Particles written to HepMC3: beams (11-13), the hard process (21-29) and the
    final state. The full shower/hadronisation history makes the HepMC3 ASCII
    writer of pyhepmc emit files that its own reader rejects (vertices missing)."""
    particle = event[i]
    return particle.isFinal() or 11 <= particle.statusAbs() <= 13 \
        or 21 <= particle.statusAbs() <= 29


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gridpack", required=True)
    ap.add_argument("-n", "--events", type=int, default=2000)
    ap.add_argument("-o", "--output", default="ppjjj.hepmc3",
                    help="output file; a .gz/.bz2/.xz suffix compresses it")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    sys.path.insert(0, os.path.join(os.path.abspath(args.gridpack), "bin"))
    from event_stream import EventStream

    with EventStream(seed=args.seed) as stream:
        lha = StreamLHAup(stream)
        pythia = pythia8mc.Pythia("", False)
        for line in [
            "Beams:frameType = 5",
            "Random:setSeed = on",
            f"Random:seed = {args.seed + 1}",
            "Next:numberCount = 0",
            "PartonLevel:MPI = on",
            "HadronLevel:Hadronize = on",
        ]:
            pythia.readString(line)
        pythia.setLHAupPtr(lha)
        pythia.init()

        convert = pythia8mc.Pythia8ToFilteredHepMC3()
        convert.selector = keep
        with pyhepmc.open(args.output, "w") as out:
            n_done = 0
            while n_done < args.events:
                if not pythia.next():
                    continue
                event = convert.hepmc3.GenEvent()
                convert.fill_next_event(event, ievnum=n_done, pythia=pythia)
                out.write(event)
                n_done += 1
                if n_done % 500 == 0:
                    print(f"  {n_done}/{args.events} events written")

        pythia.stat()
        print(f"\nwrote {n_done} events to {args.output}; parton-level cross section "
              f"{stream.integral:.4g} +- {stream.integral_error:.2g} pb")


if __name__ == "__main__":
    main()
