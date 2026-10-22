# SpaceX orbital-data-center constellation research generator

This research prototype uses parameters reviewed directly in the IRTF SPACE RG registry:
- `spacex-orbital-data-center.yaml` (Git blob SHA `b7561294fb870f5b8729e80b9a1e5e52a96851e8`): **filed** application, up to 1,000,000 satellites at altitudes 500–2000 km, 30° and sun-synchronous inclinations, shells up to 50 km thick. No plane count, phasing, or individual optical link capacity is specified. All generated ODC satellite layout and optical terminal specifications are **research assumptions**.
- `starlink-gen1.yaml` (Git blob SHA `8fbf37e48f309d51f67f91d3bc6db234fa486adb`): includes shell `gen1-s1` with 1,584 satellites, 72 planes, 550 km, 53°, Walker D, phase 39 **assumed** from literature.

Upstream records:
- https://github.com/irtf-spacerg/id-leo-constellations/blob/main/data/constellations/spacex-orbital-data-center.yaml
- https://github.com/irtf-spacerg/id-leo-constellations/blob/main/data/constellations/starlink-gen1.yaml

## Run on Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python generate.py --preset spacex-odc --satellites 60 --planes 6 --altitude 550 --inclination 30 --phasing 1 --duration 1800 --dt 30 --output output/odc
python generate.py --preset starlink-gen1 --duration 600 --dt 30 --output output/gen1
```

For a 600-satellite **synthetic ODC slice**:

```bash
python generate.py --preset spacex-odc --satellites 600 --planes 20 --altitude 550 --inclination 30 --phasing 1 --duration 3600 --dt 30 --terminals 4 --acquisition-s 30 --min-hold-s 120 --output output/odc600
```

## Output formats

- `satellites.csv`: configured orbital elements (not TLEs).
- `positions.csv`: ECI Cartesian positions in km at each sampled time.
- `candidate_contacts.csv`: geometrically feasible optical contacts within an assumed maximum range.
- `scheduled_contacts.csv`: subset chosen by a greedy nearest-neighbor scheduler with up to N optical terminals.
- `directed_trace.csv`: bidirectional contact/outage intervals for links activated at least once, suitable for **later** conversion to a Linux emulator schedule.
- `manifest.json`: provenance, numerical choices and limitations.

## Model caveats

This is a **first-pass research prototype**, not a physical reproduction of the SpaceX filing or a rigorous optical scheduler. It uses circular two-body motion in ECI (ignores J2 precession, atmospheric drag, eccentricity), spherical Earth line-of-sight and a hypothetical range limit/rate. In particular, this model does **not** correctly maintain sun-synchronous orbital planes. The optical terminal scheduler is greedy; acquisition, holding and link delay are sampled simplifications. Numerical contact start/end times have at most the chosen sample resolution, and delay is frozen at interval entry. The outputs should be validated and improved before making publication-quality contact-plan claims.

A separate subsequent module should implement event-driven continuous-time contacts, realistic laser terminal visibility/slew/acquisition, proper J2 or SGP4 orbital propagation, and explicit TVR YANG encoding. Do not present the CSV as a validated TVR-compliant schedule or as measured SpaceX operating data.
