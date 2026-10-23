
import math
from datetime import datetime, timezone
from sgp4.api import Satrec, WGS72
from sgp4.api import jday

MU = 398600.8  # km^3/s^2, WGS-72
RE = 6378.135   # km, WGS-72

def build_sgp4_satellites(nodes, epoch):
    jd, fr = jday(
        epoch.year, epoch.month, epoch.day,
        epoch.hour, epoch.minute,
        epoch.second + epoch.microsecond / 1e6
    )

    epoch_sgp4 = jd + fr - 2433281.5
    satellites = []

    for i, node in enumerate(nodes):
        radius = RE + node["altitude_km"]

        # Initial approximation for circular mean motion
        mean_motion = math.sqrt(MU / radius**3) *60 # rad/min 

        sat = Satrec()
        sat.sgp4init(
            WGS72,
            "i",
            i + 1,
            epoch_sgp4,
            0.0,  # BSTAR drag parameter
            0.0,  # ndot
            0.0,  # nddot
            0.0,  # eccentricity
            math.radians(node["raan_deg"] * 0),
            math.radians(node["inclination_deg"]),
            math.radians(node["argument_latitude_deg"]),
            mean_motion,
            math.radians(node["raan_deg"]),
        )

        satellites.append(sat)

    return satellites



import numpy as np
from datetime import timedelta

def sgp4_positions(satellites, epoch, t):
    current_time = epoch + timedelta(seconds=float(t))

    jd, fr = jday(
        current_time.year,
        current_time.month,
        current_time.day,
        current_time.hour,
        current_time.minute,
        current_time.second
        + current_time.microsecond / 1e6
    )

    positions = []

    for sat in satellites:
        error, r, v = sat.sgp4(jd, fr)

        if error != 0:
            raise RuntimeError(
                f"SGP4 propagation error: {error}"
            )

        positions.append(r)

    return np.asarray(positions)
