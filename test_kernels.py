# test_kernels.py
import os
import spiceypy as spice
from datetime import datetime, timezone

# 1. Load the meta-kernel set
kernel_dir = os.path.expanduser("~/kernels")
for k in ["naif0012.tls", "de440.bsp", "pck00010.tpc"]:
    spice.furnsh(os.path.join(kernel_dir, k))

# 2. Current UTC time -> ET
utc_now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
et = spice.utc2et(utc_now)

# 3. State vector of Mars (499) as seen from Earth (399) in J2000
# spkezr returns [x, y, z, vx, vy, vz] and light time (lt)

state, lt = spice.spkezr('MARS BARYCENTER', et, 'J2000', 'LT+S', 'EARTH')

posicion = state[:3]      # km
velocidad = state[3:]     # km/s
distancia_km = spice.vnorm(posicion)
distancia_au = distancia_km / 149597870.7

# Range rate: scalar projection of velocity onto the position vector direction
range_rate = spice.vdot(posicion, velocidad) / distancia_km

print(f"--- MISSION STATUS: OK ---")
print(f"Date (UTC):         {utc_now}")
print(f"Distance to Mars:   {distancia_km:,.2f} km ({distancia_au:.4f} AU)")
print(f"Light time (LT):    {lt/60:.2f} min")
print(f"Relative velocity:  {range_rate:+.3f} km/s "
      f"({'Moving away' if range_rate > 0 else 'Approaching'})")

# Memory cleanup
spice.kclear()