#!/usr/bin/env python3
"""test_spice_simple.py - Quick SPICE verification"""

import os
import spiceypy as spice

# Load kernels
kernel_dir = os.path.expanduser("~/kernels")
kernels = ["naif0012.tls", "de440.bsp", "gm_de440.tpc", "pck00010.tpc"]

for k in kernels:
    path = os.path.join(kernel_dir, k)
    if os.path.exists(path):
        spice.furnsh(path)
        print(f"✅ {k}")
    else:
        print(f"❌ {k} not found")

# Test getting GM values
try:
    gm_earth = spice.bodvcd(399, "GM", 1)
    print(f"✅ GM Earth: {gm_earth[1][0]:.4e} m³/s²")

    gm_moon = spice.bodvcd(301, "GM", 1)
    print(f"✅ GM Moon: {gm_moon[1][0]:.4e} m³/s²")

    # Test Moon state
    et = spice.utc2et("2026-08-27T12:00:00")
    state, lt = spice.spkezr("MOON", et, "J2000", "NONE", "EARTH")
    print(f"✅ Moon position: ({state[0]:.1f}, {state[1]:.1f}) km")

    print("\n✅ SPICE working correctly")

except Exception as e:
    print(f"❌ Error: {e}")