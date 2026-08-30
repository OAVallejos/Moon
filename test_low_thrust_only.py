#!/usr/bin/env python3
"""
test_low_thrust_only.py - L-RAIL TUG-BIG: Orbital descent 800→200 km
Uses the compiled .so kernel (tangential, 2000 segments)
TOF=6.97 days for correct descent
"""

import ctypes
import os
import sys
import math

# ============================================================
# PRELOAD CSPICE
# ============================================================
def preload_cspice():
    conda_prefix = os.environ.get("CONDA_PREFIX", os.path.expanduser("~/micromamba/envs/Astro"))
    lib_paths = [
        os.path.join(conda_prefix, "lib", "libcspice.so"),
        "/usr/local/lib/libcspice.so",
    ]
    for lib_path in lib_paths:
        if os.path.exists(lib_path):
            ctypes.CDLL(lib_path, mode=ctypes.RTLD_GLOBAL)
            print(f"✅ CSPICE loaded from: {lib_path}")
            return True
    return False

preload_cspice()

# ============================================================
# IMPORT ASTRO_TFC (already compiled kernel)
# ============================================================
sys.path = [p for p in sys.path if 'site-packages' not in p]
sys.path.insert(0, os.path.expanduser("~/rust"))
sys.path.insert(0, os.path.expanduser("~/rust/target/release"))

try:
    import astro_tfc
    print(f"✅ ASTRO-TFC imported from: {astro_tfc.__file__}")
except ImportError:
    import importlib.util
    so_path = os.path.expanduser("~/rust/target/release/libastro_tfc.so")
    if os.path.exists(so_path):
        spec = importlib.util.spec_from_file_location("astro_tfc", so_path)
        astro_tfc = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(astro_tfc)
        print(f"✅ ASTRO-TFC loaded from: {so_path}")
    else:
        raise

# ============================================================
# CRTBP CONSTANTS
# ============================================================
MU = 0.01215058560962404
X_MOON = 1.0 - MU
X_EARTH = -MU

D_CHAR = astro_tfc.D_CHAR  # meters
V_CHAR = astro_tfc.V_CHAR  # m/s
T_CHAR = astro_tfc.T_CHAR  # seconds
D_CHAR_km = D_CHAR / 1000.0

r_moon_km = 1737.0
r_moon_can = r_moon_km / D_CHAR_km

alt_800_km = 800.0
alt_800_can = alt_800_km / D_CHAR_km
r_800_can = r_moon_can + alt_800_can

print("\n" + "=" * 80)
print("🚀 L-RAIL TUG-BIG: ORBITAL DESCENT 800→200 km")
print("   Compiled .so kernel | 2000 segments | Retrograde tangential")
print("=" * 80)

# ============================================================
# MISSION PARAMETERS
# ============================================================
TOF_DESC = 6.97  # days (calculated for 800→200 km)
THRUST = 0.64    # N (8x NEXT-C)
ISP = 4190.0     # s
MASS = 1910.5    # kg
G0 = 9.80665     # m/s²

xenon_expected = THRUST / (ISP * G0) * TOF_DESC * 86400.0

print(f"\n📊 MANEUVER PARAMETERS:")
print(f"   TOF: {TOF_DESC} days")
print(f"   Thrust: {THRUST*1000:.0f} mN")
print(f"   Expected xenon: {xenon_expected:.1f} kg")

# ============================================================
# INITIAL STATE: CIRCULAR LLO ORBIT 800 km
# ============================================================
rho_x = 0.0
rho_y = r_800_can
rho_z = 0.0
r_rel = math.sqrt(rho_x**2 + rho_y**2 + rho_z**2)

v_circ = math.sqrt(MU / r_rel)
v_in_x = -v_circ
v_in_y = 0.0

omega = 1.0
vx_synodic = v_in_x + rho_y
vy_synodic = v_in_y - rho_x

state_llo = [
    X_MOON + rho_x,
    rho_y,
    rho_z,
    vx_synodic,
    vy_synodic,
    0.0
]

cj = astro_tfc.compute_jacobi_constant(state_llo)

dx = state_llo[0] - X_MOON
dy = state_llo[1]
alt_check = (math.sqrt(dx*dx + dy*dy) - r_moon_can) * D_CHAR_km

print(f"\n📐 INITIAL STATE:")
print(f"   C_J: {cj:.6f}")
print(f"   Altitude: {alt_check:.1f} km")
print(f"   Position: ({state_llo[0]:.6f}, {state_llo[1]:.6f})")
print(f"   Velocity: ({state_llo[3]:.6f}, {state_llo[4]:.6f})")

# ============================================================
# MISSION CONFIGURATION
# ============================================================
config = astro_tfc.MissionConfig()
config.verbose = True
config.spacecraft_mass_kg = MASS
config.engine_thrust_n = THRUST
config.engine_isp_s = ISP

mission = astro_tfc.AstroTFCMission(config)

# ============================================================
# MANEUVER: ORBITAL DESCENT
# ============================================================
print("\n" + "=" * 80)
print("📡 ORBITAL DESCENT (800 km → 200 km)")
print("=" * 80)

print(f"\n🔧 Executing descent...")
print(f"   TOF: {TOF_DESC} days")
print("-" * 80)

try:
    trajectory, controls = mission.optimize_low_thrust(
        state_llo,
        200.0,
        TOF_DESC
    )

    # Verify final altitude
    final = trajectory[-1]
    dx = final[0] - X_MOON
    dy = final[1]
    dz = final[2]
    r_final = math.sqrt(dx*dx + dy*dy + dz*dz)
    alt_final = (r_final - r_moon_can) * D_CHAR_km

    # Actual xenon consumption
    mass_current = MASS
    dt_seconds = TOF_DESC * 86400.0 / len(controls)
    for c in controls:
        m_dot = c[0] * THRUST / (ISP * G0)
        mass_current -= m_dot * dt_seconds
    xenon_actual = MASS - mass_current

    print(f"\n✅ DESCENT COMPLETED")
    print(f"   Final altitude: {alt_final:.1f} km (target: 200 km)")
    print(f"   Error: {abs(alt_final - 200):.1f} km")
    print(f"   Xenon consumed: {xenon_actual:.1f} kg")
    print(f"   Match: {'✅' if abs(alt_final - 200) < 10 else '⚠️'}")

    # Control analysis
    total_throttle = sum(c[0] for c in controls)
    avg_throttle = total_throttle / len(controls)
    max_throttle = max(c[0] for c in controls)
    active = sum(1 for c in controls if c[0] > 0.01)

    print(f"\n📈 CONTROL PROFILE:")
    print(f"   Average throttle: {avg_throttle:.3f}")
    print(f"   Maximum throttle: {max_throttle:.3f}")
    print(f"   Active segments: {active}/{len(controls)}")

except Exception as e:
    print(f"\n❌ Error: {e}")

print("\n✅ TEST COMPLETED")
print("=" * 80)