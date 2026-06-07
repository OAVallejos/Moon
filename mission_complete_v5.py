# mission_complete_v5.py
#!/usr/bin/env python3     
"""                        
AstroTFC Mission: Adaptive sweep at the edge of chaos.                     
Version v5 FINAL: Using CORRECT canonical C_J (~3.188).

SOLUTION: The backend expects C_J ≈ 3.188, not ≈ 162.
compute_jacobi_constant() returns dimensional value, we need the canonical one.
"""
import astro_tfc
import sys
import numpy as np
import math

# ------------------------------------------------------------
# Physical constants
# ------------------------------------------------------------
MU_M = astro_tfc.MU_M
TARGET_ALT = 200.0  # km
MASS_0 = 250.0
THRUST = 0.08
ISP = 1500.0
G0 = 9.80665
VE = ISP * G0
FUEL_INJECTION = 3.0

# CR3BP canonical constants
D_CHAR = 384400000.0  # m (Earth-Moon distance)
V_CHAR = 1024.5       # m/s (lunar orbital velocity, more precise value)
MU_EARTH = 3.986004418e14  # m³/s²

def circular_velocity(alt_km):
    """Circular velocity at a given altitude above the Moon."""
    r = astro_tfc.R_MOON + alt_km * 1000.0
    return math.sqrt(MU_M / r)

def spiral_dv(alt_initial, alt_final):
    """Delta-V for spiral between two lunar altitudes."""
    alt_i = max(abs(alt_initial), TARGET_ALT)
    return abs(circular_velocity(alt_final) - circular_velocity(alt_i))

def total_mission(capture_alt_km, transit_days, dv_extra=0.0):
    """Calculates total fuel consumption and mission duration."""
    m_after_inj = MASS_0 - FUEL_INJECTION
    t_inj = (FUEL_INJECTION / (THRUST / VE)) / 86400.0

    if dv_extra > 0:
        m_after_brake = m_after_inj / math.exp(dv_extra / VE)
        fuel_brake = m_after_inj - m_after_brake
        t_brake = fuel_brake / (THRUST / VE) / 86400.0
    else:
        m_after_brake = m_after_inj
        fuel_brake = 0.0
        t_brake = 0.0

    dv_spiral = spiral_dv(capture_alt_km, TARGET_ALT)
    m_final = m_after_brake / math.exp(dv_spiral / VE)
    fuel_spiral = m_after_brake - m_final
    t_spiral = fuel_spiral / (THRUST / VE) / 86400.0

    total_xenon = FUEL_INJECTION + fuel_brake + fuel_spiral
    total_days = t_inj + transit_days + t_brake + t_spiral

    return total_xenon, total_days

def compute_jacobi_canonical(state_canonical):
    """
    Calculates the Jacobi Constant in CANONICAL units.

    Formula: C_J = 2U - v²
    where U = (x² + y²)/2 + (1-μ)/r1 + μ/r2
    and everything in canonical units.
    """
    x, y, z, vx, vy, vz = state_canonical

    # Effective potential in canonical units
    mu = 0.01215058560962404  # Lunar mass parameter

    r1 = math.sqrt((x + mu)**2 + y**2 + z**2)  # Distance to Earth
    r2 = math.sqrt((x - 1 + mu)**2 + y**2 + z**2)  # Distance to Moon

    U = 0.5 * (x**2 + y**2) + (1 - mu) / r1 + mu / r2

    v_squared = vx**2 + vy**2 + vz**2

    cj = 2.0 * U - v_squared

    return cj

def get_l1_jacobi_canonical():
    """
    Calculates the C_J at L1 using the canonical position of L1.
    The backend reports: "Theoretical minimum: 3.188341"
    """
    xl1 = astro_tfc.get_l1_position()

    # State at L1 (zero velocity)
    state_l1 = [xl1, 0.0, 0.0, 0.0, 0.0, 0.0]

    # Calculate canonical C_J
    cj_l1_canonical = compute_jacobi_canonical(state_l1)

    return cj_l1_canonical, xl1

def test_injection_direct(cj_target, verbose=True):
    """
    Direct injection test with C_J in canonical units.

    Args:
        cj_target: Target C_J in CANONICAL units (~3.18)
    """
    try:
        config = astro_tfc.MissionConfig()
        config.target_jacobi = cj_target  # Pass canonical value
        config.verbose = verbose

        mission = astro_tfc.AstroTFCMission(config)

        if verbose:
            print(f"   🎯 Target C_J = {cj_target:.6f} (canonical)")

        s1 = mission.execute_stage_1()
        if not s1.success:
            error_msg = getattr(s1, 'error', 'Unknown')
            print(f"   ❌ Stage 1 failed: {error_msg}")
            return None

        if verbose:
            print(f"   ✅ Injection: ΔV={s1.dv_total_ms:.1f} m/s")

        s2 = mission.execute_stage_2()
        if not s2.success:
            error_msg = getattr(s2, 'error', 'Unknown')
            print(f"   ❌ Stage 2 failed: {error_msg}")
            return None

        s3 = mission.execute_stage_3()

        # Extract information
        capture_alt = 500.0
        tof = s2.time_of_flight_days if hasattr(s2, 'time_of_flight_days') else 10.0

        if hasattr(s3, 'capture_altitude_km'):
            capture_alt = s3.capture_altitude_km
        elif hasattr(s3, 'details'):
            import re
            numbers = re.findall(r'\d+\.?\d*', str(s3.details))
            if numbers:
                capture_alt = float(numbers[0])

        return {
            'cj': cj_target,
            'altitude': capture_alt,
            'tof': tof,
            'dv': s1.dv_total_ms
        }

    except Exception as e:
        if verbose:
            print(f"   💥 Exception: {e}")
        return None

def main():
    print("=" * 80)
    print("🚀 ASTRO-TFC: DIRECT SEARCH WITH CANONICAL C_J")
    print("=" * 80)

    # 1. Get canonical C_J(L1)
    cj_l1_canonical, xl1 = get_l1_jacobi_canonical()

    print(f"🌍 L1 Point: x = {xl1:.6f}")
    print(f"⚡ C_J(L1) CANONICAL calculated: {cj_l1_canonical:.6f}")
    print(f"   Backend reports minimum: 3.188341")
    print(f"   Difference: {abs(cj_l1_canonical - 3.188341):.6f}")
    print()

    # 2. If our calculation differs from backend, use backend value
    if abs(cj_l1_canonical - 3.188341) > 0.1:
        print("⚠️  Using backend value (3.188341) as reference")
        cj_l1_canonical = 3.188341
        print(f"   C_J(L1) corrected: {cj_l1_canonical:.6f}")
        print()

    # 3. Test values below C_J(L1)
    #    Start conservative and then relax
    test_cj_values = [
        3.188,    # Right at L1 (should fail)
        3.187,    # Very close
        3.185,    # Close
        3.180,    # Moderate margin
        3.175,    # Wide margin
        3.170,    # Large margin
        3.160,    # Far below
        3.150,    # Extreme
    ]

    results = []

    for cj_test in test_cj_values:
        print(f"🎯 Testing C_J = {cj_test:.3f}...", end=" ")

        result = test_injection_direct(cj_test, verbose=False)

        if result:
            print(f"✅ Alt={result['altitude']:.0f}km, TOF={result['tof']:.1f}d, ΔV={result['dv']:.1f}m/s")
            results.append(result)
        else:
            print(f"❌")

    print()

    # 4. Show results
    if results:
        print("=" * 80)
        print("📊 SWEEP RESULTS")
        print("=" * 80)
        print(f"{'C_J':<10} {'Altitude':<12} {'TOF':<10} {'ΔV':<10}")
        print("-" * 42)

        best = None
        for r in results:
            print(f"{r['cj']:<10.3f} {r['altitude']:<12.0f} {r['tof']:<10.1f} {r['dv']:<10.1f}")

            if best is None or abs(r['altitude'] - TARGET_ALT) < abs(best['altitude'] - TARGET_ALT):
                if r['altitude'] > 0:
                    best = r

        if best:
            print()
            print("=" * 80)
            print("🏆 BEST TRAJECTORY")
            print("=" * 80)
            print(f"   Canonical C_J:  {best['cj']:.6f}")
            print(f"   Altitude:       {best['altitude']:.0f} km")
            print(f"   Flight time:    {best['tof']:.1f} days")
            print(f"   Injection ΔV:   {best['dv']:.1f} m/s")

            xenon, duration = total_mission(best['altitude'], best['tof'])
            print(f"\n📊 BUDGET:")
            print(f"   Xenon: {xenon:.1f} kg")
            print(f"   Total duration: {duration:.1f} days")
    else:
        print("❌ NO TRAJECTORY FOUND")
        print()
        print("🔍 This is strange. Check:")
        print("   1. Does the backend accept C_J=3.18?")
        print("   2. Is there any other rejection criteria?")
        print("   3. Try manually:")
        print("      config = astro_tfc.MissionConfig()")
        print("      config.target_jacobi = 3.18")
        print("      mission = astro_tfc.AstroTFCMission(config)")
        print("      s1 = mission.execute_stage_1()")

if __name__ == "__main__":
    main()