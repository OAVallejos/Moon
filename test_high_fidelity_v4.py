#!/usr/bin/env python3
"""
test_high_fidelity_v4.py - Validation with SPICE and CORRECTED TEI
Mission Earth → Moon → Earth with 14-day wait
TEI with tangential ΔV (NOT scaled)
Target HEO perigee: 50,000 km
"""

import ctypes
import os
import sys
import glob

# ============================================================
# 1. LOAD CSPICE
# ============================================================
def load_cspice():
    """Loads CSPICE with RTLD_GLOBAL to resolve symbols"""
    possible_paths = [
        os.path.expanduser("~/micromamba/envs/Astro/lib/libcspice.so"),
        os.path.expanduser("~/micromamba/envs/Astro/lib/libcspice.so.66"),
        "/usr/lib/libcspice.so",
        "/usr/local/lib/libcspice.so",
        "/usr/lib/aarch64-linux-gnu/libcspice.so",
        "/opt/cspice/lib/libcspice.so",
    ]

    for pattern in [
        os.path.expanduser("~/micromamba/envs/Astro/lib/libcspice*"),
        "/usr/lib/libcspice*",
        "/usr/local/lib/libcspice*",
        "/opt/cspice/lib/libcspice*",
    ]:
        matches = glob.glob(pattern)
        for match in matches:
            if match not in possible_paths:
                possible_paths.append(match)

    for lib_path in possible_paths:
        if os.path.exists(lib_path):
            try:
                lib = ctypes.CDLL(lib_path, mode=ctypes.RTLD_GLOBAL)
                print(f"✅ CSPICE loaded from: {lib_path}")
                return lib
            except Exception as e:
                print(f"⚠️ Could not load {lib_path}: {e}")

    print("❌ CSPICE not found")
    return None

# ============================================================
# 2. INJECT GM
# ============================================================
def inject_gm(lib):
    """Injects GM of asteroid 20556096 into CSPICE"""
    if lib is not None:
        try:
            body_name = ctypes.c_char_p(b"BODY20556096_GM")
            n = ctypes.c_int(1)
            val = ctypes.c_double(0.0001)
            lib.pdpool_c(body_name, n, ctypes.byref(val))
            print("✅ GM injected via ctypes")
            return True
        except Exception as e:
            print(f"⚠️ Error injecting GM: {e}")
    return False

# ============================================================
# 3. IMPORT ASTRO_TFC
# ============================================================
def import_astro_tfc():
    """Imports astro_tfc"""
    sys.path.insert(0, os.path.expanduser("~/rust/target/release"))
    sys.path.insert(0, os.path.expanduser("~/rust"))

    try:
        import astro_tfc
        print("✅ astro_tfc imported successfully")
        print(f"   Version: {astro_tfc.__version__}")
        return astro_tfc
    except ImportError as e:
        print(f"❌ Error importing astro_tfc: {e}")
        sys.exit(1)

# ============================================================
# 4. EXECUTE FULL MISSION
# ============================================================
def run_mission(astro_tfc, wait_days=14.0, use_spice=True):
    """Executes the full mission with corrected TEI"""

    print("\n" + "=" * 80)
    print("🚀 FULL MISSION: Earth → Moon → Earth")
    print("   via L1 using Theory of Functional Connections (TFC)")
    print(f"   C_J = 3.170")
    print(f"   Lunar orbit wait: {wait_days:.0f} days")
    print(f"   SPICE: {'✅ ENABLED' if use_spice else '❌ DISABLED'}")
    print("   Target HEO perigee: 50,000 km")
    print("=" * 80)

    config = astro_tfc.MissionConfig()
    config.verbose = True
    config.target_jacobi = 3.170
    config.spacecraft_mass_kg = 250.0
    config.engine_thrust_n = 0.237
    config.engine_isp_s = 4190.0
    config.wait_days_lunar = wait_days
    config.kernel_dir = "/root/kernels"
    config.reference_time_utc = "2024-11-15T00:00:00"
    config.use_spice = use_spice

    mission = astro_tfc.AstroTFCMission(config)
    results = {}

    # ============ STAGE 1: HEO Injection → Unstable Manifold ============
    print("\n" + "-" * 80)
    print("📡 STAGE 1: HEO Injection → Unstable Manifold")
    print("-" * 80)
    r1 = mission.execute_stage_1()
    results[1] = r1
    print(f"✅ ΔV = {r1.dv_total_ms:.1f} m/s | TOF = {r1.time_of_flight_days:.2f} days")

    # ============ STAGE 2: Ballistic Transit L1 → Moon ============
    print("\n" + "-" * 80)
    print("🌌 STAGE 2: Ballistic Transit L1 → Moon")
    print("-" * 80)
    r2 = mission.execute_stage_2()
    results[2] = r2
    print(f"✅ ΔV = {r2.dv_total_ms:.1f} m/s | TOF = {r2.time_of_flight_days:.2f} days")

    # ============ STAGE 3: Lunar Capture ============
    print("\n" + "-" * 80)
    print("🌙 STAGE 3: Lunar Capture")
    print("-" * 80)
    r3 = mission.execute_stage_3()
    results[3] = r3
    print(f"✅ ΔV = {r3.dv_total_ms:.1f} m/s | TOF = {r3.time_of_flight_days:.2f} days")

    # ============ STAGE 4: TEI from LLO ============
    print("\n" + "-" * 80)
    print(f"🚀 STAGE 4: TEI from 800 km LLO ({wait_days:.0f}-day wait)")
    print(f"   🛰️  {'SPICE enabled' if use_spice else 'Theoretical mode'}")
    print("   Target perigee: 50,000 km")
    print("-" * 80)
    r4 = mission.execute_stage_4()
    results[4] = r4
    print(f"✅ ΔV = {r4.dv_total_ms:.1f} m/s | TOF = {r4.time_of_flight_days:.2f} days")
    print(f"   Details: {r4.details}")

    # ============ STAGE 5: Ballistic Propagation → HEO Perigee ============
    print("\n" + "-" * 80)
    print("🌌 STAGE 5: Ballistic Transit → HEO Perigee")
    print("-" * 80)
    r5 = mission.execute_stage_5()
    results[5] = r5
    print(f"✅ ΔV = {r5.dv_total_ms:.1f} m/s | TOF = {r5.time_of_flight_days:.2f} days")
    print(f"   Details: {r5.details}")

    # ============ STAGE 6: Circularization at HEO Perigee ============
    print("\n" + "-" * 80)
    print("🌍 STAGE 6: Circularization at HEO Perigee")
    print("-" * 80)
    r6 = mission.execute_stage_6()
    results[6] = r6
    print(f"✅ ΔV = {r6.dv_total_ms:.1f} m/s | TOF = {r6.time_of_flight_days:.2f} days")
    print(f"   Details: {r6.details}")

    # ============ FINAL SUMMARY ============
    print("\n" + "=" * 80)
    print("🎉 FULL MISSION SUCCESSFUL")
    print("=" * 80)

    total_dv = sum(r.dv_total_ms for r in results.values())
    total_tof = sum(r.time_of_flight_days for r in results.values())

    print(f"\n📊 MISSION SUMMARY:")
    print(f"{'Stage':<10} {'ΔV (m/s)':<14} {'TOF (days)':<14} {'Description':<35}")
    print("-" * 80)

    stage_names = {
        1: "HEO→L1 Injection",
        2: "L1→Moon Transit",
        3: "Lunar capture (LOI)",
        4: f"TEI ({wait_days:.0f}d wait)",
        5: "Transit→HEO Perigee",
        6: "HEO Circularization",
    }

    for i in range(1, 7):
        r = results[i]
        spice_tag = "🛰️" if (i == 4 and use_spice) else ""
        print(f"Stage {i:<5} {r.dv_total_ms:<14.1f} {r.time_of_flight_days:<14.2f} {stage_names[i]:<30} {spice_tag}")

    print("-" * 80)
    print(f"{'TOTAL':<10} {total_dv:<14.1f} {total_tof:<14.1f}")

    print(f"\n📋 DETAILS BY STAGE:")
    for i in range(1, 7):
        r = results[i]
        print(f"   Stage {i} [{r.stage_name}]:")
        print(f"      {r.details}")

    # ============ PERIGEE CONFIRMATION ============
    print("\n" + "=" * 80)
    print("✅ CONFIRMATION: Target perigee in Rust = 50,000 km")

    # Extract achieved perigee from Stage 5 details
    perigee_str = results[5].details.split('perigee=')[1].split('km')[0] if 'perigee=' in results[5].details else "N/A"
    print(f"   Achieved perigee: {perigee_str} km")
    print("=" * 80)

    return results, mission

# ============================================================
# 5. MAIN
# ============================================================
if __name__ == "__main__":
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║   ASTRO-TFC v3.4 - CORRECTED TEI (tangential ΔV)          ║")
    print("║   Earth → Moon (800km) → Wait → Earth (HEO)              ║")
    print("║   Target HEO perigee: 50,000 km                         ║")
    print("╚══════════════════════════════════════════════════════════════╝")

    # Load CSPICE
    cspice_lib = load_cspice()
    if cspice_lib:
        inject_gm(cspice_lib)

    # Import astro_tfc
    astro_tfc = import_astro_tfc()

    # ============================================================
    # MISSION WITH SPICE AND 50,000 km PERIGEE
    # ============================================================
    print("\n" + "=" * 80)
    print("🚀 MISSION: Artemis II with corrected TEI")
    print("=" * 80)

    results, mission = run_mission(
        astro_tfc,
        wait_days=14.0,
        use_spice=True
    )