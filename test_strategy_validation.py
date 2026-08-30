#!/usr/bin/env python3
"""
test_strategy_validation.py - Validation compatible with ctypes
Does NOT break existing scripts
"""

import ctypes
import os
import sys
import glob

def main():
    print("🚀 MISSION STRATEGY VALIDATION (ctypes)")
    print("=" * 60)

    # 1. Load CSPICE first (same as test_high_fidelity_v4.py)
    print("\n📡 Loading CSPICE...")
    cspice_paths = [
        os.path.expanduser("~/micromamba/envs/Astro/lib/libcspice.so"),
        "/usr/lib/libcspice.so",
        "/usr/local/lib/libcspice.so",
    ]

    cspice_lib = None
    for path in cspice_paths:
        if os.path.exists(path):
            try:
                cspice_lib = ctypes.CDLL(path, mode=ctypes.RTLD_GLOBAL)
                print(f"✅ CSPICE loaded from: {path}")
                break
            except Exception as e:
                print(f"⚠️ Error with {path}: {e}")

    if not cspice_lib:
        print("❌ CSPICE not found")
        return 1

    # 2. Load libastro_tfc.so (same as before)
    print("\n🔧 Loading ASTRO-TFC...")
    lib_path = os.path.expanduser("~/rust/target/release/libastro_tfc.so")

    if not os.path.exists(lib_path):
        print(f"❌ Not found: {lib_path}")
        print("   Run first: cargo build --release")
        return 1

    try:
        lib = ctypes.CDLL(lib_path, mode=ctypes.RTLD_GLOBAL)
        print(f"✅ ASTRO-TFC loaded from: {lib_path}")
    except Exception as e:
        print(f"❌ Error loading ASTRO-TFC: {e}")
        return 1

    # 3. Configure the validation function
    # First verify it exists
    try:
        validate_func = lib.validate_all_strategies
    except AttributeError:
        print("❌ Function validate_all_strategies does not exist in the library")
        print("   You need to recompile with the function added")
        return 1

    validate_func.restype = ctypes.c_int
    validate_func.argtypes = []

    # 4. Execute validation
    print("\n🔧 Executing strategy validation...")
    print("=" * 60)
    result = validate_func()
    print("=" * 60)

    if result == 0:
        print("\n✅ Validation completed successfully")
        print("   Our strategies have been validated")
        print("   against real historical missions")
        return 0
    else:
        print("\n❌ Error during validation")
        return 1

if __name__ == "__main__":
    sys.exit(main())