# test_high_fidelity.py
#!/usr/bin/env python3
import ctypes
import os
import sys

# ============================================================
# PRELOAD CSPICE BEFORE IMPORTING ASTRO_TFC
# ============================================================
def preload_cspice():
    """Loads libcspice.so into global space before importing astro_tfc"""
    # Search in several locations
    conda_prefix = os.environ.get("CONDA_PREFIX", os.path.expanduser("~/micromamba/envs/Astro"))
    lib_paths = [
        os.path.join(conda_prefix, "lib", "libcspice.so"),
        "/usr/local/lib/libcspice.so",
        "/usr/lib/libcspice.so",
        os.path.expanduser("~/micromamba/envs/Astro/lib/libcspice.so"),
    ]

    for lib_path in lib_paths:
        if os.path.exists(lib_path):
            try:
                ctypes.CDLL(lib_path, mode=ctypes.RTLD_GLOBAL)
                print(f"✅ CSPICE loaded from: {lib_path}")
                return True
            except Exception as e:
                print(f"   Failed {lib_path}: {e}")

    print("⚠️ libcspice.so not found. SPICE will not be available.")
    return False

# Preload before importing astro_tfc
preload_cspice()

# Now import astro_tfc
sys.path.insert(0, os.path.expanduser("~/rust"))
import astro_tfc