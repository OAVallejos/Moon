//! Constantes físicas del sistema Tierra-Luna para CRTBP adimensional
//! Basado en Almeida Jr. et al. (2026), Sección 4.2

/// Parámetro gravitacional de la Tierra [m³/s²]
pub const MU_E: f64 = 3.975837768911438e14;

/// Parámetro gravitacional de la Luna [m³/s²]
pub const MU_M: f64 = 4.890329364450684e12;

/// Distancia media Tierra-Luna [m]
pub const L_DISTANCE: f64 = 3.84405000e8;

/// Radio medio de la Tierra [m]
pub const R_EARTH: f64 = 6.378e6;

/// Radio medio de la Luna [m]
pub const R_MOON: f64 = 1.738e6;

/// Velocidad angular del sistema [rad/s] = sqrt((μ_E + μ_M) / L³)
pub const OMEGA: f64 = 2.661861e-6;

/// Tiempo característico para adimensionalización [s] = 1/ω
pub const T_CHAR: f64 = 1.0 / OMEGA; // ~4.36 días

/// Distancia característica = L [m]
pub const D_CHAR: f64 = L_DISTANCE;

/// Velocidad característica = L*ω [m/s]
pub const V_CHAR: f64 = L_DISTANCE * OMEGA; // ~1023.5 m/s

/// Distancia de la Tierra al baricentro [m]
pub const D1: f64 = L_DISTANCE * MU_M / (MU_E + MU_M);

/// Distancia de la Luna al baricentro [m]
pub const D2: f64 = L_DISTANCE * MU_E / (MU_E + MU_M);

/// Parámetro de masa normalizado μ = M_Luna/(M_Tierra + M_Luna)
pub const MU_NORMALIZED: f64 = MU_M / (MU_E + MU_M);