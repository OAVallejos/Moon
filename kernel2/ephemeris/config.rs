//! ephemeris/config.rs - Configuración de efemérides

use super::spice::SpiceKernelConfig;
use crate::integration::IntegratorConfig;

/// Modo de operación de efemérides
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum EphemerisMode {
    Spice,
    Analytical,
    Hybrid,
}

/// Configuración de efemérides (sin duplicados)
#[derive(Clone)]
pub struct EphemerisConfig {
    pub mode: EphemerisMode,
    pub spice: SpiceKernelConfig,
    pub integrator: IntegratorConfig,
    pub bodies: Vec<String>,
    pub frame: String,
    pub include_aberration: bool,
    pub scale_factor: f64,
    pub verbose: bool,
}

impl Default for EphemerisConfig {
    fn default() -> Self {
        EphemerisConfig {
            mode: EphemerisMode::Spice,
            spice: SpiceKernelConfig::default(),
            integrator: IntegratorConfig::default(),
            bodies: vec![
                "SUN".to_string(),
                "EARTH".to_string(),
                "MOON".to_string(),
                "MARS".to_string(),
                "JUPITER".to_string(),
                "VENUS".to_string(),
            ],
            frame: "J2000".to_string(),
            include_aberration: true,
            scale_factor: 1.0,
            verbose: false,
        }
    }
}
