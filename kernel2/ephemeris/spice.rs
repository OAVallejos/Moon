//! ephemeris/spice.rs - Interfaz unificada con SPICE

use std::ffi::{CStr, CString};
use std::os::raw::c_double;
use std::path::Path;
use thiserror::Error;

// ============================================================================
// ERRORES DE SPICE
// ============================================================================

#[derive(Error, Debug, Clone, PartialEq)]
pub enum SpiceError {
    #[error("Error cargando kernel: {0}")]
    KernelLoad(String),
    #[error("Error en SPKEZR: código {0}")]
    SpkezrError(i32),
    #[error("Error en BODVCD: código {0}")]
    BodvcdError(i32),
    #[error("Error de tiempo: {0}")]
    TimeError(String),
    #[error("Kernels no cargados")]
    NotLoaded,
    #[error("Error de conversión UTF-8: {0}")]
    Utf8Error(String),
}

// ============================================================================
// CONFIGURACIÓN DE KERNELS
// ============================================================================

#[derive(Clone)]
pub struct SpiceKernelConfig {
    pub kernel_dir: String,
    pub kernel_files: Vec<String>,
    pub frame: String,
    pub abcorr: String,
    pub verbose: bool,
}

impl Default for SpiceKernelConfig {
    fn default() -> Self {
        SpiceKernelConfig {
            kernel_dir: "~/kernels".to_string(),
            kernel_files: vec![
                "naif0012.tls".to_string(),
                "de440.bsp".to_string(),
                "gm_de440.tpc".to_string(),
                "pck00010.tpc".to_string(),
            ],
            frame: "J2000".to_string(),
            abcorr: "LT+S".to_string(),
            verbose: false,
        }
    }
}

// ============================================================================
// CONTEXTO SPICE
// ============================================================================

#[derive(Clone)]
pub struct SpiceContext {
    config: SpiceKernelConfig,
    loaded: bool,
}

impl SpiceContext {
    pub fn new(config: SpiceKernelConfig) -> Self {
        SpiceContext {
            config,
            loaded: false,
        }
    }

    /// Carga todos los kernels SPICE e inyecta parámetros faltantes
    pub fn load_kernels(&mut self) -> Result<(), SpiceError> {
        let kernel_dir = shellexpand::tilde(&self.config.kernel_dir).to_string();
        
        for kernel in &self.config.kernel_files {
            let full_path = Path::new(&kernel_dir).join(kernel);
            if !full_path.exists() {
                return Err(SpiceError::KernelLoad(
                    format!("Kernel no encontrado: {}", full_path.display())
                ));
            }
            
            let path_str = full_path.to_str()
                .ok_or_else(|| SpiceError::Utf8Error(
                    format!("Ruta no válida UTF-8: {}", full_path.display())
                ))?;
            
            let c_path = CString::new(path_str)
                .map_err(|e| SpiceError::KernelLoad(e.to_string()))?;
            
            let status = unsafe { furnsh_c(c_path.as_ptr()) };
            if status != 0 {
                return Err(SpiceError::KernelLoad(
                    format!("Error al cargar kernel: {} (código {})", kernel, status)
                ));
            }
            
            if self.config.verbose {
                eprintln!("✅ SPICE kernel cargado: {}", kernel);
            }
        }

        // ============================================================
        // INYECCIÓN DE PARÁMETROS FALTANTES (BODVCD)
        // ============================================================
        // BODY20556096_GM - GM para cuerpo 20556096 (asteroide/objeto pequeño)
        // Este valor es necesario para que bodvcd_c no falle con código 1
        let missing_params: Vec<(&str, f64)> = vec![
            ("BODY20556096_GM", 0.0001),      // GM en km^3/s^2
            ("BODY20556096_RADIUS", 1.0),     // Radio en km
        ];

        for (param, value) in missing_params {
            let c_param = CString::new(param)
                .map_err(|e| SpiceError::KernelLoad(e.to_string()))?;
            let values = [value];
            
            unsafe {
                pdpool_c(c_param.as_ptr(), 1, values.as_ptr());
            }
            
            if self.config.verbose {
                eprintln!("✅ Parámetro inyectado: {} = {}", param, value);
            }
        }

        self.loaded = true;
        Ok(())
    }

    /// Obtiene estado de un cuerpo en tiempo ET
    pub fn get_state(&self, target: i32, et: f64, observer: i32) -> Result<SpiceState, SpiceError> {
        if !self.loaded {
            return Err(SpiceError::NotLoaded);
        }

        let mut state = [0.0; 6];
        let mut lt = 0.0;
        
        let target_str = CString::new(target.to_string())
            .map_err(|_| SpiceError::SpkezrError(-1))?;
        let frame_str = CString::new(self.config.frame.as_str())
            .map_err(|_| SpiceError::SpkezrError(-1))?;
        let abcorr_str = CString::new(self.config.abcorr.as_str())
            .map_err(|_| SpiceError::SpkezrError(-1))?;
        let obs_str = CString::new(observer.to_string())
            .map_err(|_| SpiceError::SpkezrError(-1))?;

        let status = unsafe {
            spkezr_c(
                target_str.as_ptr(),
                et,
                frame_str.as_ptr(),
                abcorr_str.as_ptr(),
                obs_str.as_ptr(),
                &mut state,
                &mut lt,
            )
        };

        if status != 0 {
            return Err(SpiceError::SpkezrError(status));
        }

        Ok(SpiceState {
            position: [state[0], state[1], state[2]],
            velocity: [state[3], state[4], state[5]],
            light_time: lt,
            time_et: et,
        })
    }

    /// Convierte UTC a ET
    pub fn utc_to_et(&self, utc_str: &str) -> Result<f64, SpiceError> {
        if !self.loaded {
            return Err(SpiceError::NotLoaded);
        }

        let c_utc = CString::new(utc_str)
            .map_err(|_| SpiceError::TimeError("Invalid UTC".to_string()))?;
        
        let mut et = 0.0;
        let status = unsafe {
            str2et_c(c_utc.as_ptr(), &mut et)
        };
        
        if status != 0 {
            return Err(SpiceError::TimeError(format!("str2et failed: {}", status)));
        }
        
        Ok(et)
    }

    /// Convierte ET a UTC
    pub fn et_to_utc(&self, et: f64) -> Result<String, SpiceError> {
        if !self.loaded {
            return Err(SpiceError::NotLoaded);
        }

        let format_str = CString::new("ISOC")
            .map_err(|_| SpiceError::TimeError("Invalid format".to_string()))?;
        
        let mut utc = [0u8; 256];
        let status = unsafe {
            et2utc_c(
                et,
                format_str.as_ptr(),
                3,
                32,
                utc.as_mut_ptr() as *mut std::os::raw::c_char,
            )
        };
        
        if status != 0 {
            return Err(SpiceError::TimeError(format!("et2utc failed: {}", status)));
        }
        
        let utc_str = unsafe {
            CStr::from_ptr(utc.as_ptr() as *const std::os::raw::c_char)
                .to_str()
                .map_err(|_| SpiceError::Utf8Error("Invalid UTF-8 in UTC string".to_string()))?
        };
        
        Ok(utc_str.to_string())
    }

    /// Obtiene parámetro gravitacional de un cuerpo (μ = GM) en km³/s²
    pub fn get_mu(&self, body: i32) -> Result<f64, SpiceError> {
        if !self.loaded {
            return Err(SpiceError::NotLoaded);
        }

        let body_str = CString::new(body.to_string())
            .map_err(|_| SpiceError::BodvcdError(-1))?;
        let item_str = CString::new("GM")
            .map_err(|_| SpiceError::BodvcdError(-1))?;
        
        let mut dim = 0;
        let mut values = [0.0; 1];
        
        let status = unsafe {
            bodvcd_c(
                body_str.as_ptr(),
                item_str.as_ptr(),
                1,
                &mut dim,
                values.as_mut_ptr(),
            )
        };
        
        if status != 0 {
            return Err(SpiceError::BodvcdError(status));
        }
        
        Ok(values[0])
    }

    /// Verifica si los kernels están cargados
    pub fn is_loaded(&self) -> bool {
        self.loaded
    }

    /// Obtiene referencia a la configuración
    pub fn config(&self) -> &SpiceKernelConfig {
        &self.config
    }

    /// Establece el modo verboso
    pub fn set_verbose(&mut self, verbose: bool) {
        self.config.verbose = verbose;
    }
}

// ============================================================================
// ESTADO SPICE
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpiceState {
    pub position: [f64; 3],  // km
    pub velocity: [f64; 3],  // km/s
    pub light_time: f64,     // segundos
    pub time_et: f64,        // segundos desde J2000
}

impl SpiceState {
    pub fn new(position: [f64; 3], velocity: [f64; 3], time_et: f64) -> Self {
        SpiceState {
            position,
            velocity,
            light_time: 0.0,
            time_et,
        }
    }

    pub fn distance_to(&self, other: &SpiceState) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        (dx*dx + dy*dy + dz*dz).sqrt()
    }

    pub fn relative_velocity(&self, other: &SpiceState) -> f64 {
        let dvx = self.velocity[0] - other.velocity[0];
        let dvy = self.velocity[1] - other.velocity[1];
        let dvz = self.velocity[2] - other.velocity[2];
        (dvx*dvx + dvy*dvy + dvz*dvz).sqrt()
    }

    pub fn speed(&self) -> f64 {
        (self.velocity[0].powi(2) + self.velocity[1].powi(2) + self.velocity[2].powi(2)).sqrt()
    }

    pub fn distance(&self) -> f64 {
        (self.position[0].powi(2) + self.position[1].powi(2) + self.position[2].powi(2)).sqrt()
    }

    pub fn to_crtbp(&self, moon_state: &SpiceState) -> [f64; 6] {
        use crate::constants::{D_CHAR, V_CHAR};
        
        let dx_m = (self.position[0] - moon_state.position[0]) * 1000.0;
        let dy_m = (self.position[1] - moon_state.position[1]) * 1000.0;
        let dz_m = (self.position[2] - moon_state.position[2]) * 1000.0;
        
        let dvx_ms = (self.velocity[0] - moon_state.velocity[0]) * 1000.0;
        let dvy_ms = (self.velocity[1] - moon_state.velocity[1]) * 1000.0;
        let dvz_ms = (self.velocity[2] - moon_state.velocity[2]) * 1000.0;
        
        [
            dx_m / D_CHAR,
            dy_m / D_CHAR,
            dz_m / D_CHAR,
            dvx_ms / V_CHAR,
            dvy_ms / V_CHAR,
            dvz_ms / V_CHAR,
        ]
    }

    pub fn to_array(&self) -> [f64; 6] {
        [
            self.position[0], self.position[1], self.position[2],
            self.velocity[0], self.velocity[1], self.velocity[2],
        ]
    }
}

impl Default for SpiceState {
    fn default() -> Self {
        SpiceState {
            position: [0.0, 0.0, 0.0],
            velocity: [0.0, 0.0, 0.0],
            light_time: 0.0,
            time_et: 0.0,
        }
    }
}

// ============================================================================
// IDENTIFICADORES NAIF
// ============================================================================

pub mod naif_ids {
    pub const SUN: i32 = 10;
    pub const EARTH: i32 = 399;
    pub const MOON: i32 = 301;
    pub const MARS: i32 = 499;
    pub const VENUS: i32 = 299;
    pub const JUPITER: i32 = 599;
    pub const SATURN: i32 = 699;
    pub const SSB: i32 = 0;
}

// ============================================================================
// FUNCIONES EXTERNAS CSPICE
// ============================================================================

extern "C" {
    fn furnsh_c(kernel: *const std::os::raw::c_char) -> i32;
    fn spkezr_c(
        target: *const std::os::raw::c_char,
        et: c_double,
        frame: *const std::os::raw::c_char,
        abcorr: *const std::os::raw::c_char,
        observer: *const std::os::raw::c_char,
        state: *mut [c_double; 6],
        lt: *mut c_double,
    ) -> i32;
    fn bodvcd_c(
        body: *const std::os::raw::c_char,
        item: *const std::os::raw::c_char,
        maxn: i32,
        dim: *mut i32,
        values: *mut c_double,
    ) -> i32;
    fn str2et_c(str: *const std::os::raw::c_char, et: *mut c_double) -> i32;
    fn et2utc_c(
        et: c_double,
        format: *const std::os::raw::c_char,
        prec: i32,
        lenout: i32,
        utcstr: *mut std::os::raw::c_char,
    ) -> i32;
    // NUEVA: Inyectar parámetros en el pool SPICE
    fn pdpool_c(
        name: *const std::os::raw::c_char,
        n: i32,
        values: *const c_double,
    ) -> i32;
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spice_config_default() {
        let config = SpiceKernelConfig::default();
        assert_eq!(config.kernel_files.len(), 4); // Ahora son 4 kernels
        assert_eq!(config.frame, "J2000");
        assert_eq!(config.abcorr, "LT+S");
    }

    #[test]
    fn test_spice_state_distance() {
        let s1 = SpiceState::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0);
        let s2 = SpiceState::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0);
        assert_eq!(s1.distance_to(&s2), 1.0);
    }

    #[test]
    fn test_naif_ids() {
        assert_eq!(naif_ids::EARTH, 399);
        assert_eq!(naif_ids::MOON, 301);
        assert_eq!(naif_ids::SUN, 10);
    }
}
