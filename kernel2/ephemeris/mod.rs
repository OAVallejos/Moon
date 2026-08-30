//! ephemeris/mod.rs - Sistema completo de efemérides para ASTRO-TFC

pub mod config;
pub mod spice;
pub mod integrator;
pub mod planner;
pub mod validation;

// Re-exportar tipos principales
pub use config::{EphemerisConfig, EphemerisMode};
pub use spice::{SpiceContext, SpiceState, SpiceError, SpiceKernelConfig};
pub use integrator::{NBodyIntegrator, NBodyState};
pub use planner::{LaunchWindowPlanner, LaunchWindow};

// CORREGIDO: Usar los nombres correctos de validation.rs
pub use validation::{
    MissionStrategyValidator,
    StrategyValidation,
    MissionBenchmark,
    MissionStrategyType,
    Smart1Benchmark,
    GrailBenchmark,
    CapstoneBenchmark,
    validate_all_strategies,
};

use crate::integration::IntegratorConfig;

// ============================================================================
// CONFIGURACIÓN (compatibilidad)
// ============================================================================

/// Alias para compatibilidad con código existente
pub type EphemerisFullConfig = EphemerisConfig;

// ============================================================================
// MODELO ANALÍTICO
// ============================================================================

use std::collections::HashMap;

pub struct AnalyticalModel {
    mu_table: HashMap<String, f64>,
}

impl AnalyticalModel {
    pub fn new() -> Self {
        let mut mu_table = HashMap::new();
        mu_table.insert("SUN".to_string(), 1.32712440018e20);
        mu_table.insert("EARTH".to_string(), 3.986004418e14);
        mu_table.insert("MOON".to_string(), 4.902800118e12);
        mu_table.insert("JUPITER".to_string(), 1.26686534e17);
        mu_table.insert("VENUS".to_string(), 3.24858592e14);
        mu_table.insert("MARS".to_string(), 4.28283758e13);

        AnalyticalModel { mu_table }
    }

    pub fn get_mu(&self, body: &str) -> Option<f64> {
        self.mu_table.get(body).copied()
    }

    pub fn solve_kepler(&self, M: f64, e: f64) -> f64 {
        let mut E = M;
        for _ in 0..100 {
            let dE = (E - e * E.sin() - M) / (1.0 - e * E.cos());
            E -= dE;
            if dE.abs() < 1e-12 {
                break;
            }
        }
        E
    }

    pub fn earth_state(&self, t: f64) -> spice::SpiceState {
        let a: f64 = 1.495978707e11;
        let e: f64 = 0.0167086;
        let omega: f64 = 2.0 * std::f64::consts::PI / 365.25;
        let long_peri: f64 = 1.796;

        let M = omega * t / 86400.0;
        let E = self.solve_kepler(M, e);
        let nu = 2.0 * ((1.0 + e)/(1.0 - e)).sqrt().atan2((E/2.0).tan());

        let r = a * (1.0 - e * E.cos());

        spice::SpiceState {
            position: [r * nu.cos(), r * nu.sin(), 0.0],
            velocity: [
                -r * nu.sin() * omega / 86400.0,
                r * nu.cos() * omega / 86400.0,
                0.0,
            ],
            light_time: 0.0,
            time_et: t,
        }
    }

    pub fn moon_state(&self, t: f64) -> spice::SpiceState {
        let earth = self.earth_state(t);

        let a: f64 = 3.844e8;
        let e: f64 = 0.0549;
        let omega: f64 = 2.0 * std::f64::consts::PI / 27.32166;
        let incl: f64 = 0.0898;
        let long_peri: f64 = 1.627;

        let M = omega * t / 86400.0;
        let E = self.solve_kepler(M, e);
        let nu = 2.0 * ((1.0 + e)/(1.0 - e)).sqrt().atan2((E/2.0).tan());

        let r = a * (1.0 - e * E.cos());

        let x_orbit = r * nu.cos();
        let y_orbit = r * nu.sin();

        let x = x_orbit * long_peri.cos() - y_orbit * long_peri.sin();
        let y = (x_orbit * long_peri.sin() + y_orbit * long_peri.cos()) * incl.cos();
        let z = (x_orbit * long_peri.sin() + y_orbit * long_peri.cos()) * incl.sin();

        spice::SpiceState {
            position: [earth.position[0] + x, earth.position[1] + y, earth.position[2] + z],
            velocity: earth.velocity,
            light_time: 0.0,
            time_et: t,
        }
    }

    pub fn sun_state(&self, _t: f64) -> spice::SpiceState {
        spice::SpiceState {
            position: [0.0, 0.0, 0.0],
            velocity: [0.0, 0.0, 0.0],
            light_time: 0.0,
            time_et: _t,
        }
    }

    pub fn mars_state(&self, t: f64) -> spice::SpiceState {
        let a: f64 = 2.279e11;
        let e: f64 = 0.0934;
        let omega: f64 = 2.0 * std::f64::consts::PI / 686.97;

        let M = omega * t / 86400.0;
        let E = self.solve_kepler(M, e);
        let nu = 2.0 * ((1.0 + e)/(1.0 - e)).sqrt().atan2((E/2.0).tan());
        let r = a * (1.0 - e * E.cos());

        spice::SpiceState {
            position: [r * nu.cos(), r * nu.sin(), 0.0],
            velocity: [0.0, 0.0, 0.0],
            light_time: 0.0,
            time_et: t,
        }
    }

    pub fn get_state(&self, body: &str, t: f64) -> Option<spice::SpiceState> {
        match body {
            "EARTH" => Some(self.earth_state(t)),
            "MOON" => Some(self.moon_state(t)),
            "SUN" => Some(self.sun_state(t)),
            "MARS" => Some(self.mars_state(t)),
            _ => None,
        }
    }
}

// ============================================================================
// SIMULADOR DE N-CUERPOS
// ============================================================================

pub struct NBodySimulator {
    config: EphemerisConfig,
    spice: Option<SpiceContext>,
    analytical: AnalyticalModel,
    body_cache: HashMap<String, spice::SpiceState>,
}

impl NBodySimulator {
    pub fn new(config: EphemerisConfig) -> Self {
        NBodySimulator {
            config: config.clone(),
            spice: None,
            analytical: AnalyticalModel::new(),
            body_cache: HashMap::new(),
        }
    }

    pub fn load_spice_kernels(&mut self) -> Result<(), SpiceError> {
        if self.config.mode == EphemerisMode::Analytical {
            return Ok(());
        }

        let mut spice = SpiceContext::new(self.config.spice.clone());
        spice.load_kernels()?;
        self.spice = Some(spice);
        Ok(())
    }

    pub fn get_body_state(&mut self, body: &str, t: f64) -> Result<spice::SpiceState, SpiceError> {
        // Verificar caché
        if let Some(cached) = self.body_cache.get(body) {
            if (cached.time_et - t).abs() < 60.0 {
                return Ok(*cached);
            }
        }

        let state = match self.config.mode {
            EphemerisMode::Spice => {
                if let Some(ref spice) = self.spice {
                    let body_id = self.body_name_to_id(body)?;
                    spice.get_state(body_id, t, 0)?
                } else {
                    return Err(SpiceError::NotLoaded);
                }
            }
            EphemerisMode::Analytical => {
                self.analytical.get_state(body, t)
                    .ok_or_else(|| SpiceError::KernelLoad(format!("Cuerpo no encontrado: {}", body)))?
            }
            EphemerisMode::Hybrid => {
                match body {
                    "EARTH" | "MOON" | "SUN" => {
                        if let Some(ref spice) = self.spice {
                            let body_id = self.body_name_to_id(body)?;
                            spice.get_state(body_id, t, 0)?
                        } else {
                            self.analytical.get_state(body, t)
                                .ok_or_else(|| SpiceError::KernelLoad(format!("Cuerpo no encontrado: {}", body)))?
                        }
                    }
                    _ => {
                        self.analytical.get_state(body, t)
                            .ok_or_else(|| SpiceError::KernelLoad(format!("Cuerpo no encontrado: {}", body)))?
                    }
                }
            }
        };

        self.body_cache.insert(body.to_string(), state);
        Ok(state)
    }

    fn body_name_to_id(&self, body: &str) -> Result<i32, SpiceError> {
        match body {
            "SUN" => Ok(10),
            "EARTH" => Ok(399),
            "MOON" => Ok(301),
            "MARS" => Ok(499),
            "JUPITER" => Ok(599),
            "VENUS" => Ok(299),
            _ => Err(SpiceError::KernelLoad(format!("ID no encontrado: {}", body))),
        }
    }

    pub fn nbody_acceleration(
        &mut self,
        spacecraft_pos: &[f64; 3],
        time_utc: f64,
    ) -> Result<[f64; 3], SpiceError> {
        let mut acc = [0.0; 3];

        // Usar un clon de los bodies para evitar problemas de borrow
        let bodies: Vec<String> = self.config.bodies.clone();

        for body in &bodies {
            let state = self.get_body_state(body, time_utc)?;
            let mu = if let Some(ref spice) = self.spice {
                let id = self.body_name_to_id(body)?;
                spice.get_mu(id)?
            } else {
                self.analytical.get_mu(body)
                    .ok_or_else(|| SpiceError::KernelLoad(format!("μ no encontrado: {}", body)))?
            };

            let dx = spacecraft_pos[0] - state.position[0];
            let dy = spacecraft_pos[1] - state.position[1];
            let dz = spacecraft_pos[2] - state.position[2];
            let r = (dx*dx + dy*dy + dz*dz).sqrt();

            if r < 1.0 {
                continue;
            }

            let g = mu / (r * r * r);
            acc[0] += g * dx;
            acc[1] += g * dy;
            acc[2] += g * dz;
        }

        Ok(acc)
    }
}

// ============================================================================
// FUNCIÓN DE ALTO NIVEL
// ============================================================================

pub fn simulate_mission_with_ephemeris(
    initial_state: [f64; 6],
    launch_time_utc: f64,
    duration_days: f64,
    config: &EphemerisConfig,
) -> Result<Vec<NBodyState>, SpiceError> {
    use crate::integration::IntegratorConfig;

    let mut simulator = NBodySimulator::new(config.clone());
    simulator.load_spice_kernels()?;

    let integrator_config = IntegratorConfig {
        rtol: 1e-12,
        atol: 1e-15,
        max_step: 3600.0 / 86400.0,
        ..IntegratorConfig::default()
    };

    let mut integrator = NBodyIntegrator::new(
        simulator.spice.ok_or(SpiceError::NotLoaded)?,
        integrator_config,
        config.bodies.clone(),
    );

    let duration_seconds = duration_days * 86400.0;
    let max_steps = (duration_seconds / 3600.0) as usize + 10;

    integrator.propagate(
        &initial_state,
        launch_time_utc,
        duration_seconds,
        max_steps,
    )
}
