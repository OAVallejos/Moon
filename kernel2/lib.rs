//! lib4 astro_tfc: Simulation engine for Earth-Moon transfers via L1
//! using the Theory of Functional Connections (TFC).
//!
//! VERSIÓN 3.5: Integración SPICE + Low Thrust + Validación de Estrategias

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use pyo3::types::PyDict;

mod constants;
mod crtbp;
mod tfc;
mod geocentric;
mod lunar;
mod manifold;
mod integration;
mod r#return;
mod propulsion;

pub mod ephemeris;
pub mod low_thrust;

use constants::{
    MU_E, MU_M, L_DISTANCE, R_EARTH, R_MOON, OMEGA, D1, D2,
    T_CHAR, D_CHAR, V_CHAR, MU_NORMALIZED,
};
use crtbp::{StateVector, crtbp_derivatives, jacobi_constant};
use geocentric::{OrbitalState, compute_injection};
use lunar::{CaptureConfig, analyze_capture, CaptureDetection};
use manifold::{
    l1_position, unstable_eigenvalue, eigenvector_y_factor, ManifoldType,
};
use integration::{
    EarthMoon, IntegratorConfig,
    propagate_manifold, lunar_altitude_km_from_dimensional,
};
use r#return::{
    ReturnConfig,
    compute_tei_from_llo,
    propagate_ballistic_to_earth_perigee,
    analyze_earth_capture,
    TeiResult,
    ReturnPropagationResult,
    EarthCapture,
};
use propulsion::{IonEngineConfig, compute_mission_budget};

// ============================================================================
// UNIT CONVERSION LAYER
// ============================================================================

fn normalize_state(dimensional: &[f64]) -> StateVector {
    StateVector::new(
        dimensional[0] / D_CHAR, dimensional[1] / D_CHAR, dimensional[2] / D_CHAR,
        dimensional[3] / V_CHAR, dimensional[4] / V_CHAR, dimensional[5] / V_CHAR,
    )
}

fn denormalize_array(canonical: &[f64; 6]) -> Vec<f64> {
    vec![
        canonical[0] * D_CHAR, canonical[1] * D_CHAR, canonical[2] * D_CHAR,
        canonical[3] * V_CHAR, canonical[4] * V_CHAR, canonical[5] * V_CHAR,
    ]
}

// ============================================================================
// PYTHON CLASSES
// ============================================================================

#[pyclass]
#[derive(Clone)]
struct MissionConfig {
    #[pyo3(get, set)]
    target_jacobi: f64,
    #[pyo3(get, set)]
    integration_rtol: f64,
    #[pyo3(get, set)]
    max_iterations: usize,
    #[pyo3(get, set)]
    verbose: bool,
    #[pyo3(get, set)]
    spacecraft_mass_kg: f64,
    #[pyo3(get, set)]
    engine_thrust_n: f64,
    #[pyo3(get, set)]
    engine_isp_s: f64,
    #[pyo3(get, set)]
    wait_days_lunar: f64,
    
    // ============================================================
    // NUEVOS: Configuración SPICE
    // ============================================================
    #[pyo3(get, set)]
    kernel_dir: String,
    #[pyo3(get, set)]
    reference_time_utc: String,
    #[pyo3(get, set)]
    use_spice: bool,
}

#[pymethods]
impl MissionConfig {
    #[new]
    fn new() -> Self {
        MissionConfig {
            target_jacobi: 3.170,
            integration_rtol: 1e-12,
            max_iterations: 500,
            verbose: false,
            spacecraft_mass_kg: 250.0,
            engine_thrust_n: 0.08,
            engine_isp_s: 1500.0,
            wait_days_lunar: 14.0,
            // Valores por defecto para SPICE
            kernel_dir: "~/kernels".to_string(),
            reference_time_utc: "2026-01-01T00:00:00".to_string(),
            use_spice: false,
        }
    }
}

#[pyclass]
#[derive(Clone)]
struct StageResult {
    #[pyo3(get)]
    stage_name: String,
    #[pyo3(get)]
    dv_total_ms: f64,
    #[pyo3(get)]
    time_of_flight_days: f64,
    #[pyo3(get)]
    final_state: Vec<f64>,
    #[pyo3(get)]
    canonical_state: Vec<f64>,
    #[pyo3(get)]
    jacobi_error: f64,
    #[pyo3(get)]
    success: bool,
    #[pyo3(get)]
    details: String,
}

#[pymethods]
impl StageResult {
    fn to_dict(&self, py: Python) -> PyResult<PyObject> {
        let dict = PyDict::new(py);
        dict.set_item("stage_name", &self.stage_name)?;
        dict.set_item("dv_total_ms", self.dv_total_ms)?;
        dict.set_item("time_of_flight_days", self.time_of_flight_days)?;
        dict.set_item("final_state", &self.final_state)?;
        dict.set_item("jacobi_error", self.jacobi_error)?;
        dict.set_item("success", self.success)?;
        dict.set_item("details", &self.details)?;
        Ok(dict.into())
    }
}

// ============================================================================
// ASTROTFC MISSION
// ============================================================================

#[pyclass]
struct AstroTFCMission {
    config: MissionConfig,
    #[pyo3(get)]
    stage_1: Option<StageResult>,
    #[pyo3(get)]
    stage_2: Option<StageResult>,
    #[pyo3(get)]
    stage_3: Option<StageResult>,
    #[pyo3(get)]
    stage_4: Option<StageResult>,
    #[pyo3(get)]
    stage_5: Option<StageResult>,
    #[pyo3(get)]
    stage_6: Option<StageResult>,
    #[pyo3(get)]
    total_dv_ms: f64,
    #[pyo3(get)]
    total_tof_days: f64,
    #[pyo3(get)]
    mission_complete: bool,
}

#[pymethods]
impl AstroTFCMission {
    #[new]
    fn new(config: MissionConfig) -> Self {
        AstroTFCMission {
            config,
            stage_1: None, stage_2: None, stage_3: None,
            stage_4: None, stage_5: None, stage_6: None,
            total_dv_ms: 0.0, total_tof_days: 0.0, mission_complete: false,
        }
    }

    // ========================================================================
    // STAGE 1: Injection from HEO to Unstable Manifold
    // ========================================================================
    fn execute_stage_1(&mut self) -> PyResult<StageResult> {
        if self.config.verbose {
            println!("🚀 STAGE 1: HEO → Unstable Manifold injection");
        }

        let xl1 = l1_position();
        let lambda = unstable_eigenvalue();
        let k = eigenvector_y_factor();

        let orbital = OrbitalState::new_canonical(
            xl1 - 0.003, -0.001, 0.0, 0.0, 0.0, 0.0,
        );

        if self.config.verbose {
            println!("   HEO apogee: x={:.6} (xL1={:.6})", orbital.position[0], xl1);
            println!("   λ={:.6}, k={:.6}", lambda, k);
        }

        let injection = compute_injection(
            &orbital,
            self.config.target_jacobi,
            0.005,
            1e-10,
            self.config.max_iterations,
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        if injection.state[0] <= xl1 {
            return Err(pyo3::exceptions::PyRuntimeError::new_err(
                format!("Stage 1 failed: x={:.6} ≤ xL1={:.6}", injection.state[0], xl1)
            ));
        }

        let canonical_state = vec![
            injection.state[0], injection.state[1], injection.state[2],
            injection.state[3], injection.state[4], injection.state[5],
        ];

        let burn_time_days = if injection.dv_total_ms > 0.0 {
            self.config.spacecraft_mass_kg * injection.dv_total_ms
                / self.config.engine_thrust_n / 86400.0
        } else {
            0.0
        };

        let v_canonical = (injection.state[3].powi(2) + injection.state[4].powi(2)).sqrt();
        let v_ms = v_canonical * V_CHAR;

        if self.config.verbose {
            println!("   ε (computed) = {:.6}", injection.epsilon);
            println!("   C_J achieved  = {:.8}", injection.jacobi_achieved);
            println!("   ΔV = {:.1} m/s", injection.dv_total_ms);
            println!("   |v| = {:.1} m/s (canonical: {:.4})", v_ms, v_canonical);
            println!("   x  = {:.6} > xL1 = {:.6} ✓ lunar realm", injection.state[0], xl1);
            println!("   Burn time = {:.2} days", burn_time_days);
        }

        let result = StageResult {
            stage_name: "Injection_HEO_to_Unstable_Manifold".to_string(),
            dv_total_ms: injection.dv_total_ms,
            time_of_flight_days: burn_time_days,
            final_state: injection.to_dimensional().to_vec(),
            canonical_state,
            jacobi_error: injection.jacobi_error,
            success: true,
            details: format!(
                "ε={:.6}, C_J={:.8}, ΔV={:.1}m/s, |v|={:.1}m/s, burn={:.2}d",
                injection.epsilon, injection.jacobi_achieved,
                injection.dv_total_ms, v_ms, burn_time_days,
            ),
        };

        self.total_dv_ms += injection.dv_total_ms;
        self.total_tof_days += burn_time_days;
        self.stage_1 = Some(result.clone());
        Ok(result)
    }

    // ========================================================================
    // STAGE 2: Ballistic Transit L1 → Moon
    // ========================================================================
    fn execute_stage_2(&mut self) -> PyResult<StageResult> {
        if self.config.verbose {
            println!("🌌 STAGE 2: Ballistic transit L1 → Moon");
        }

        let initial_state_canonical = match &self.stage_1 {
            Some(r) => {
                let s = &r.canonical_state;
                StateVector::new(s[0], s[1], s[2], s[3], s[4], s[5])
            }
            None => return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Must execute stage_1 first"
            )),
        };

        let xl1 = l1_position();
        if initial_state_canonical[0] <= xl1 {
            return Err(pyo3::exceptions::PyRuntimeError::new_err(
                format!("Stage 2 requires lunar realm: x={:.6} ≤ xL1={:.6}",
                        initial_state_canonical[0], xl1)
            ));
        }

        let config = IntegratorConfig {
            rtol: self.config.integration_rtol,
            verbose: self.config.verbose,
            max_step: 0.01,
            max_step_near_body: 0.0005,
            min_step: 1e-12,
            sampling_interval: 0.02,
            jacobi_monitor_freq: 500,
            ..IntegratorConfig::default()
        };

        let propagation = propagate_manifold::<EarthMoon>(
            initial_state_canonical.as_slice(),
            ManifoldType::TransitToMoon,
            self.config.target_jacobi,
            &config,
            60.0,
            Some(800.0),
            None,
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        let final_state_dim = denormalize_array(&propagation.final_state);
        let alt_km = lunar_altitude_km_from_dimensional(&final_state_dim);
        let tof_days = propagation.tof * T_CHAR / 86400.0;

        if self.config.verbose {
            println!("   TOF = {:.2} days", tof_days);
            println!("   Altitude = {:.0} km", alt_km);
        }

        let result = StageResult {
            stage_name: "Ballistic_Transit_to_Moon".to_string(),
            dv_total_ms: 0.0,
            time_of_flight_days: tof_days,
            final_state: final_state_dim,
            canonical_state: propagation.final_state.to_vec(),
            jacobi_error: propagation.jacobi_error,
            success: propagation.target_reached,
            details: format!(
                "alt={:.0}km, TOF={:.2}d, {}", alt_km, tof_days, propagation.termination_reason,
            ),
        };

        self.total_tof_days += tof_days;
        self.stage_2 = Some(result.clone());
        Ok(result)
    }

    // ========================================================================
    // STAGE 3: Lunar Capture Analysis
    // ========================================================================
    fn execute_stage_3(&mut self) -> PyResult<StageResult> {
        if self.config.verbose {
            println!("🌙 STAGE 3: Lunar capture analysis");
        }

        let final_state_dim = match &self.stage_2 {
            Some(r) => r.final_state.clone(),
            None => return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Must execute stage_2 first"
            )),
        };

        let capture: CaptureDetection = analyze_capture(
            &final_state_dim,
            self.config.target_jacobi,
            self.stage_2.as_ref().unwrap().time_of_flight_days,
            &CaptureConfig { verbose: self.config.verbose, ..CaptureConfig::default() },
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        let is_viable = capture.altitude_km > 0.0 && 
                        capture.dv_circularization_ms < 5000.0 &&
                        capture.eccentricity_estimate < 2.0;

        let result = StageResult {
            stage_name: "Lunar_Capture".to_string(),
            dv_total_ms: capture.dv_circularization_ms,
            time_of_flight_days: 0.0,
            final_state: capture.state.to_vec(),
            canonical_state: vec![0.0; 6],
            jacobi_error: (capture.jacobi_at_capture - self.config.target_jacobi).abs(),
            success: is_viable,
            details: format!(
                "alt={:.0}km, ΔV_circ={:.1}m/s, v={:.1}m/s, e={:.4}",
                capture.altitude_km, capture.dv_circularization_ms,
                capture.velocity_ms, capture.eccentricity_estimate,
            ),
        };

        self.total_dv_ms += capture.dv_circularization_ms;
        self.stage_3 = Some(result.clone());
        Ok(result)
    }

    // ========================================================================
    // STAGE 4: TEI desde LLO con SPICE
    // ========================================================================
    fn execute_stage_4(&mut self) -> PyResult<StageResult> {
        if self.config.verbose {
            println!("🚀 STAGE 4: TEI desde LLO 800 km (con SPICE)");
            println!("   Espera: {:.0} días", self.config.wait_days_lunar);
            println!("   SPICE: {}", if self.config.use_spice { "✅" } else { "❌" });
        }

        let llo_state = match &self.stage_3 {
            Some(r) => r.final_state.clone(),
            None => return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Must execute stage_3 first"
            )),
        };

        let return_config = ReturnConfig {
            target_jacobi: self.config.target_jacobi,
            wait_days: self.config.wait_days_lunar,
            target_perigee_km: 50000.0,
            llo_altitude_km: 800.0,
            verbose: self.config.verbose,
            kernel_dir: self.config.kernel_dir.clone(),
            reference_time_utc: self.config.reference_time_utc.clone(),
            use_spice: self.config.use_spice,
            ..ReturnConfig::default()
        };

        let tei: TeiResult = compute_tei_from_llo(
            &llo_state,
            &return_config,
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        let burn_time_days = if tei.dv_tei_ms > 0.0 {
            self.config.spacecraft_mass_kg * tei.dv_tei_ms
                / self.config.engine_thrust_n / 86400.0
        } else {
            0.0
        };

        if self.config.verbose {
            println!("   ΔV TEI = {:.1} m/s", tei.dv_tei_ms);
            println!("   C_J = {:.8}", tei.jacobi_achieved);
            println!("   TOF = {:.2} días", tei.tof_return_days);
            println!("   Burn = {:.2} days", burn_time_days);
            if tei.used_spice {
                println!("   🛰️  Usando SPICE para posición de la Tierra");
            }
        }

        let result = StageResult {
            stage_name: "Lunar_TEI".to_string(),
            dv_total_ms: tei.dv_tei_ms,
            time_of_flight_days: tei.tof_return_days + burn_time_days,
            final_state: tei.tei_state_dim.to_vec(),
            canonical_state: tei.tei_state.to_vec(),
            jacobi_error: tei.jacobi_error,
            success: tei.converged,
            details: format!(
                "ΔV={:.1}m/s, TOF={:.2}d, C_J={:.8}, perigeo={:.0}km{}",
                tei.dv_tei_ms, tei.tof_return_days, tei.jacobi_achieved,
                tei.expected_perigee_km,
                if tei.used_spice { " 🛰️" } else { " ⚠️" },
            ),
        };

        self.total_dv_ms += tei.dv_tei_ms;
        self.total_tof_days += tei.tof_return_days + burn_time_days;
        self.stage_4 = Some(result.clone());
        Ok(result)
    }

    // ========================================================================
    // STAGE 5: Propagación balística hasta perigeo HEO
    // ========================================================================
    fn execute_stage_5(&mut self) -> PyResult<StageResult> {
        if self.config.verbose {
            println!("🌌 STAGE 5: Tránsito Balístico → Perigeo HEO");
        }

        let tei_state = match &self.stage_4 {
            Some(r) => {
                let s = &r.canonical_state;
                [s[0], s[1], s[2], s[3], s[4], s[5]]
            }
            None => return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Must execute stage_4 first"
            )),
        };

        let tei_dv = self.stage_4.as_ref().unwrap().dv_total_ms;
        let return_config = ReturnConfig {
            target_jacobi: self.config.target_jacobi,
            verbose: self.config.verbose,
            ..ReturnConfig::default()
        };

        let prop: ReturnPropagationResult = propagate_ballistic_to_earth_perigee(
            &tei_state,
            tei_dv,
            20.0,
            &return_config,
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        let tof_days = prop.tof_return_days;

        if self.config.verbose {
            println!("   TOF = {:.2} days", tof_days);
            println!("   ΔV = 0.00 m/s (ballistic)");
            println!("   Perigeo HEO: {:.0} km", prop.perigee_altitude_km);
            println!("   Apogeo HEO: {:.0} km", prop.apogee_altitude_km);
            println!("   Excentricidad: {:.4}", prop.eccentricity);
            println!("   Captura: {}", if prop.captured { "✅" } else { "❌" });
        }

        let result = StageResult {
            stage_name: "Ballistic_Return_to_HEO".to_string(),
            dv_total_ms: 0.0,
            time_of_flight_days: tof_days,
            final_state: prop.perigee_state_dim.to_vec(),
            canonical_state: prop.perigee_state_canonical.to_vec(),
            jacobi_error: 0.0,
            success: prop.captured,
            details: format!(
                "TOF={:.2}d, perigee={:.0}km, apogee={:.0}km, e={:.4}, captured={}",
                tof_days, prop.perigee_altitude_km, prop.apogee_altitude_km,
                prop.eccentricity, prop.captured,
            ),
        };

        self.total_tof_days += tof_days;
        self.stage_5 = Some(result.clone());
        Ok(result)
    }

    // ========================================================================
    // STAGE 6: Circularización en perigeo HEO
    // ========================================================================
    fn execute_stage_6(&mut self) -> PyResult<StageResult> {
        if self.config.verbose {
            println!("🌍 STAGE 6: Circularización en Perigeo HEO");
        }

        let final_state_dim = match &self.stage_5 {
            Some(r) => r.final_state.clone(),
            None => return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Must execute stage_5 first"
            )),
        };

        let tof_days = self.stage_5.as_ref().unwrap().time_of_flight_days;
        let capture: EarthCapture = analyze_earth_capture(
            &final_state_dim,
            tof_days,
            &ReturnConfig::default(),
        );

        let burn_time_days = if capture.dv_circularization_ms > 0.0 {
            self.config.spacecraft_mass_kg * capture.dv_circularization_ms
                / self.config.engine_thrust_n / 86400.0
        } else {
            0.0
        };

        if self.config.verbose {
            println!("   Altitud perigeo: {:.0} km", capture.perigee_altitude_km);
            println!("   Velocidad: {:.1} m/s", capture.velocity_ms);
            println!("   Velocidad circular: {:.1} m/s", capture.circular_velocity_ms);
            println!("   ΔV circularización: {:.1} m/s", capture.dv_circularization_ms);
            println!("   Excentricidad: {:.4}", capture.eccentricity);
            println!("   Captura balística: {}", if capture.is_ballistic { "✅" } else { "❌" });
            println!("   Viable: {}", if capture.is_viable { "✅" } else { "❌" });
        }

        let result = StageResult {
            stage_name: "HEO_Circularization".to_string(),
            dv_total_ms: capture.dv_circularization_ms,
            time_of_flight_days: burn_time_days,
            final_state: capture.state.to_vec(),
            canonical_state: vec![0.0; 6],
            jacobi_error: 0.0,
            success: capture.is_viable,
            details: format!(
                "alt={:.0}km, v={:.1}m/s, v_circ={:.1}m/s, ΔV={:.1}m/s, e={:.4}",
                capture.perigee_altitude_km, capture.velocity_ms,
                capture.circular_velocity_ms, capture.dv_circularization_ms,
                capture.eccentricity,
            ),
        };

        self.total_dv_ms += capture.dv_circularization_ms;
        self.total_tof_days += burn_time_days;
        self.stage_6 = Some(result.clone());
        self.mission_complete = true;
        Ok(result)
    }

    // ========================================================================
    // PROPULSION BUDGET
    // ========================================================================
    fn compute_propulsion_budget(&self) -> PyResult<PyObject> {
        let dv1 = self.stage_1.as_ref().map(|s| s.dv_total_ms).unwrap_or(0.0);
        let dv3 = self.stage_3.as_ref().map(|s| s.dv_total_ms).unwrap_or(0.0);
        let dv4 = self.stage_4.as_ref().map(|s| s.dv_total_ms).unwrap_or(0.0);
        let dv6 = self.stage_6.as_ref().map(|s| s.dv_total_ms).unwrap_or(0.0);

        let tof_ida = self.stage_2.as_ref().map(|s| s.time_of_flight_days).unwrap_or(0.0);
        let tof_ret = self.stage_5.as_ref().map(|s| s.time_of_flight_days).unwrap_or(0.0);

        let budget = compute_mission_budget(
            dv1, dv3, dv4, dv6, tof_ida, tof_ret,
            &IonEngineConfig {
                initial_mass_kg: self.config.spacecraft_mass_kg,
                verbose: self.config.verbose,
                ..IonEngineConfig::default()
            },
        );

        Python::with_gil(|py| {
            let dict = PyDict::new(py);
            dict.set_item("initial_mass_kg", budget.initial_mass_kg)?;
            dict.set_item("final_mass_kg", budget.final_mass_kg)?;
            dict.set_item("total_propellant_kg", budget.total_propellant_kg)?;
            dict.set_item("total_mission_days", budget.total_mission_days)?;
            dict.set_item("mass_efficiency_pct", budget.mass_efficiency_pct)?;
            
            dict.set_item("dv_stage_1", dv1)?;
            dict.set_item("dv_stage_3", dv3)?;
            dict.set_item("dv_stage_4", dv4)?;
            dict.set_item("dv_stage_6", dv6)?;
            
            Ok(dict.into())
        })
    }

    // ========================================================================
    // FULL MISSION
    // ========================================================================
    fn execute_full_mission(&mut self) -> PyResult<PyObject> {
        if self.config.verbose {
            println!("╔══════════════════════════════════════════╗");
            println!("║   ASTRO-TFC: MISIÓN COMPLETA           ║");
            println!("║   Tierra → Luna → Tierra              ║");
            println!("║   C_J = {:.4}                         ║", self.config.target_jacobi);
            println!("║   Espera lunar: {:.0} días            ║", self.config.wait_days_lunar);
            println!("║   SPICE: {}                           ║", if self.config.use_spice { "✅" } else { "❌" });
            println!("╚══════════════════════════════════════════╝");
        }

        self.execute_stage_1()?;
        self.execute_stage_2()?;
        self.execute_stage_3()?;
        self.execute_stage_4()?;
        self.execute_stage_5()?;
        self.execute_stage_6()?;

        self.mission_complete = true;
        let budget = self.compute_propulsion_budget()?;

        Python::with_gil(|py| {
            let summary = PyDict::new(py);
            summary.set_item("total_dv_ms", self.total_dv_ms)?;
            summary.set_item("total_tof_days", self.total_tof_days)?;
            summary.set_item("mission_complete", self.mission_complete)?;
            summary.set_item("target_jacobi", self.config.target_jacobi)?;
            summary.set_item("propulsion_budget", budget)?;
            summary.set_item("used_spice", self.config.use_spice)?;
            summary.set_item("wait_days_lunar", self.config.wait_days_lunar)?;

            if let Some(ref s1) = self.stage_1 {
                summary.set_item("dv_stage_1", s1.dv_total_ms)?;
                summary.set_item("tof_stage_1", s1.time_of_flight_days)?;
            }
            if let Some(ref s2) = self.stage_2 {
                summary.set_item("tof_stage_2", s2.time_of_flight_days)?;
            }
            if let Some(ref s3) = self.stage_3 {
                summary.set_item("dv_stage_3", s3.dv_total_ms)?;
            }
            if let Some(ref s4) = self.stage_4 {
                summary.set_item("dv_stage_4", s4.dv_total_ms)?;
                summary.set_item("tof_stage_4", s4.time_of_flight_days)?;
            }
            if let Some(ref s5) = self.stage_5 {
                summary.set_item("tof_stage_5", s5.time_of_flight_days)?;
            }
            if let Some(ref s6) = self.stage_6 {
                summary.set_item("dv_stage_6", s6.dv_total_ms)?;
                summary.set_item("tof_stage_6", s6.time_of_flight_days)?;
            }

            if self.config.verbose {
                println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
                println!("🎉 MISIÓN COMPLETA EXITOSA");
                println!("   Total ΔV: {:.2} m/s", self.total_dv_ms);
                println!("   Total TOF: {:.1} días", self.total_tof_days);
                println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
            }

            Ok(summary.into())
        })
    }

    // ========================================================================
    // LOW THRUST OPTIMIZATION
    // ========================================================================
    fn optimize_low_thrust(
        &self,
        initial: Vec<f64>,
        target_altitude_km: f64,
        time_of_flight_days: f64,
    ) -> PyResult<(Vec<Vec<f64>>, Vec<Vec<f64>>)> {
        use low_thrust::{optimize_low_thrust_transfer, LowThrustConfig};

        if initial.len() != 6 {
            return Err(pyo3::exceptions::PyValueError::new_err("6 components required"));
        }

        // Usar la configuración de la misión
        let config = LowThrustConfig {
            thrust_max: self.config.engine_thrust_n,
            isp: self.config.engine_isp_s,
            initial_mass: self.config.spacecraft_mass_kg,
            ..LowThrustConfig::default()
        };

        if self.config.verbose {
            println!("🚀 Optimización de bajo empuje:");
            println!("   Estado inicial: {:?}", initial);
            println!("   Altitud objetivo: {:.1} km", target_altitude_km);
            println!("   TOF: {:.2} días", time_of_flight_days);
            println!("   Empuje máx: {:.3} N", self.config.engine_thrust_n);
            println!("   ISP: {:.0} s", self.config.engine_isp_s);
            println!("   Masa inicial: {:.1} kg", self.config.spacecraft_mass_kg);
        }

        let (trajectory, controls) = optimize_low_thrust_transfer(
            &[initial[0], initial[1], initial[2],
              initial[3], initial[4], initial[5]],
            target_altitude_km,
            time_of_flight_days,
            config,
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        if self.config.verbose {
            println!("   ✅ Optimización completada:");
            println!("   Puntos de trayectoria: {}", trajectory.len());
            println!("   Puntos de control: {}", controls.len());
            if !trajectory.is_empty() {
                println!("   Posición final: ({:.4}, {:.4})", 
                         trajectory.last().unwrap()[0], 
                         trajectory.last().unwrap()[1]);
            }
        }

        let traj_vec: Vec<Vec<f64>> = trajectory.iter().map(|s| s.to_vec()).collect();
        let ctrl_vec: Vec<Vec<f64>> = controls.iter()
            .map(|c| vec![c.throttle, c.alpha, c.beta]).collect();

        Ok((traj_vec, ctrl_vec))
    }
}

// ============================================================================
// STANDALONE PYTHON FUNCTIONS
// ============================================================================

#[pyfunction]
fn compute_jacobi_constant(state: Vec<f64>) -> PyResult<f64> {
    if state.len() != 6 {
        return Err(pyo3::exceptions::PyValueError::new_err("6 components required"));
    }

    let state_vec = if state[0].abs() > 10.0 {
        let d_char_km = D_CHAR / 1000.0;
        let v_char_kms = V_CHAR / 1000.0;
        StateVector::new(
            state[0] / d_char_km,
            state[1] / d_char_km,
            state[2] / d_char_km,
            state[3] / v_char_kms,
            state[4] / v_char_kms,
            state[5] / v_char_kms,
        )
    } else {
        StateVector::new(state[0], state[1], state[2], state[3], state[4], state[5])
    };

    Ok(jacobi_constant(&state_vec))
}

#[pyfunction]
fn compute_crtbp_derivatives(t: f64, state: Vec<f64>) -> PyResult<Vec<f64>> {
    if state.len() != 6 {
        return Err(pyo3::exceptions::PyValueError::new_err("6 components required"));
    }
    let s_can = normalize_state(&state);
    let d_can = crtbp_derivatives(t / T_CHAR, &s_can);
    Ok(vec![
        d_can[0] * V_CHAR, d_can[1] * V_CHAR, d_can[2] * V_CHAR,
        d_can[3] * D_CHAR / T_CHAR.powi(2),
        d_can[4] * D_CHAR / T_CHAR.powi(2),
        d_can[5] * D_CHAR / T_CHAR.powi(2),
    ])
}

#[pyfunction]
fn circular_velocity(r: f64, mu: Option<f64>) -> f64 {
    (mu.unwrap_or(MU_E) / r).sqrt()
}

#[pyfunction]
fn get_l1_position() -> f64 {
    l1_position()
}

#[pyfunction]
fn get_unstable_eigenvalue() -> f64 {
    unstable_eigenvalue()
}

#[pyfunction]
fn get_eigenvector_y_factor() -> f64 {
    eigenvector_y_factor()
}



// ============================================================================
// FUNCIÓN PARA CTYPES (compatible con scripts existentes)
// ============================================================================

#[no_mangle]
pub extern "C" fn validate_all_strategies() -> i32 {
    use ephemeris::validation::{
        MissionStrategyValidator,
        Smart1Benchmark,
        GrailBenchmark,
        CapstoneBenchmark,
    };

    println!("\n🚀 VALIDACIÓN DE ESTRATEGIAS DE MISIÓN");
    println!("   Comparando topologías de espacio de fases\n");
    
    let validator = match MissionStrategyValidator::new() {
        Ok(v) => v,
        Err(e) => {
            eprintln!("❌ Error inicializando validador: {}", e);
            return 1;
        }
    };
    
    // Validar SMART-1
    let smart1 = Smart1Benchmark;
    match validator.validate(&smart1) {
        Ok(_) => println!("✅ SMART-1 validada\n"),
        Err(e) => eprintln!("❌ Error SMART-1: {}\n", e),
    }
    
    // Validar GRAIL
    let grail = GrailBenchmark;
    match validator.validate(&grail) {
        Ok(_) => println!("✅ GRAIL validada\n"),
        Err(e) => eprintln!("❌ Error GRAIL: {}\n", e),
    }
    
    // Validar CAPSTONE
    let capstone = CapstoneBenchmark;
    match validator.validate(&capstone) {
        Ok(_) => println!("✅ CAPSTONE validada\n"),
        Err(e) => eprintln!("❌ Error CAPSTONE: {}\n", e),
    }
    
    println!("🎉 VALIDACIÓN COMPLETADA");
    0
}

// ============================================================================
// PYTHON MODULE
// ============================================================================

#[pymodule]
fn astro_tfc(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add("MU_E", MU_E)?;
    m.add("MU_M", MU_M)?;
    m.add("R_EARTH", R_EARTH)?;
    m.add("R_MOON", R_MOON)?;
    m.add("L_DISTANCE", L_DISTANCE)?;
    m.add("OMEGA", OMEGA)?;
    m.add("D1", D1)?;
    m.add("D2", D2)?;
    m.add("T_CHAR", T_CHAR)?;
    m.add("V_CHAR", V_CHAR)?;
    m.add("D_CHAR", D_CHAR)?;
    m.add("MU_NORMALIZED", MU_NORMALIZED)?;

    m.add_class::<MissionConfig>()?;
    m.add_class::<StageResult>()?;
    m.add_class::<AstroTFCMission>()?;

    m.add_function(wrap_pyfunction!(compute_jacobi_constant, m)?)?;
    m.add_function(wrap_pyfunction!(compute_crtbp_derivatives, m)?)?;
    m.add_function(wrap_pyfunction!(circular_velocity, m)?)?;
    m.add_function(wrap_pyfunction!(get_l1_position, m)?)?;
    m.add_function(wrap_pyfunction!(get_unstable_eigenvalue, m)?)?;
    m.add_function(wrap_pyfunction!(get_eigenvector_y_factor, m)?)?;

    m.add("__version__", "3.5.0")?;
    m.add("__author__", "Omar Ariel Vallejos")?;
    m.add("__reference__", "Almeida Jr. et al. (2026) - Astrodynamics")?;
    m.add("__validation__", "Vallejos (2026) DOI: 10.5281/zenodo.20584996")?;

    Ok(())
}
