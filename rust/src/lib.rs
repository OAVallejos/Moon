//! AstroTFC: Simulation engine for Earth-Moon transfers via L1
//! using the Theory of Functional Connections (TFC).
//! 
//! Based on: Almeida Jr. et al. (2026) - Astrodynamics
//! "Earth–Moon transfer via the L1 Lagrangian point using 
//!  the theory of functional connections"
//! https://doi.org/10.1007/s42064-025-0297-x
//!
//! COMPLETE ARCHITECTURE (v2.4.1 - June 2026):
//! - Stage 1: Injection to L1 STABLE manifold (geocentric.rs)
//! - Stage 2: Ballistic transit Earth→Moon (integration.rs) [CORRECTED]
//! - Stage 3: Ballistic lunar capture (lunar.rs)
//! - Stage 4: Lunar escape → L1 UNSTABLE manifold (return.rs)
//! - Stage 5: Ballistic transit Moon→Earth (integration.rs)
//! - Stage 6: Earth capture (return.rs)
//! - Budget: ion propulsion (propulsion.rs)
//!
//! CORRECTIONS v2.4.1:
//! - Stage 2: ManifoldType::UNSTABLE (direction +1, toward the Moon)
//! - Stage 2: Initial state in CANONICAL units
//! - Stage 5: ManifoldType::STABLE (direction -1, toward Earth)
//! - Consistent canonical↔dimensional conversion

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use pyo3::types::PyDict;

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// MÓDULOS INTERNOS
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
mod constants;
mod crtbp;
mod tfc;
mod geocentric;
mod lunar;
mod manifold;
mod integration;
mod r#return;
mod propulsion;

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// IMPORTACIONES
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
use constants::{
    MU_E, MU_M, L_DISTANCE, R_EARTH, R_MOON, OMEGA, D1, D2,
    T_CHAR, D_CHAR, V_CHAR, MU_NORMALIZED,
};
use crtbp::{StateVector, crtbp_derivatives, jacobi_constant};
use geocentric::{
    OrbitalState, compute_injection,
};
use lunar::{
    CaptureConfig, analyze_capture,
};
use manifold::{
    l1_position, unstable_eigenvalue, eigenvector_y_factor,
    ManifoldType,
};
use integration::{
    EarthMoon, IntegratorConfig,
    propagate_manifold, lunar_altitude_km,
};
use r#return::{
    ReturnConfig,
    compute_escape_to_unstable_manifold, analyze_earth_capture,
};
use propulsion::{
    IonEngineConfig,
    compute_mission_budget,
};

// ═══════════════════════════════════════════════════════════════
// CAPA DE NORMALIZACIÓN (UNIDADES CANÓNICAS ↔ DIMENSIONALES)
// ═══════════════════════════════════════════════════════════════

/// Convierte estado dimensional [m, m, m, m/s, m/s, m/s] a canónico.
fn normalize_state(dimensional: &[f64]) -> StateVector {
    StateVector::new(
        dimensional[0] / D_CHAR,
        dimensional[1] / D_CHAR,
        dimensional[2] / D_CHAR,
        dimensional[3] / V_CHAR,
        dimensional[4] / V_CHAR,
        dimensional[5] / V_CHAR,
    )
}

/// Convierte estado canónico a dimensional [m, m, m, m/s, m/s, m/s].
fn denormalize_state(canonical: &StateVector) -> Vec<f64> {
    vec![
        canonical.x * D_CHAR,
        canonical.y * D_CHAR,
        canonical.z * D_CHAR,
        canonical[3] * V_CHAR,
        canonical[4] * V_CHAR,
        canonical[5] * V_CHAR,
    ]
}

/// Convierte array canónico [f64; 6] a vector dimensional.
fn denormalize_array(canonical: &[f64; 6]) -> Vec<f64> {
    vec![
        canonical[0] * D_CHAR,
        canonical[1] * D_CHAR,
        canonical[2] * D_CHAR,
        canonical[3] * V_CHAR,
        canonical[4] * V_CHAR,
        canonical[5] * V_CHAR,
    ]
}

fn _normalize_radius(r_meters: f64) -> f64 { r_meters / D_CHAR }
fn _normalize_velocity(v_ms: f64) -> f64 { v_ms / V_CHAR }
fn _denormalize_dv(dv_canonical: f64) -> f64 { dv_canonical * V_CHAR }
fn _denormalize_time(t_canonical: f64) -> f64 { t_canonical * T_CHAR }

// ═══════════════════════════════════════════════════════════════
// CLASES PYTHON
// ═══════════════════════════════════════════════════════════════

#[pyclass]
#[derive(Clone)]
struct MissionConfig {
    #[pyo3(get, set)]
    leo_altitude_m: f64,
    #[pyo3(get, set)]
    target_jacobi: f64,
    #[pyo3(get, set)]
    n_chebyshev_points: usize,
    #[pyo3(get, set)]
    integration_rtol: f64,
    #[pyo3(get, set)]
    max_iterations: usize,
    #[pyo3(get, set)]
    verbose: bool,
}

#[pymethods]
impl MissionConfig {
    #[new]
    fn new() -> Self {
        MissionConfig {
            leo_altitude_m: 200_000.0,
            target_jacobi: 3.187,  // CORREGIDO: valor canónico por defecto
            n_chebyshev_points: 50,
            integration_rtol: 1e-12,
            max_iterations: 500,
            verbose: false,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "MissionConfig(LEO={:.0}km, C_J={:.4}, n_cheb={})",
            self.leo_altitude_m / 1000.0,
            self.target_jacobi,
            self.n_chebyshev_points,
        )
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
    jacobi_error: f64,
    #[pyo3(get)]
    success: bool,
    #[pyo3(get)]
    details: String,
}

#[pymethods]
impl StageResult {
    fn __repr__(&self) -> String {
        format!(
            "StageResult({}, ΔV={:.1}m/s, TOF={:.1}d, success={})",
            self.stage_name, self.dv_total_ms, self.time_of_flight_days, self.success
        )
    }

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
            stage_1: None,
            stage_2: None,
            stage_3: None,
            stage_4: None,
            stage_5: None,
            stage_6: None,
            total_dv_ms: 0.0,
            total_tof_days: 0.0,
            mission_complete: false,
        }
    }

    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    // ETAPA 1: Inyección a variedad ESTABLE L1
    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    fn execute_stage_1(&mut self) -> PyResult<StageResult> {
        if self.config.verbose {
            println!("🚀 ETAPA 1: Inyección a variedad ESTABLE L1...");
        }

        let xl1 = l1_position();

        // Estado orbital post-lanzamiento en cuenca terrestre (UNIDADES CANÓNICAS)
        let orbital = OrbitalState::new_canonical(
            xl1 - 0.01, -0.001, 0.0,
            0.0, 0.0, 0.0,
        );

        let injection = compute_injection(
            &orbital,
            self.config.target_jacobi,
            0.005,
            1e-10,
            30,
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        // Convertir a dimensional para Python
        let final_state_dim = injection.to_dimensional().to_vec();

        let stage_result = StageResult {
            stage_name: "Injection_to_Stable_Manifold".to_string(),
            dv_total_ms: injection.dv_total_ms,
            time_of_flight_days: 0.0,
            final_state: final_state_dim,
            jacobi_error: injection.jacobi_error,
            success: true,
            details: format!(
                "ε={:.6}, C_J={:.8}, ΔV_radial={:.1}m/s, ΔV_tang={:.1}m/s",
                injection.epsilon, injection.jacobi_achieved,
                injection.dv_radial_ms, injection.dv_tangential_ms,
            ),
        };

        if self.config.verbose {
            println!("   ✅ ΔV = {:.2} m/s", injection.dv_total_ms);
            println!("   ✅ C_J = {:.8}", injection.jacobi_achieved);
        }

        self.stage_1 = Some(stage_result.clone());
        self.total_dv_ms += injection.dv_total_ms;
        Ok(stage_result)
    }

    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    // ETAPA 2: Tránsito balístico Tierra → Luna
    // CORREGIDO: Usar estado de Etapa 1, propagar Unstable
    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    fn execute_stage_2(&mut self) -> PyResult<StageResult> {
        if self.config.verbose {
            println!("🌌 ETAPA 2: Tránsito balístico L1 → Luna (variedad INESTABLE)...");
        }

        let initial_state_dim = match &self.stage_1 {
            Some(r) => r.final_state.clone(),
            None => return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "¡Ejecuta execute_stage_1 primero!"
            )),
        };

        // Convertir a canónico
        let initial_state_canonical = normalize_state(&initial_state_dim);

        if self.config.verbose {
            let cj_check = jacobi_constant(&initial_state_canonical);
            println!("   C_J inicial (canónico): {:.10}", cj_check);
            println!("   Estado: x={:.6}, y={:.6}, vx={:.6}, vy={:.6}",
                initial_state_canonical[0], initial_state_canonical[1],
                initial_state_canonical[3], initial_state_canonical[4]);
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

        // Propagar en dirección Unstable (hacia adelante en el tiempo)
        // El estado ya está en la cuenca terrestre (x < xL1) con vx > 0
                let propagation = propagate_manifold::<EarthMoon>(
            initial_state_canonical.as_slice(),
            ManifoldType::TransitToMoon,
            self.config.target_jacobi,
            &config,
            60.0,
            Some(500.0),
            None,  // Sin umbral de escape L1
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        let final_state_dim = denormalize_array(&propagation.final_state);
        let alt_km = lunar_altitude_km(&final_state_dim);

        let stage_result = StageResult {
            stage_name: "Ballistic_Transit_to_Moon".to_string(),
            dv_total_ms: 0.0,
            time_of_flight_days: propagation.tof * T_CHAR / 86400.0,
            final_state: final_state_dim,
            jacobi_error: propagation.jacobi_error,
            success: propagation.target_reached,
            details: format!(
                "Altitud lunar={:.0}km, pasos={}/{}, deriva_max_CJ={:.2e}, {}",
                alt_km, propagation.steps_accepted, propagation.steps_attempted,
                propagation.jacobi_error_max, propagation.termination_reason,
            ),
        };

        if self.config.verbose {
            println!("   ✅ ΔV = 0.00 m/s (balístico)");
            println!("   ✅ TOF = {:.1} días", stage_result.time_of_flight_days);
            println!("   ✅ Altitud = {:.0} km", alt_km);
            println!("   ✅ {}", propagation.termination_reason);
        }

        self.stage_2 = Some(stage_result.clone());
        self.total_tof_days += stage_result.time_of_flight_days;
        Ok(stage_result)
    }

    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    // ETAPA 3: Análisis de captura lunar
    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    fn execute_stage_3(&mut self) -> PyResult<StageResult> {
        if self.config.verbose {
            println!("🌙 ETAPA 3: Análisis de captura balística lunar...");
        }

        let final_state_dim = match &self.stage_2 {
            Some(r) => r.final_state.clone(),
            None => return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "¡Ejecuta execute_stage_2 primero!"
            )),
        };

        let capture_config = CaptureConfig {
            verbose: self.config.verbose,
            ..CaptureConfig::default()
        };

        let capture = analyze_capture(
            &final_state_dim,
            self.config.target_jacobi,
            self.stage_2.as_ref().unwrap().time_of_flight_days,
            &capture_config,
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        let stage_result = StageResult {
            stage_name: "Lunar_Capture_Analysis".to_string(),
            dv_total_ms: capture.dv_circularization_ms,
            time_of_flight_days: 0.0,
            final_state: capture.state.to_vec(),
            jacobi_error: (capture.jacobi_at_capture - self.config.target_jacobi).abs(),
            success: capture.is_viable,
            details: format!(
                "Altitud={:.0}km, v={:.1}m/s, v_circ={:.1}m/s, e={:.4}, {}",
                capture.altitude_km, capture.velocity_ms,
                capture.circular_velocity_ms, capture.eccentricity_estimate,
                capture.diagnosis,
            ),
        };

        if self.config.verbose {
            println!("   ✅ Altitud captura: {:.0} km", capture.altitude_km);
            println!("   ✅ ΔV circularización: {:.1} m/s", capture.dv_circularization_ms);
            println!("   ✅ Viable: {}", capture.is_viable);
        }

        self.stage_3 = Some(stage_result.clone());
        self.total_dv_ms += capture.dv_circularization_ms;
        Ok(stage_result)
    }

    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    // ETAPA 4: Escape lunar → variedad INESTABLE L1
    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    fn execute_stage_4(&mut self) -> PyResult<StageResult> {
        if self.config.verbose {
            println!("🚀 ETAPA 4: Escape lunar → variedad INESTABLE L1...");
        }

        let llo_state = match &self.stage_3 {
            Some(r) => r.final_state.clone(),
            None => return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "¡Ejecuta execute_stage_3 primero!"
            )),
        };

        let return_config = ReturnConfig {
            target_jacobi: self.config.target_jacobi,
            verbose: self.config.verbose,
            ..ReturnConfig::default()
        };

        let escape = compute_escape_to_unstable_manifold(&llo_state, &return_config)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        let stage_result = StageResult {
            stage_name: "Lunar_Escape_to_Unstable_Manifold".to_string(),
            dv_total_ms: escape.dv_escape_ms,
            time_of_flight_days: 0.0,
            final_state: escape.state.to_vec(),
            jacobi_error: escape.jacobi_error,
            success: true,
            details: format!(
                "ε={:.6}, C_J={:.8}",
                escape.epsilon, escape.jacobi_achieved,
            ),
        };

        if self.config.verbose {
            println!("   ✅ ΔV escape = {:.2} m/s", escape.dv_escape_ms);
            println!("   ✅ C_J = {:.8}", escape.jacobi_achieved);
        }

        self.stage_4 = Some(stage_result.clone());
        self.total_dv_ms += escape.dv_escape_ms;
        Ok(stage_result)
    }

    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    // ETAPA 5: Tránsito balístico Luna → Tierra
    // CORREGIDO: Estado CANÓNICO
    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    fn execute_stage_5(&mut self) -> PyResult<StageResult> {
        if self.config.verbose {
            println!("🌌 ETAPA 5: Tránsito balístico Luna→Tierra (variedad ESTABLE)...");
        }

        let initial_state_dim = match &self.stage_4 {
            Some(r) => r.final_state.clone(),
            None => return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "¡Ejecuta execute_stage_4 primero!"
            )),
        };

        // Convertir a canónico
        let initial_state_canonical = normalize_state(&initial_state_dim);

        let config = IntegratorConfig {
            rtol: self.config.integration_rtol,
            verbose: self.config.verbose,
            max_step: 0.01,
            sampling_interval: 0.02,
            ..IntegratorConfig::default()
        };

        // Para retorno a Tierra: Stable (dirección -1)
                let propagation = propagate_manifold::<EarthMoon>(
            initial_state_canonical.as_slice(),
            ManifoldType::TransitToEarth,        // CORRECTO: Stable para ir hacia la Tierra
            self.config.target_jacobi,
            &config,
            40.0,
            None,
            Some(0.05),
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        let final_state_dim = denormalize_array(&propagation.final_state);

        let stage_result = StageResult {
            stage_name: "Ballistic_Return_to_Earth".to_string(),
            dv_total_ms: 0.0,
            time_of_flight_days: propagation.tof * T_CHAR / 86400.0,
            final_state: final_state_dim,
            jacobi_error: propagation.jacobi_error,
            success: propagation.target_reached,
            details: format!(
                "pasos={}/{}, deriva_max_CJ={:.2e}, {}",
                propagation.steps_accepted, propagation.steps_attempted,
                propagation.jacobi_error_max, propagation.termination_reason,
            ),
        };

        if self.config.verbose {
            println!("   ✅ ΔV = 0.00 m/s (balístico)");
            println!("   ✅ TOF = {:.1} días", stage_result.time_of_flight_days);
        }

        self.stage_5 = Some(stage_result.clone());
        self.total_tof_days += stage_result.time_of_flight_days;
        Ok(stage_result)
    }

    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    // ETAPA 6: Captura terrestre
    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    fn execute_stage_6(&mut self) -> PyResult<StageResult> {
        if self.config.verbose {
            println!("🌍 ETAPA 6: Captura terrestre...");
        }

        let final_state_dim = match &self.stage_5 {
            Some(r) => r.final_state.clone(),
            None => return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "¡Ejecuta execute_stage_5 primero!"
            )),
        };

        let return_config = ReturnConfig::default();
        let tof_days = self.stage_5.as_ref().unwrap().time_of_flight_days;

        let capture = analyze_earth_capture(&final_state_dim, tof_days, &return_config);

        let stage_result = StageResult {
            stage_name: "Earth_Capture".to_string(),
            dv_total_ms: capture.dv_circularization_ms,
            time_of_flight_days: 0.0,
            final_state: capture.state.to_vec(),
            jacobi_error: 0.0,
            success: capture.is_ballistic,
            details: format!(
                "Altitud={:.0}km, v={:.1}m/s, v_circ={:.1}m/s, balístico={}",
                capture.perigee_altitude_km, capture.velocity_ms,
                capture.circular_velocity_ms, capture.is_ballistic,
            ),
        };

        if self.config.verbose {
            println!("   ✅ Altitud captura: {:.0} km", capture.perigee_altitude_km);
            println!("   ✅ ΔV circularización: {:.3} m/s", capture.dv_circularization_ms);
            println!("   ✅ Captura balística: {}", if capture.is_ballistic { "✅ SÍ" } else { "❌ NO" });
        }

        self.stage_6 = Some(stage_result.clone());
        self.total_dv_ms += capture.dv_circularization_ms;
        Ok(stage_result)
    }

    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    // PRESUPUESTO DE PROPULSIÓN
    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    fn compute_propulsion_budget(&self) -> PyResult<PyObject> {
        let dv1 = self.stage_1.as_ref().map(|s| s.dv_total_ms).unwrap_or(0.0);
        let dv3 = self.stage_3.as_ref().map(|s| s.dv_total_ms).unwrap_or(0.0);
        let dv4 = self.stage_4.as_ref().map(|s| s.dv_total_ms).unwrap_or(0.0);
        let dv6 = self.stage_6.as_ref().map(|s| s.dv_total_ms).unwrap_or(0.0);

        let tof_ida = self.stage_2.as_ref().map(|s| s.time_of_flight_days).unwrap_or(0.0);
        let tof_ret = self.stage_5.as_ref().map(|s| s.time_of_flight_days).unwrap_or(0.0);

        let engine_config = IonEngineConfig {
            initial_mass_kg: 250.0,
            verbose: self.config.verbose,
            ..IonEngineConfig::default()
        };

        let budget = compute_mission_budget(
            dv1, dv3, dv4, dv6, tof_ida, tof_ret, &engine_config,
        );

        Python::with_gil(|py| {
            let dict = PyDict::new(py);
            dict.set_item("initial_mass_kg", budget.initial_mass_kg)?;
            dict.set_item("final_mass_kg", budget.final_mass_kg)?;
            dict.set_item("total_propellant_kg", budget.total_propellant_kg)?;
            dict.set_item("total_mission_days", budget.total_mission_days)?;
            dict.set_item("mass_efficiency_pct", budget.mass_efficiency_pct)?;
            dict.set_item("stage1_propellant_kg", budget.stage1.propellant_with_margins_kg)?;
            dict.set_item("stage3_propellant_kg", budget.stage3.propellant_with_margins_kg)?;
            dict.set_item("stage4_propellant_kg", budget.stage4.propellant_with_margins_kg)?;
            dict.set_item("stage6_propellant_kg", budget.stage6.propellant_with_margins_kg)?;
            Ok(dict.into())
        })
    }

    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    // MISIÓN COMPLETA
    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    fn execute_full_mission(&mut self) -> PyResult<PyObject> {
        if self.config.verbose {
            println!("╔══════════════════════════════════════════╗");
            println!("║   ASTRO-TFC: MISIÓN COMPLETA            ║");
            println!("║   Tierra → Luna → Tierra                ║");
            println!("╚══════════════════════════════════════════╝");
            println!("   C_J objetivo: {:.4}", self.config.target_jacobi);
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

            if let Some(ref s1) = self.stage_1 {
                summary.set_item("dv_stage_1", s1.dv_total_ms)?;
            }
            if let Some(ref s3) = self.stage_3 {
                summary.set_item("dv_stage_3", s3.dv_total_ms)?;
            }
            if let Some(ref s4) = self.stage_4 {
                summary.set_item("dv_stage_4", s4.dv_total_ms)?;
            }
            if let Some(ref s6) = self.stage_6 {
                summary.set_item("dv_stage_6", s6.dv_total_ms)?;
            }

            if self.config.verbose {
                println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
                println!("🎉 MISIÓN COMPLETADA");
                println!("   ΔV Total: {:.2} m/s", self.total_dv_ms);
                println!("   Tiempo total: {:.1} días", self.total_tof_days);
            }

            Ok(summary.into())
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "AstroTFCMission(ΔV={:.1}m/s, TOF={:.1}d, complete={})",
            self.total_dv_ms, self.total_tof_days, self.mission_complete
        )
    }
}

// ═══════════════════════════════════════════════════════════════
// FUNCIONES PYTHON INDEPENDIENTES
// ═══════════════════════════════════════════════════════════════

#[pyfunction]
fn compute_jacobi_constant(state: Vec<f64>) -> PyResult<f64> {
    if state.len() != 6 {
        return Err(pyo3::exceptions::PyValueError::new_err("6 componentes requeridas"));
    }
    let s_can = normalize_state(&state);
    Ok(jacobi_constant(&s_can))
}

#[pyfunction]
fn compute_crtbp_derivatives(t: f64, state: Vec<f64>) -> PyResult<Vec<f64>> {
    if state.len() != 6 {
        return Err(pyo3::exceptions::PyValueError::new_err("6 componentes requeridas"));
    }
    let s_can = normalize_state(&state);
    let t_can = t / T_CHAR;
    let d_can = crtbp_derivatives(t_can, &s_can);
    Ok(vec![
        d_can.x * V_CHAR,
        d_can.y * V_CHAR,
        d_can.z * V_CHAR,
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

// ═══════════════════════════════════════════════════════════════
// MÓDULO PYTHON
// ═══════════════════════════════════════════════════════════════

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

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__author__", "Omar Ariel Vallejos")?;
    m.add("__reference__", "Almeida Jr. et al. (2026) - Astrodynamics")?;

    Ok(())
}
