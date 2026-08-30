//! ephemeris/planner.rs - Planificador de ventanas de lanzamiento

use super::*;
use crate::constants::D_CHAR;

/// Ventana de lanzamiento encontrada
#[derive(Debug, Clone)]
pub struct LaunchWindow {
    pub time_et: f64,
    pub time_utc: String,
    pub distance_km: f64,
    pub relative_velocity_kms: f64,
    pub v_inf_kms: f64,
    pub phase_angle_deg: f64,
    pub merit_score: f64,
    pub synodic_period_days: f64,
}

/// Planificador de ventanas de lanzamiento
pub struct LaunchWindowPlanner {
    spice: SpiceContext,
    config: EphemerisConfig,
}

impl LaunchWindowPlanner {
    pub fn new(spice: SpiceContext, config: EphemerisConfig) -> Self {
        LaunchWindowPlanner { spice, config }
    }

    /// Encuentra ventanas de lanzamiento en un rango de tiempo
    pub fn find_windows(
        &self,
        target_body: i32,
        start_utc: &str,
        end_utc: &str,
        step_days: f64,
    ) -> Result<Vec<LaunchWindow>, SpiceError> {
        let start_et = self.spice.utc_to_et(start_utc)?;
        let end_et = self.spice.utc_to_et(end_utc)?;
        
        let mut windows = Vec::new();
        let mut et = start_et;
        let step_seconds = step_days * 86400.0;
        
        while et < end_et {
            // Obtener estados
            let target = self.spice.get_state(target_body, et, 0)?;
            let earth = self.spice.get_state(399, et, 0)?; // Earth = 399
            
            // Calcular geometría
            let distance = target.distance_to(&earth);
            let rel_vel = target.relative_velocity(&earth);
            
            // V_inf estimada
            let v_inf = self.estimate_v_inf(&target, &earth, et)?;
            
            // Ángulo de fase
            let phase = self.calculate_phase(&target, &earth, et)?;
            
            // Mérito
            let merit = (distance / 1.5e8) + (rel_vel / 20.0);
            
            // Período sinódico
            let synodic = self.calculate_synodic_period(target_body)?;
            
            windows.push(LaunchWindow {
                time_et: et,
                time_utc: self.spice.et_to_utc(et)?,
                distance_km: distance,
                relative_velocity_kms: rel_vel,
                v_inf_kms: v_inf,
                phase_angle_deg: phase,
                merit_score: merit,
                synodic_period_days: synodic,
            });
            
            et += step_seconds;
        }
        
        // Ordenar por mérito (mejor primero)
        windows.sort_by(|a, b| a.merit_score.partial_cmp(&b.merit_score).unwrap());
        
        Ok(windows)
    }

    fn estimate_v_inf(
        &self,
        target: &SpiceState,
        earth: &SpiceState,
        et: f64,
    ) -> Result<f64, SpiceError> {
        let rel_vel = target.relative_velocity(earth);
        let mu = self.spice.get_mu(399)?; // Earth GM
        
        // Velocidad de escape de la Tierra en la posición actual
        let r_earth = earth.distance() * 1000.0; // metros
        let v_esc = (2.0 * mu / r_earth).sqrt();
        
        // V_inf² = V_rel² - V_esc²
        let v_inf_sq = (rel_vel * 1000.0).powi(2) - v_esc.powi(2);
        
        if v_inf_sq <= 0.0 {
            Ok(0.0)
        } else {
            Ok(v_inf_sq.sqrt() / 1000.0)
        }
    }

    fn calculate_phase(
        &self,
        target: &SpiceState,
        earth: &SpiceState,
        et: f64,
    ) -> Result<f64, SpiceError> {
        let sun = self.spice.get_state(10, et, 0)?; // Sun = 10
        
        let v1 = [
            target.position[0] - sun.position[0],
            target.position[1] - sun.position[1],
            target.position[2] - sun.position[2],
        ];
        let v2 = [
            earth.position[0] - sun.position[0],
            earth.position[1] - sun.position[1],
            earth.position[2] - sun.position[2],
        ];
        
        let dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
        let n1 = (v1[0].powi(2) + v1[1].powi(2) + v1[2].powi(2)).sqrt();
        let n2 = (v2[0].powi(2) + v2[1].powi(2) + v2[2].powi(2)).sqrt();
        
        let cos_phase = dot / (n1 * n2);
        let phase = cos_phase.acos() * 180.0 / std::f64::consts::PI;
        
        Ok(phase)
    }

    fn calculate_synodic_period(&self, target_body: i32) -> Result<f64, SpiceError> {
    let mu = self.spice.get_mu(10)?; // Sun GM
    
    // Períodos orbitales (ley de Kepler)
    let a_earth: f64 = 1.495978707e11;
    let a_target: f64 = match target_body {
        499 => 2.279e11,  // Mars
        299 => 1.082e11,  // Venus
        599 => 7.785e11,  // Jupiter
        _ => return Ok(365.25), // Default
    };
    
    let t_earth = 2.0 * std::f64::consts::PI * (a_earth.powi(3) / mu).sqrt();
    let t_target = 2.0 * std::f64::consts::PI * (a_target.powi(3) / mu).sqrt();
    
    // Período sinódico: 1/T_syn = |1/T_earth - 1/T_target|
    let t_syn = 1.0 / (1.0/t_earth - 1.0/t_target).abs();
    
    Ok(t_syn / 86400.0) // Convertir a días
 }
}
