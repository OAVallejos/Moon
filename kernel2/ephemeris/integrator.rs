//! ephemeris/integrator.rs - Integrador de n-cuerpos con SPICE
//! Versión mejorada que usa el SpiceContext unificado

use super::spice::{SpiceContext, SpiceState, SpiceError};
use crate::integration::IntegratorConfig;
use crate::constants::D_CHAR;

pub type NBodyState = [f64; 7]; // [x, y, z, vx, vy, vz, time_et]

/// Integrador de n-cuerpos con efemérides SPICE
pub struct NBodyIntegrator {
    spice: SpiceContext,
    config: IntegratorConfig,
    bodies: Vec<String>,
    /// Factor de escala para convertir km ↔ m (SPICE devuelve km)
    km_to_m: f64,
}

impl NBodyIntegrator {
    pub fn new(spice: SpiceContext, config: IntegratorConfig, bodies: Vec<String>) -> Self {
        NBodyIntegrator {
            spice,
            config,
            bodies,
            km_to_m: 1000.0,
        }
    }

    /// Propaga el estado usando efemérides SPICE
    /// 
    /// # Argumentos
    /// * `initial_state` - Estado inicial [x, y, z, vx, vy, vz] en **metros** y **m/s**
    /// * `initial_et` - Tiempo inicial en **Ephemeris Time (ET)** [segundos desde J2000]
    /// * `duration_seconds` - Duración de la propagación [segundos]
    /// * `max_steps` - Número máximo de pasos
    /// 
    /// # Retorna
    /// * `Vec<NBodyState>` - Trayectoria con estados [x, y, z, vx, vy, vz, time_et]
    pub fn propagate(
        &mut self,
        initial_state: &[f64; 6],
        initial_et: f64,
        duration_seconds: f64,
        max_steps: usize,
    ) -> Result<Vec<NBodyState>, SpiceError> {
        // Verificar que SPICE está cargado
        if !self.spice.is_loaded() {
            return Err(SpiceError::NotLoaded);
        }

        let mut trajectory = Vec::with_capacity(max_steps);
        
        // Estado: [x(m), y(m), z(m), vx(m/s), vy(m/s), vz(m/s), time_et(s)]
        let mut state = [
            initial_state[0], initial_state[1], initial_state[2],
            initial_state[3], initial_state[4], initial_state[5],
            initial_et,
        ];
        
        let mut t = initial_et;
        // max_step en días → segundos
        let dt = self.config.max_step * 86400.0;
        let mut steps = 0;

        // Validar paso
        if dt <= 0.0 || dt > 86400.0 {
            // Si el paso es inválido, usar un valor por defecto (1 hora)
            let dt_default = 3600.0;
            println!("⚠️ Integrator: max_step inválido ({}), usando {:.0}s", dt, dt_default);
            let dt = dt_default;
        }

        trajectory.push(state);

        while steps < max_steps && (t - initial_et) < duration_seconds {
            // Derivadas con SPICE (RK4)
            let k1 = self.derivatives(&state)?;
            let k2 = self.derivatives(&self.step_state(&state, &k1, dt/2.0)?)?;
            let k3 = self.derivatives(&self.step_state(&state, &k2, dt/2.0)?)?;
            let k4 = self.derivatives(&self.step_state(&state, &k3, dt)?)?;

            // RK4
            for i in 0..6 {
                state[i] += dt * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
            }
            state[6] += dt;
            t += dt;
            steps += 1;

            // Verificar límite del sistema (escape)
            // Si la nave se aleja más de 2.0e9 m (~5.2 veces la distancia Tierra-Luna)
            // consideramos que ha escapado
            let r = (state[0].powi(2) + state[1].powi(2) + state[2].powi(2)).sqrt();
            let escape_limit = 2.0e9; // 2,000,000,000 m
            if r > escape_limit {
                // Añadir estado actual y terminar
                trajectory.push(state);
                break;
            }

            // Verificar impacto con cuerpos (cada 10 pasos para eficiencia)
            if steps % 10 == 0 {
                let mut impact_detected = false;
                for body in &self.bodies {
                    let body_id = match body.as_str() {
                        "EARTH" => 399,
                        "MOON" => 301,
                        "SUN" => 10,
                        "MARS" => 499,
                        "JUPITER" => 599,
                        "VENUS" => 299,
                        _ => continue,
                    };
                    
                    let body_state = self.spice.get_state(body_id, t, 0)?;
                    let dist = self.distance_to_body(&state, &body_state);
                    let radius = self.get_body_radius(body)?;
                    
                    // Impacto si la distancia es menor que el radio + margen de seguridad
                    if dist < radius * 1.05 {
                        impact_detected = true;
                        break;
                    }
                }
                
                if impact_detected {
                    trajectory.push(state);
                    break;
                }
            }

            trajectory.push(state);
        }

        Ok(trajectory)
    }

    /// Calcula las derivadas del estado
    fn derivatives(&self, state: &[f64; 7]) -> Result<[f64; 6], SpiceError> {
        let pos = [state[0], state[1], state[2]];
        let vel = [state[3], state[4], state[5]];
        let et = state[6];

        let acc = self.nbody_acceleration(&pos, et)?;

        Ok([vel[0], vel[1], vel[2], acc[0], acc[1], acc[2]])
    }

    /// Aceleración de n-cuerpos usando SPICE
    fn nbody_acceleration(&self, pos: &[f64; 3], et: f64) -> Result<[f64; 3], SpiceError> {
        let mut acc = [0.0; 3];
        let km_to_m = self.km_to_m;

        for body in &self.bodies {
            let body_id = match body.as_str() {
                "SUN" => 10,
                "EARTH" => 399,
                "MOON" => 301,
                "MARS" => 499,
                "JUPITER" => 599,
                "VENUS" => 299,
                _ => continue,
            };

            // Obtener estado del cuerpo (posición en km, velocidad en km/s)
            let body_state = self.spice.get_state(body_id, et, 0)?;
            let mu = self.spice.get_mu(body_id)?; // m³/s²

            // Convertir posición de km a m
            let bx = body_state.position[0] * km_to_m;
            let by = body_state.position[1] * km_to_m;
            let bz = body_state.position[2] * km_to_m;

            // Vector desde el cuerpo a la nave (metros)
            let dx = pos[0] - bx;
            let dy = pos[1] - by;
            let dz = pos[2] - bz;
            
            let r = (dx*dx + dy*dy + dz*dz).sqrt();
            
            // Evitar singularidades (menos de 1 metro)
            if r < 1.0 {
                continue;
            }
            
            // Aceleración gravitacional: a = -μ * r / r³
            let g = mu / (r * r * r);
            acc[0] += g * dx;
            acc[1] += g * dy;
            acc[2] += g * dz;
        }

        Ok(acc)
    }

    /// Avanza el estado un paso
    fn step_state(&self, state: &[f64; 7], deriv: &[f64; 6], dt: f64) -> Result<[f64; 7], SpiceError> {
        let mut new_state = *state;
        for i in 0..6 {
            new_state[i] += deriv[i] * dt;
        }
        new_state[6] += dt;
        Ok(new_state)
    }

    /// Distancia entre la nave y un cuerpo (metros)
    fn distance_to_body(&self, state: &[f64; 7], body_state: &SpiceState) -> f64 {
        let dx = state[0] - body_state.position[0] * self.km_to_m;
        let dy = state[1] - body_state.position[1] * self.km_to_m;
        let dz = state[2] - body_state.position[2] * self.km_to_m;
        (dx*dx + dy*dy + dz*dz).sqrt()
    }

    /// Radio del cuerpo (metros)
    fn get_body_radius(&self, body: &str) -> Result<f64, SpiceError> {
        match body {
            "EARTH" => Ok(6.3781366e6),
            "MOON" => Ok(1.7374e6),
            "SUN" => Ok(6.957e8),
            "MARS" => Ok(3.3962e6),
            "VENUS" => Ok(6.0518e6),
            "JUPITER" => Ok(6.9911e7),
            _ => Ok(1.0),
        }
    }

    /// Obtiene referencia al contexto SPICE
    pub fn spice(&self) -> &SpiceContext {
        &self.spice
    }

    /// Obtiene referencia mutable al contexto SPICE
    pub fn spice_mut(&mut self) -> &mut SpiceContext {
        &mut self.spice
    }
}

// ============================================================================
// FUNCIONES DE ALTO NIVEL
// ============================================================================

/// Simula una misión completa con efemérides reales
/// 
/// # Argumentos
/// * `initial_state` - Estado inicial [x, y, z, vx, vy, vz] en metros y m/s
/// * `launch_time_et` - Tiempo de lanzamiento en ET (segundos desde J2000)
/// * `duration_days` - Duración de la simulación en días
/// * `config` - Configuración de efemérides
/// 
/// # Retorna
/// * `Vec<NBodyState>` - Trayectoria completa
pub fn simulate_mission_with_ephemeris(
    initial_state: [f64; 6],
    launch_time_et: f64,
    duration_days: f64,
    config: &super::EphemerisFullConfig,
) -> Result<Vec<NBodyState>, SpiceError> {
    use super::SpiceContext;
    
    // Crear contexto SPICE
    let mut spice = SpiceContext::new(config.spice.clone());
    spice.load_kernels()?;
    
    let integrator_config = crate::integration::IntegratorConfig {
        max_step: 3600.0 / 86400.0, // 1 hora en días
        ..config.integrator.clone()
    };
    
    let mut integrator = NBodyIntegrator::new(
        spice,
        integrator_config,
        config.bodies.clone(),
    );
    
    let duration_seconds = duration_days * 86400.0;
    let max_steps = (duration_seconds / 3600.0) as usize + 10; // 1 hora de paso
    
    integrator.propagate(
        &initial_state,
        launch_time_et,
        duration_seconds,
        max_steps,
    )
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integrator_creation() {
        let config = super::super::EphemerisFullConfig::default();
        let spice = SpiceContext::new(config.spice.clone());
        let integrator = NBodyIntegrator::new(
            spice,
            config.integrator,
            config.bodies.clone(),
        );
        assert!(!integrator.spice().is_loaded());
    }

    #[test]
    fn test_body_radius() {
        let config = super::super::EphemerisFullConfig::default();
        let spice = SpiceContext::new(config.spice.clone());
        let integrator = NBodyIntegrator::new(
            spice,
            config.integrator,
            config.bodies.clone(),
        );
        assert_eq!(integrator.get_body_radius("EARTH").unwrap(), 6.3781366e6);
        assert_eq!(integrator.get_body_radius("MOON").unwrap(), 1.7374e6);
        assert_eq!(integrator.get_body_radius("SUN").unwrap(), 6.957e8);
    }

    #[test]
    fn test_distance_to_body() {
        let config = super::super::EphemerisFullConfig::default();
        let spice = SpiceContext::new(config.spice.clone());
        let integrator = NBodyIntegrator::new(
            spice,
            config.integrator,
            config.bodies.clone(),
        );
        
        let state = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let body = SpiceState {
            position: [0.0, 0.0, 0.0],
            velocity: [0.0, 0.0, 0.0],
            light_time: 0.0,
            time_et: 0.0,
        };
        let dist = integrator.distance_to_body(&state, &body);
        assert_eq!(dist, 0.0);
    }
}
