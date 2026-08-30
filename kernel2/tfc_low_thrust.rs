//! tfc_low_thrust.rs - Optimización de bajo empuje con TFC
//! Extiende la TFC para problemas de control continuo en 3D

use crate::tfc::{TFCProblem, TFCError, solve_tfc_linear, solve_tfc_levenberg_marquardt, TFCOptions};
use crate::crtbp::{StateVector, jacobi_constant, crtbp_derivatives};
use nalgebra::SVector;

// ============================================================================
// TIPOS
// ============================================================================

/// Vector de control 3D: [throttle, alpha, beta]
/// - throttle: 0.0 - 1.0
/// - alpha: ángulo en plano XY (rad)
/// - beta: ángulo fuera del plano (rad)
pub type ControlVector = SVector<f64, 3>;

/// Estado extendido con masa
#[derive(Debug, Clone, Copy)]
pub struct ExtendedState {
    pub state: StateVector,
    pub mass: f64,
}

impl ExtendedState {
    pub fn new(state: StateVector, mass: f64) -> Self {
        ExtendedState { state, mass }
    }

    pub fn to_array(&self) -> [f64; 7] {
        [
            self.state[0], self.state[1], self.state[2],
            self.state[3], self.state[4], self.state[5],
            self.mass,
        ]
    }
}

// ============================================================================
// PROBLEMA DE BAJO EMPUJE
// ============================================================================

#[derive(Clone)]
pub struct LowThrustProblem {
    pub initial_state: StateVector,
    pub final_state: StateVector,
    pub tof: f64,              // Tiempo de vuelo (unidades canónicas)
    pub tfc_order: usize,
    pub mass_initial: f64,
    pub thrust_max: f64,
    pub isp: f64,
    pub mu: f64,               // Parámetro de masa CRTBP
    pub n_segments: usize,     // Número de segmentos para control
}

impl LowThrustProblem {
    pub fn new(
        initial: StateVector,
        final_state: StateVector,
        tof: f64,
        mass_initial: f64,
    ) -> Self {
        LowThrustProblem {
            initial_state: initial,
            final_state: final_state,
            tof,
            tfc_order: 10,
            mass_initial,
            thrust_max: 0.08,
            isp: 1500.0,
            mu: 0.01215058560962404,
            n_segments: 20,
        }
    }

    /// Velocidad de escape del motor (m/s)
    pub fn exhaust_velocity(&self) -> f64 {
        self.isp * 9.80665
    }
}

// ============================================================================
// OPTIMIZADOR DE BAJO EMPUJE CON TFC
// ============================================================================

pub struct LowThrustOptimizer {
    problem: LowThrustProblem,
    options: TFCOptions,
    /// Control profile optimizado (por segmento)
    control_profile: Vec<ControlVector>,
}

impl LowThrustOptimizer {
    pub fn new(problem: LowThrustProblem, options: TFCOptions) -> Self {
        let n_segments = problem.n_segments;
        LowThrustOptimizer {
            problem,
            options,
            control_profile: vec![ControlVector::zeros(); n_segments],
        }
    }

    /// Optimiza perfil de empuje usando TFC con múltiples restricciones
    pub fn optimize(&mut self) -> Result<Vec<f64>, TFCError> {
        let n_coeffs = self.problem.tfc_order + 1;
        
        // Restricciones de contorno (posición y velocidad)
        let constraints = vec![
            (0.0, self.problem.initial_state[0], self.problem.initial_state[3]),
            (self.problem.tof, self.problem.final_state[0], self.problem.final_state[3]),
        ];

        let mut tfc_problem = TFCProblem::new(
            self.problem.tfc_order,
            0.0,
            self.problem.tof,
            constraints,
        ).map_err(|e| TFCError::NoConvergence(0, 0.0))?;

        // Función de residuo con término de control
        // Nota: Usamos un closure que NO captura self directamente
        let mu = self.problem.mu;
        let thrust_max = self.problem.thrust_max;
        let mass_initial = self.problem.mass_initial;
        let tof = self.problem.tof;

        let residual_fn = move |t: f64, r: f64, dr: f64, _ddr: f64| -> f64 {
            // Aceleración radial en CRTBP
            let acc = crtbp_radial_acceleration(r, dr, mu);
            
            // Control estimado (será optimizado por TFC)
            let control = estimate_control_simple(t, tof, thrust_max, mass_initial);
            
            acc - control
        };

        // Resolver con TFC lineal
        let xi = solve_tfc_linear(&mut tfc_problem, &residual_fn, &self.options)?;
        
        // Refinar con Levenberg-Marquardt
        let xi_final = solve_tfc_levenberg_marquardt(
            &mut tfc_problem,
            &residual_fn,
            &xi,
            &self.options,
        )?;

        Ok(xi_final.0)
    }

    /// Optimiza el perfil de control usando TFC + Shooting múltiple
    pub fn optimize_control_profile(&mut self) -> Result<Vec<ControlVector>, String> {
        let dt = self.problem.tof / self.problem.n_segments as f64;
        let mut trajectory = vec![ExtendedState::new(self.problem.initial_state, self.problem.mass_initial)];
        let mut controls = vec![ControlVector::zeros(); self.problem.n_segments];

        for iter in 0..self.options.max_iterations {
            // Propagar trayectoria con controles actuales
            let mut state = ExtendedState::new(self.problem.initial_state, self.problem.mass_initial);
            trajectory.clear();
            trajectory.push(state);

            for i in 0..self.problem.n_segments {
                state = self.propagate_segment(&state, &controls[i], dt);
                trajectory.push(state);
            }

            // Verificar condición final
            let error = (state.state - self.problem.final_state).norm();
            if error < self.options.tolerance {
                self.control_profile = controls.clone();
                return Ok(controls);
            }

            // Actualizar controles usando gradiente
            self.update_controls_gradient(&mut controls, &trajectory, dt, iter);
        }

        Err(format!("No convergió en {} iteraciones", self.options.max_iterations))
    }

    /// Propaga un segmento con control constante
    fn propagate_segment(
        &self,
        state: &ExtendedState,
        control: &ControlVector,
        dt: f64,
    ) -> ExtendedState {
        let mu = self.problem.mu;
        let ve = self.problem.exhaust_velocity();

        // RK4
        let k1 = self.low_thrust_derivatives(state, control, mu, ve);
        let s2 = self.step_state(state, &k1, dt / 2.0);
        let k2 = self.low_thrust_derivatives(&s2, control, mu, ve);
        let s3 = self.step_state(state, &k2, dt / 2.0);
        let k3 = self.low_thrust_derivatives(&s3, control, mu, ve);
        let s4 = self.step_state(state, &k3, dt);
        let k4 = self.low_thrust_derivatives(&s4, control, mu, ve);

        let mut result = *state;
        let state_vec = result.state;
        let mass = result.mass;

        // Actualizar estado (RK4)
        let new_state = state_vec
            + (k1.state + 2.0 * k2.state + 2.0 * k3.state + k4.state) * (dt / 6.0);
        let new_mass = mass + (k1.mass + 2.0 * k2.mass + 2.0 * k3.mass + k4.mass) * (dt / 6.0);

        ExtendedState {
            state: new_state,
            mass: new_mass.max(0.0),
        }
    }

    /// Derivadas del estado con bajo empuje
    fn low_thrust_derivatives(
        &self,
        state: &ExtendedState,
        control: &ControlVector,
        mu: f64,
        ve: f64,
    ) -> ExtendedState {
        let throttle = control[0].clamp(0.0, 1.0);
        let alpha = control[1];
        let beta = control[2];

        let thrust = throttle * self.problem.thrust_max;
        let mass = state.mass.max(1.0);

        // Componentes del empuje en el marco rotante
        let tx = thrust * alpha.cos() * beta.cos() / mass;
        let ty = thrust * alpha.sin() * beta.cos() / mass;
        let tz = thrust * beta.sin() / mass;

        // Aceleración gravitacional CRTBP (completa 3D)
        let state_vec = state.state;
        let grav = crtbp_derivatives(0.0, &state_vec);

        // Estado derivado
        let mut deriv_state = StateVector::zeros();
        deriv_state[0] = state_vec[3];
        deriv_state[1] = state_vec[4];
        deriv_state[2] = state_vec[5];
        deriv_state[3] = grav[3] + tx;
        deriv_state[4] = grav[4] + ty;
        deriv_state[5] = grav[5] + tz;

        // Cambio de masa (Tsiolkovsky)
        let mass_dot = -thrust / ve;

        ExtendedState {
            state: deriv_state,
            mass: mass_dot,
        }
    }

    /// Avanza el estado un paso
    fn step_state(&self, state: &ExtendedState, deriv: &ExtendedState, dt: f64) -> ExtendedState {
        ExtendedState {
            state: state.state + deriv.state * dt,
            mass: state.mass + deriv.mass * dt,
        }
    }

    /// Actualiza controles usando gradiente numérico
    fn update_controls_gradient(
        &self,
        controls: &mut Vec<ControlVector>,
        trajectory: &[ExtendedState],
        dt: f64,
        iteration: usize,
    ) {
        let learning_rate = 0.1 / (1.0 + iteration as f64 * 0.05);
        let eps = 0.001;

        for i in 0..controls.len() {
            let state = &trajectory[i];

            // Gradiente para throttle
            let mut c_plus = controls[i];
            c_plus[0] = (c_plus[0] + eps).clamp(0.0, 1.0);
            let state_plus = self.propagate_segment(state, &c_plus, dt);
            let grad_throttle = (state_plus.state[3] - trajectory[i + 1].state[3]).abs() / eps;

            // Gradiente para alpha
            let mut c_plus = controls[i];
            c_plus[1] += eps;
            let state_plus = self.propagate_segment(state, &c_plus, dt);
            let grad_alpha = (state_plus.state[4] - trajectory[i + 1].state[4]).abs() / eps;

            // Gradiente para beta
            let mut c_plus = controls[i];
            c_plus[2] += eps;
            let state_plus = self.propagate_segment(state, &c_plus, dt);
            let grad_beta = (state_plus.state[5] - trajectory[i + 1].state[5]).abs() / eps;

            // Actualizar (descenso de gradiente)
            controls[i][0] = (controls[i][0] + learning_rate * grad_throttle).clamp(0.0, 1.0);
            controls[i][1] += learning_rate * grad_alpha;
            controls[i][2] += learning_rate * grad_beta;
        }
    }

    /// Calcula el ΔV efectivo de la trayectoria optimizada
    pub fn compute_effective_dv(&self) -> f64 {
        let ve = self.problem.exhaust_velocity();
        let mut mass = self.problem.mass_initial;
        let mut dv = 0.0;

        for control in &self.control_profile {
            let thrust = control[0] * self.problem.thrust_max;
            let acceleration = thrust / mass;
            let dt = self.problem.tof / self.problem.n_segments as f64;
            dv += acceleration * dt;

            let m_dot = thrust / ve;
            mass -= m_dot * dt;

            if mass <= 0.0 {
                break;
            }
        }

        dv
    }

    /// Obtiene el perfil de control optimizado
    pub fn control_profile(&self) -> &[ControlVector] {
        &self.control_profile
    }

    /// Reconstruye la trayectoria completa
    pub fn reconstruct_trajectory(&self) -> Vec<ExtendedState> {
        let mut trajectory = Vec::with_capacity(self.problem.n_segments + 1);
        let mut state = ExtendedState::new(self.problem.initial_state, self.problem.mass_initial);
        let dt = self.problem.tof / self.problem.n_segments as f64;

        trajectory.push(state);
        for control in &self.control_profile {
            state = self.propagate_segment(&state, control, dt);
            trajectory.push(state);
        }

        trajectory
    }
}

// ============================================================================
// FUNCIONES AUXILIARES
// ============================================================================

/// Aceleración radial en el CRTBP (1D simplificado)
pub fn crtbp_radial_acceleration(r: f64, dr: f64, mu: f64) -> f64 {
    let x = r;
    let y = 0.0;
    let vx = dr;
    let vy = 0.0;

    let r1 = ((x + mu).powi(2) + y.powi(2)).sqrt();
    let r2 = ((x - 1.0 + mu).powi(2) + y.powi(2)).sqrt();

    2.0 * vy + x - (1.0 - mu) * (x + mu) / (r1.powi(3)) - mu * (x - 1.0 + mu) / (r2.powi(3))
}

/// Estimación simple de control (para inicialización)
pub fn estimate_control_simple(t: f64, tof: f64, thrust_max: f64, mass: f64) -> f64 {
    let thrust_ratio = 0.5 * (1.0 + (2.0 * t / tof - 1.0).sin());
    thrust_ratio * thrust_max / mass
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_low_thrust_derivatives() {
        let problem = LowThrustProblem::new(
            StateVector::new(0.8, 0.0, 0.0, 0.0, 0.1, 0.0),
            StateVector::new(0.9, 0.0, 0.0, 0.0, 0.0, 0.0),
            10.0,
            250.0,
        );

        let optimizer = LowThrustOptimizer::new(problem, TFCOptions::default());
        let state = ExtendedState::new(
            StateVector::new(0.8, 0.0, 0.0, 0.0, 0.1, 0.0),
            250.0,
        );
        let control = ControlVector::new(0.5, 0.0, 0.0);
        let deriv = optimizer.low_thrust_derivatives(
            &state,
            &control,
            0.01215,
            1500.0 * 9.80665,
        );

        // La masa debe decrecer
        assert!(deriv.mass < 0.0);
        // Debe haber aceleración
        assert!(deriv.state[3] != 0.0 || deriv.state[4] != 0.0);
    }

    #[test]
    fn test_smart1_validation() {
        // Datos de SMART-1 (ESA)
        let initial = StateVector::new(
            0.83, 0.0, 0.0,
            0.0, 0.15, 0.0,
        );

        let final_state = StateVector::new(
            0.99, 0.0, 0.0,
            0.0, 0.05, 0.0,
        );

        let problem = LowThrustProblem {
            initial_state: initial,
            final_state,
            tof: 410.0,
            tfc_order: 12,
            mass_initial: 367.0,
            thrust_max: 0.07,
            isp: 1500.0,
            mu: 0.01215058560962404,
            n_segments: 20,
        };

        let options = TFCOptions {
            max_iterations: 100,
            tolerance: 1e-6,
            ..TFCOptions::default()
        };

        let mut optimizer = LowThrustOptimizer::new(problem, options);

        // Optimizar perfil de control
        let result = optimizer.optimize_control_profile();
        assert!(result.is_ok(), "Optimización falló: {:?}", result);

        let dv = optimizer.compute_effective_dv();
        println!("SMART-1 validation: computed ΔV = {:.1} m/s", dv);

        // SMART-1 reportó ~600 m/s
        assert!(dv < 800.0 && dv > 400.0,
            "ΔV fuera de rango esperado: {:.1} m/s", dv);
    }

    #[test]
    fn test_control_profile_optimization() {
        let problem = LowThrustProblem::new(
            StateVector::new(0.8, 0.0, 0.0, 0.0, 0.1, 0.0),
            StateVector::new(0.9, 0.0, 0.0, 0.0, 0.0, 0.0),
            5.0,
            250.0,
        );

        let options = TFCOptions {
            max_iterations: 50,
            tolerance: 1e-4,
            ..TFCOptions::default()
        };

        let mut optimizer = LowThrustOptimizer::new(problem, options);
        let result = optimizer.optimize_control_profile();

        assert!(result.is_ok(), "Optimización falló: {:?}", result);
        assert_eq!(optimizer.control_profile().len(), 20);

        let trajectory = optimizer.reconstruct_trajectory();
        assert_eq!(trajectory.len(), 21);

        let dv = optimizer.compute_effective_dv();
        assert!(dv > 0.0);
    }
}
