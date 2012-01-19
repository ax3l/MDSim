#ifndef SIMULATION_DEFINES_HPP
#define SIMULATION_DEFINES_HPP


namespace MDSIM {
  namespace simParams {
  // ------------- Simulation parameters ----------------------------------- //

  // PARAMETERS in SI Units ------------------------------------------------ //
  // Timestep in Seconds
  const double dt_SI = 2.592e6;                      // One Month
  // Simulation Time in Seconds
  const double simTime_SI = 1000.0 * 12.0 * dt_SI;
  
  // Number of Cells (no unit) - global
  const int cellsX = 400;
  const int cellsY = 500;
  
  // Standard Mass of a Particle in kg
  const double mass_SI = 1.98892e30;                 // Mass of the Sun
  // Force Cut-Off in Meters = Cell Size
  const double distance_Sun_Neptun_SI = 4.5e12;
  const double cutoff_SI = distance_Sun_Neptun_SI;
  
  // Gravitational Constant in m^3 / kg / s^2
  const double G_SI = 6.6738e-11;


  // UNITS ----------------------------------------------------------------- //
  const double UNIT_TIME = dt_SI;
  const double UNIT_MASS = mass_SI;
  const double UNIT_LENGTH = cutoff_SI;
  
  // PARAMETERS - UNITLESS ------------------------------------------------- //
  const double simTime = simTime_SI / UNIT_TIME;
  const double dt      = dt         / UNIT_TIME;
  const double mass    = mass_SI    / UNIT_MASS;
  const double cutoff  = cutoff_SI  / UNIT_LENGTH;
  
  const double G       = G_SI       * UNIT_MASS
                                    * UNIT_TIME * UNIT_TIME
                                    / UNIT_LENGTH / UNIT_LENGTH / UNIT_LENGTH;


  } // namespace simParams
} // namespace MDSIM

#endif //SIMULATION_DEFINES_HPP