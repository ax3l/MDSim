#ifndef SIMULATION_DEFINES_HPP
#define SIMULATION_DEFINES_HPP


namespace MDSIM
{
  namespace simParams
  {
    // ------------- Simulation parameters ----------------------------------- //

    // PARAMETERS in SI Units ------------------------------------------------ //
    // Timestep in Seconds
    //const double dt_SI = 2.592e6;                      // One Month
    const double dt_SI = 86400.0;               // One Day
    // Simulation Time in Seconds
    //const double simTime_SI = 1000.0 * 12.0 * dt_SI;
    const double simTime_SI = 365.0 * 86400.0;

    // Number of Cells (no unit) - global
    const int cellsX = 1;
    const int cellsY = 1;

    // Standard Mass of a Particle in kg
    const double mass_SI = 1.98892e30; // Mass of the Sun

    const double distance_Sun_Neptun_SI = 4.5e12;
    const double distance_Sun_Earth_SI  = 1.496e11;
    // Force Cut-Off in Meters = Cell Size
    const double cutoff_SI = distance_Sun_Earth_SI * 4.0;

    // Gravitational Constant in m^3 / kg / s^2
    const double G_SI = 6.6738e-11;

    // Average radial speed of the Earth around the Sun in m/s
    const double EarthSpeed_SI = 2.978e4;


    // UNITS ----------------------------------------------------------------- //
    const double UNIT_TIME = dt_SI;
    const double UNIT_MASS = mass_SI;
    const double UNIT_LENGTH = cutoff_SI;
    const double UNIT_SPEED = UNIT_LENGTH / UNIT_TIME;

    // PARAMETERS - UNITLESS ------------------------------------------------- //
    const double simTime = simTime_SI / UNIT_TIME;
    const double dt = dt_SI / UNIT_TIME;
    const double mass = mass_SI / UNIT_MASS;
    const double cutoff = cutoff_SI / UNIT_LENGTH;

    // Gravitation Constant
    const double G = G_SI * UNIT_MASS
      * UNIT_TIME * UNIT_TIME
      / UNIT_LENGTH / UNIT_LENGTH / UNIT_LENGTH;
    const double EarthSpeed = EarthSpeed_SI / UNIT_SPEED;
    
    const double distance_Sun_Neptun = distance_Sun_Neptun_SI / UNIT_LENGTH;
    const double distance_Sun_Earth  = distance_Sun_Earth_SI  / UNIT_LENGTH;

    // Mass fo the Earth in Partials of Mass of the Sun
    const double partialEarthSun = 1.0 / 332946.0;


  } // namespace simParams
} // namespace MDSIM

#endif //SIMULATION_DEFINES_HPP