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
    //const double dt_SI = 86400.0 / 3600.0;                        // One Day
    
    // Simulation Time in Seconds
    //const double simTime_SI = 1000.0 * 12.0 * dt_SI;
    const double simTime_SI = 2.0e13;
    // Timestep in Seconds
    const double dt_SI = simTime_SI / 1.0e6;            // steps

    // one lightyear in meters
    const double LJ_SI = 9.46e15;
    
    // output each nth timestep
    const int output = 1e3;
    
    // Number of Cells (no unit) - global
    const int cellsX = 5;
    const int cellsY = 8;

    // Standard Mass of a Particle in kg
    const double mass_SI = 1.98892e30 * 50.0; // Mass of the Sun

    const double distance_Sun_Neptun_SI = 4.5e12;
    const double distance_Sun_Earth_SI  = 1.496e11;
    // Force Cut-Off in Meters = Cell Size
    //const double cutoff_SI = distance_Sun_Earth_SI * 4.0;
    const double cutoff_SI = 20.0 * LJ_SI;
    const double cutoffIn_SI = 0.1 * cutoff_SI;

    // Gravitational Constant in m^3 / kg / s^2
    const double G_SI = 6.6738e-11;

    // Average radial speed of the Earth around the Sun in m/s
    const double EarthSpeed_SI = 2.978e4;
    // Average movement of Hyaden in m/s
    const double HyadenSpeed_SI = 43.0e3;
    
    // output process each n percent
    const double outPercent = 10;


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
    const double cutoffIn = cutoffIn_SI / UNIT_LENGTH;

    // Gravitation Constant
    const double G = G_SI * UNIT_MASS
      * UNIT_TIME * UNIT_TIME
      / UNIT_LENGTH / UNIT_LENGTH / UNIT_LENGTH;
    const double EarthSpeed = EarthSpeed_SI / UNIT_SPEED;
    const double HyadenSpeed = HyadenSpeed_SI / UNIT_SPEED;
    
    const double distance_Sun_Neptun = distance_Sun_Neptun_SI / UNIT_LENGTH;
    const double distance_Sun_Earth  = distance_Sun_Earth_SI  / UNIT_LENGTH;

    // Mass fo the Earth in Partials of Mass of the Sun
    const double partialEarthSun = 1.0 / 332946.0;
    
    // one lightyear in meters
    const double LJ = LJ_SI / UNIT_LENGTH;;


  } // namespace simParams
} // namespace MDSIM

#endif //SIMULATION_DEFINES_HPP