#ifndef SIMULATION_DEFINES_HPP
#define SIMULATION_DEFINES_CPP

class units;

// ------------- Simulation parameters ---------------- //
class simulationParameters {
  friend class units;
private:
  // Timestep in Seconds
  constexpr static double simTime_SI = 400.0;
  constexpr static double dt_SI = 2.0;
public:
  constexpr inline double get_dt() const;// { return dt_SI / units::UNIT_TIME; } // / units::UNIT_TIME;
  //constexpr static double dt2 = dt_SI / units::UNIT_TIME;
} simGalaxy;


// ------------- Units -------------------------------- //
class units {
  friend class simulationParameters;
public:
  constexpr static double UNIT_TIME = simulationParameters::dt_SI;
  
} unitsGalaxy;

constexpr inline double simulationParameters::get_dt() const { return dt_SI / units::UNIT_TIME; }



// ------------- Unitless variables ------------------- //
//const double dt = dt_SI / UNIT_TIME;




#endif //SIMULATION_DEFINES_CPP