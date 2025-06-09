// globals.hpp
#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <mutex>
#include <atomic>

extern std::mutex dataMutex;
extern  std::atomic<bool> simulationRunning;


#endif // GLOBALS_HPP