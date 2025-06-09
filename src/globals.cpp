// globals.cpp
#include "globals.hpp"

std::mutex dataMutex;
std::atomic<bool> simulationRunning(true);