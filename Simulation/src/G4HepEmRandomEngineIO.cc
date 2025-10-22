// Simulation/src/G4HepEmRandomEngineIO.cc
#include "G4HepEmRandomEngine.hh"
#include <fstream>
#include <stdexcept>
#include <string>
#include <random>

void G4HepEmRandomEngine::SaveStateToFile(const std::string& path) {
  auto* eng = static_cast<std::mt19937_64*>(fObject);   // or std::mt19937*
  std::ofstream os(path);
  if (!os) throw std::runtime_error("Failed to open file: " + path);
  os << *eng;   // text format, portable
}

//void G4HepEmRandomEngine::LoadStateToFile(const std::string& path); // (typo guard)
void G4HepEmRandomEngine::LoadStateFromFile(const std::string& path) {
  auto* eng = static_cast<std::mt19937_64*>(fObject);   // or std::mt19937*
  std::ifstream is(path);
  if (!is) throw std::runtime_error("Failed to open file: " + path);
  is >> *eng;
}

void G4HepEmRandomEngine::Seed(uint64_t seed) {
  auto* eng = static_cast<std::mt19937_64*>(fObject);   // or std::mt19937*
  eng->seed(seed);
}