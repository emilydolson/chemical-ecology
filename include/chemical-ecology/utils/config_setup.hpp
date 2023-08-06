#pragma once

#include <filesystem>

#include "emp/config/ArgManager.hpp"
#include "chemical-ecology/Config.hpp"

namespace chemical_ecology::utils {

void use_existing_config_file(chemical_ecology::Config & config, emp::ArgManager & am) {
  if(std::filesystem::exists("chemical-ecology.cfg")) {
    std::cout << "Configuration read from chemical-ecology.cfg" << "\n";
    config.Read("chemical-ecology.cfg");
  }
  am.UseCallbacks();
  if (am.HasUnused())
    std::exit(EXIT_FAILURE);
}

void setup_config_native(chemical_ecology::Config & config, int argc, char* argv[]) {
  auto specs = emp::ArgManager::make_builtin_specs(&config);
  emp::ArgManager am(argc, argv, specs);
  use_existing_config_file(config, am);
}

} // End of chemical_ecology::utils namespace

