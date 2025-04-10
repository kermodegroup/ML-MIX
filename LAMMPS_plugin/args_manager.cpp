#include "args_manager.h"
#include <iostream>

// Store a new set of arguments
void ArgsManager::storeArgs(const std::string &style, int narg, char **arg) {
    stored_args_map.emplace(style, StoredArgs(narg, arg)); // Store deep copy
}

// Retrieve a stored set by style name
const StoredArgs &ArgsManager::getArgs(const std::string &style) const {
    return stored_args_map.at(style);
}

// Print all stored argument sets
void ArgsManager::printAll() const {
    for (const auto &pair : stored_args_map) {
        pair.second.print();
    }
}
