#ifndef ARGSMANAGER_H
#define ARGSMANAGER_H

#include "arg_store.h"
#include <map>
#include <string>

class ArgsManager {
public:
    std::map<std::string, StoredArgs> stored_args_map;

    // Store a new set of arguments with a key (e.g., pair style name)
    void storeArgs(const std::string &style, int narg, char **arg);

    // Retrieve a stored set by style name
    const StoredArgs &getArgs(const std::string &style) const;

    // Print all stored argument sets
    void printAll() const;
};

#endif // ARGSMANAGER_H
