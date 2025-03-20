#include "arg_store.h"
#include <cstring>
#include <iostream>
StoredArgs::StoredArgs(int n, char **arg) : nargs(n), args(n) {
    for (int i = 0; i < n; ++i) {
        args[i] = arg[i]; // Automatically handles deep copy
    }
}

void StoredArgs::print() const {
    std::cout << "Number of arguments: " << nargs << std::endl;
    for (const auto &s : args) {
        std::cout << s << " ";
    }
    std::cout << std::endl;
}

char** StoredArgs::toArgArray() const {
    char **argArray = new char*[nargs];
    for (int i = 0; i < nargs; ++i) {
        argArray[i] = new char[args[i].size() + 1];
        std::strcpy(argArray[i], args[i].c_str());
    }
    return argArray;
}