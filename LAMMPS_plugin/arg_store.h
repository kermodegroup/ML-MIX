#ifndef ARG_STORE_H
#define ARG_STORE_H

#include <vector>
#include <string>
#include <iostream>

class StoredArgs {
public:
    StoredArgs(int n, char **arg);
    StoredArgs(const StoredArgs &other) = default;
    StoredArgs(StoredArgs &&other) noexcept = default;
    StoredArgs &operator=(const StoredArgs &other) = default;
    StoredArgs &operator=(StoredArgs &&other) noexcept = default;
    ~StoredArgs() = default;
    char** toArgArray() const;

    void print() const;
    int nargs;
    std::vector<std::string> args;

private:

};

#endif // ARG_STORE_H