#include <iostream>

class Thing {
public:
    // 1. Static Member Variable (Declaration)
    static int count;

    Thing() {
        ++count; std::cout << "Thing created! (Total things: " << count << ")" << std::endl;
    }
    ~Thing() {
        --count; std::cout << "Thing destroyed! (Total things: " << count << ")" << std::endl;
    }

    // 2. Static Member Function
    static int howMany() {
        return count; // Accesses the static 'count'
    }
};

// Definition and Initialization of the static member variable outside the class
int Thing::count = 0;

int main() {
    std::cout << "Initial count: " << Thing::howMany() << std::endl; // Call static function
    Thing t1; Thing t2;
    std::cout << "Current count: " << Thing::howMany() << std::endl;
    {Thing t3; // Create third object in a limited scope
     std::cout << "Count inside scope: " << Thing::howMany() << std::endl;} // t3 is destroyed here, destructor runs
    std::cout << "Count after scope: " << Thing::count << std::endl; // Can also access directly (if public)
    return 0;
}