#include <iostream>
#include <string>   // For using std::string

int main() { // Main entry point
    int count = 5; // integer
    if (count > 0) { // if condition
        std::cout << "Count is positive: " << count << std::endl; // Print to console
    }

    for (int i = 0; i < count; ++i) { // for i in range(count)
         std::cout << "i = " << i << std::endl;
    }

    // prefer std::string for simplicity
    std::string message = "mystring";
    std::cout << message << std::endl;

    // division between integer is "//" in C++
    std::cout << 1/2 << std::endl; // 0, not 0.5

    return 0; // Indicates successful execution
}