// g++ -c extern.cpp -o extern.o
// nm extern.o

// Will be mangled by C++ compiler
void process_data(int data) {
    // Dummy implementation
    volatile int x = data;
}

// Overloaded version - will have a DIFFERENT mangled name
void process_data(double data) {
    // Dummy implementation
    volatile double y = data;
}

// *** Using extern "C" ***
// Tells C++ compiler NOT to mangle this name
extern "C" void process_data_for_c(int data) {
    // Dummy implementation
    volatile int z = data; //ignore volatile
}
