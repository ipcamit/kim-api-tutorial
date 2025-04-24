#include <iostream>
#include <string>

// Function template to find the maximum of two values
template <typename T> // "T" is a placeholder for any type
T maximum(T a, T b) {
  return (a > b) ? a : b;
}

int main() {
  std::cout << "Max(5, 10): " << maximum(5, 10) << std::endl; // T is int
  std::cout << "Max(3.14, 2.71): " << maximum(3.14, 2.71) << std::endl; // T is double
  std::cout << "Max('a', 'z'): " << maximum('a', 'z') << std::endl; // T is char
  // std::cout << "Max(5, 3.14): " << maximum(5, 3.14); // Error! T cannot be both int and double
  return 0;
}
