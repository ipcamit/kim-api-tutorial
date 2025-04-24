#include <iostream>

struct Data { int a; double b; };

int main() {
  Data myData = {10, 3.14};
  // Treat the Data object's memory as raw bytes
  char* bytePtr = reinterpret_cast<char*>(&myData);

  std::cout << "First few bytes of Data object:\n";
  for (int i = 0; i < sizeof(Data); ++i) {
    std::cout << std::hex << (int)(unsigned char)bytePtr[i] << " ";
  }
  std::cout << std::dec << std::endl;

  // DANGEROUS: Converting unrelated pointer types
  long addr = reinterpret_cast<long>(&myData);
  std::cout << "Address stored as long: " << addr << std::endl;

  return 0;
}
