#include <iostream>
#include <vector> 

double calculate_average(const std::vector<double>& numbers) {// const = do not change numbers
    if (numbers.empty()) {
        return 0.0;
    }

    double sum = 0.0;
    for (int i = 0; i < numbers.size(); i++){
        sum += numbers[i];
    }

    return sum / numbers.size();
}

int main(){
    std::vector<double> numbers; // use vectors for arrays, like list 
    numbers.reserve(10); // save space for 10 elements
    for (int  i = 0; i < 10; i++){
        numbers[i] = 0;
    }
    double avg = calculate_average(numbers);
    std::cout << avg << std::endl;
    return 0;
}

// take aways, use passby reference wereever possible