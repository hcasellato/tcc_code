#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <map>

struct TimeData {
    std::vector<double> times_100;
    std::vector<double> times_1000;
    std::vector<double> times_10000;
    std::vector<double> times_100000;
};

void parseOutput(const std::string& filename, TimeData& data) {
    std::ifstream file(filename);
    std::string line;
    
    while (std::getline(file, line)) {
        size_t time_pos = line.find(" em ");
        if (time_pos == std::string::npos) continue;
        
        // Extract time value (comes after " em ")
        double time = std::stod(line.substr(time_pos + 4));
        
        // Determine which step this is
        size_t step_pos = line.find("(");
        if (step_pos == std::string::npos) continue;
        
        std::string step = line.substr(step_pos + 1, 6);
        if (step == "000100") data.times_100.push_back(time);
        else if (step == "001000") data.times_1000.push_back(time);
        else if (step == "010000") data.times_10000.push_back(time);
        else if (step == "100000") data.times_100000.push_back(time);
    }
    file.close();
}

double calculateAverage(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    double sum = 0.0;
    for (double val : values) sum += val;
    return sum / values.size();
}

void runProgramAnalysis(const std::string& program_name, int runs) {
    const std::string temp_output = "temp_output.txt";
    const std::string results_file = "time_averages_" + program_name + ".txt";
    TimeData timeData;
    
    std::cout << "Starting analysis for " << program_name << " (" << runs << " runs)\n";
    
    for (int i = 0; i < runs; ++i) {
        // Execute program and redirect output to temp file
        std::string command = "./" + program_name + " > " + temp_output;
        int ret = system(command.c_str());
        
        if (ret != 0) {
            std::cerr << "Run " << i+1 << " failed with code " << ret << std::endl;
            continue;
        }
        
        // Parse the output for time values
        parseOutput(temp_output, timeData);
        
        // Delete temp file
        remove(temp_output.c_str());
        
        if ((i+1) % 100 == 0) {
            std::cout << "Completed " << i+1 << " runs for " << program_name << "\n";
        }
    }
    
    // Calculate and save time averages
    std::ofstream out(results_file);
    out << std::fixed << std::setprecision(15);
    
    out << "Average Time Values (" << runs << " runs) - " << program_name << "\n";
    out << "=========================================\n\n";
    out << "100 steps:    " << calculateAverage(timeData.times_100) << " s\n";
    out << "1000 steps:   " << calculateAverage(timeData.times_1000) << " s\n";
    out << "10000 steps:  " << calculateAverage(timeData.times_10000) << " s\n";
    out << "100000 steps: " << calculateAverage(timeData.times_100000) << " s\n";
    
    out.close();
    std::cout << "\nTime averages saved to " << results_file << std::endl;
}

int main() {
    const std::vector<std::string> programs = {
        "1D_EDP_ChLDL_ex1",
        "1D_EDP_Ch_ex1",
        "1D_EDP_CF_ex1"
    };
    
    const int runs = 1000;
    
    for (const auto& program : programs) {
        runProgramAnalysis(program, runs);
    }
    
    std::cout << "\nAll program analyses completed!\n";
    return 0;
}
// g++ -o execute execute.cpp && ./execute