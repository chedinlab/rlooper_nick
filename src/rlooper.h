//
// Created by Robert Stolz on 6/27/17.
// Modified & maintained by Stella Hartono 2023-2024
//

#include <iostream>
#include <cstdio> 
#include <chrono>
#include "simulation.h"
#include <string.h>
#include <unistd.h>
#include "mytime.h"
#include "math.h"
#include <stdexcept>
#include <cstring>
#include <string>
#include <unistd.h> // For fork, exec, pipe, and dup2
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
using namespace std;
using namespace std;

typedef std::chrono::high_resolution_clock::time_point TimeVar;
typedef std::string String;

int main(int argc, char* argv[]);
void printUsage(std::string& myscript);
void printUsageLong();
bool myverbose2 = false;
std::string getDefaultOutfile(int argc, char* argv[], std::map<std::string, std::string>& myargs, std::map<int, string>& myorder, string &runScript, stringstream& debugprint);
std::string turnSigmaIntoNiceStringForOutput(double sigma, stringstream& debugprint);
int parse_argument(int argc, char* argv[], Simulation& sim, Rloop_model& model, bool &dynamic_flag, bool &orig_flag, bool &sandbox);
void splitPath(const std::string& fullPath, std::string& folder, std::string& fileName);
int countCharInString(String s, char delim);
std::string exec(const char* cmd);
std::string chomp(string& text);
void parseSettingsFile(const std::string& filename, 
					   std::map<int, std::string>& myorder, std::map<std::string, std::string>& myargs, 
					   std::stringstream& debugprint);
// std::string getFilename(const std::string& filePath);
// std::string getFolder(const std::string& filePath);
// void getFolderfull(const std::string& filePath, std::string& folder, std::string& fileName);
// std::string getExtension(const std::string& filePath);
std::string getFilename(std::string& filePath);
std::string getFolder(std::string& filePath);
void getFolderfull(std::string& filePath, std::string& folder, std::string& fileName);
std::string getExtension(std::string& filePath);
std::string whichFullPath(string& filePath);
std::string getFullPath(string& filePath);
std::string trimSpace(std::string& mystring);
inline bool fileexists (const std::string& file);
bool atob(const string& str);