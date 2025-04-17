#include <iostream>
#include <cstdio> 
#include <chrono>
#include "simulation.h"
#include <string.h>
#include <unistd.h>
#include "rlooper.h"
#include "math.h"
#include <stdexcept>
#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>
#include <unistd.h> // For fork, exec, pipe, and dup2
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <algorithm>
#include <cctype>
using namespace std;

typedef std::chrono::high_resolution_clock::time_point TimeVar;
typedef std::string String;



int main(int argc, char* argv[]) {

	const char *LPR = "\033[1;95m";
	const char *LBU = "\033[1;94m";
	const char *LRD = "\033[1;31m";
	const char *LCY = "\033[1;36m";
	const char *YW = "\033[1;33m";
	const char *LGN = "\033[1;32m";
	const char *NN = "\033[1;0m";


	TimeVar t1=timeNow();
	Simulation sim;
	Rloop_model model;
	bool dynamic_flag = false;
	bool orig_flag = false;
	bool sandbox = false;
	sim.set_seed(42);
	
	if (parse_argument(argc, argv, sim, model, dynamic_flag, orig_flag, sandbox) == 1) {return 1;}
	// Simulation
	//if (dynamic_flag) { // NOT IMPLEMENTED
	//	cout << YW << "\n>Runnning simulation A (dynamic): \n";
	//	sim.simulation_A();
	//}
	cout << LPR << "\nrlooper.cpp::" << "main" << NN << "\n";
	cout << LPR << "\nrlooper.cpp::" << "parse_argument" << NN << "\n";
	cout << YW << "\n>Running simulation: \n";
	sim.simulation_A();

	cout << YW << "\n>Simulation done!\n" << NN << "\n";
	cout << YW << "\nTotal runtime: " << LGN << duration(timeNow()-t1)/1e9 << YW << " seconds\n" << NN << "\n";
	return 0;
}

std::string getDefaultOutfile(int argc, char* argv[], std::map<std::string, std::string>& myargs, std::map<int, std::string>& myorder, string &runScript, stringstream& debugprint) {
	const char *LPR = "\033[1;95m";
	const char *LBU = "\033[1;94m";
	const char *LRD = "\033[1;31m";
	const char *LCY = "\033[1;36m";
	const char *YW = "\033[1;33m";
	const char *LGN = "\033[1;32m";
	const char *NN = "\033[1;0m";
	
	
	std::map<std::string, std::string> currargs;
	std::stringstream ssoutfile;
	std::string myscriptname = string(argv[0]);
	for (int i = 1; i < argc; i++) {
		// checkArgumentName(argv[i]);
		std::string argkey = argv[i];
		argkey = argkey.substr(2, argkey.length());
		if (argkey == "h") {
			printUsage(myscriptname);
			exit(1);
		}
		else if (argkey == "H") {
			printUsage(myscriptname);
			printUsageLong();
			exit(1);
		}
		else if (argkey == "i") {
			argkey = "fasta";
		}
		else if (argkey == "o") {
			argkey = "outname";
		}
		else if (argkey == "D") {
			argkey = "dynamic";
		}
		else if (argkey == "w") {
			argkey = "window";
		}
		else if (argkey == "n") {
			argkey = "nick";
		}
		
			
		if (i == argc - 1 || (i < argc - 1 && string(argv[i+1]).substr(0, 2) == "--")) {
			currargs[argkey] = "True";
			runScript += " --" + argkey;
			debugprint << LPR << "rlooper::getDefaultOutfile:" << LBU << argkey << "=" << LGN << "True" << NN << endl;
			//ssoutfile << "_" << argkey;
		}
		else {
			string argval = argv[i+1];
			chomp(argval);
			currargs[argkey] = argval;
			debugprint << LPR << "rlooper::getDefaultOutfile:" << LCY << argkey << "=" << LGN << argval << NN << endl;
			if (argkey == "sigma") {
				argval = turnSigmaIntoNiceStringForOutput(atof(argv[i + 1]),debugprint);
			}
			else if (argkey == "fasta") {
				currargs["fasta"] = getFullPath(currargs["fasta"]);
				chomp(currargs["fasta"]);
			}
			else if (argkey != "fasta" && argkey != "outname" && argkey != "outdir") {
				currargs[argkey] = string(argv[i+1]);
				argval = string(argv[i+1]);
			}
			if (argkey != "fasta" && argkey != "outname" && argkey != "outdir") {
				debugprint << LRD << "added" << argkey << "," << argval << NN << "\n";
				// ssoutfile << "_" << argkey << "," << argval;
			}
			runScript += " --" + argkey + " " + argval;

			i ++;
		}
	}

	
	for (map<string,string>::iterator it = currargs.begin(); it != currargs.end(); it++) {
		const string& argkey = it->first;
		const string& argval = it->second;
		debugprint << argkey << ":" << argval << endl;
	}	
	string argUsed = "";
	string argDefaultUsed = "";
	string argDefaultNotUsed = "";
	for (const auto& pair : myorder) {
		const string& argkey = pair.second;
		// const string& val = myargs[argkey];
		std::string argval = myargs[argkey];
		if (currargs.count(argkey) > 0) {
			string currval = currargs[argkey];
			if (argkey == "sigma") {
				currval = turnSigmaIntoNiceStringForOutput(atof(currval.c_str()),debugprint);
			}
			if (argkey != "fasta" && argkey != "outname" && argkey != "outdir" && argkey != "verbose" && argkey != "debug") {
				if (currval == "True") {
					// cout << LPR << "rlooper::getDefaultOutfile:" << LCY << argkey << "=" << LGN << currval << NN << endl;
					// cout << argkey << "=" << "true" << endl;
					ssoutfile << "_" << argkey;
				}
				else {
					// cout << LPR << "rlooper::getDefaultOutfile:" << LCY << argkey << "=" << LGN << currval << NN << endl;
					// cout << argkey << "=" << currval << endl;
					ssoutfile << "_" << argkey << "," << currval;
				}
			}
			if (currval == "True") {
				argUsed += " --" + argkey;
			}
			else {
				argUsed += " --" + argkey + " " + currargs[argkey];
			}
			cout << LPR << "rlooper::getDefaultOutfile:" << LCY << argkey << "=" << LGN << currargs[argkey] << NN << " (default=" << LPR << myargs[argkey] << NN << ")" << endl;

			myargs[argkey] = currargs[argkey];
		}
		else if (argval == "False") {
			argDefaultNotUsed += " " + argkey;		
		}
		else if (argval == "True") {
			argDefaultUsed += " " + argkey;
		}
		else {
			argDefaultUsed += " " + argkey + " " + argval;
		}
	}
	
	
	if (!atob(myargs["dynamic"].c_str()) && !atob(myargs["orig"].c_str())) {
		myargs["naive"] = "True";
	}

	cout << "verbose=" << myargs["verbose"] << endl;
	if (myargs["fasta"] == "False" || !fileexists(myargs["fasta"])) {
		printUsage(myscriptname);
		cout << LRD << "Error: Fasta file doesn't exist! " << NN << "(" << LCY << myargs["fasta"] << NN << ")" << endl << endl;
		if (atob(myargs["verbose"].c_str())) {
			cout << "Runscript:\n" << runScript << "\n" << debugprint.str() << "\n\n";
		}
		exit(1);
	}

	if (!atob(myargs["outdir"])) {
		myargs["outdir"] = getFolder(myargs["fasta"]);
		cout << LPR << "rlooper::getDefaultOutfile:" << LCY << "outdir" << "=" << LGN << myargs["outdir"] <<  NN << endl;
	}

	if (!atob(myargs["outname"])) {
		myargs["outname"] = getFilename(myargs["fasta"]);
		cout << LPR << "rlooper::getDefaultOutfile:" << LCY << "outname" << "=" << LGN << myargs["outname"] <<  NN << endl;
	}

	
	myargs["outfile"] = myargs["outdir"] + "/" + myargs["outname"] + ssoutfile.str();
	cout << LPR << "rlooper::getDefaultOutfile:" << LCY << "outfile" << "=" << LGN << myargs["outfile"] <<  NN << endl;
	
	debugprint << LPR << "rlooper::getDefaultOutfile:" << LCY << "outfile_addition" << "=" << LGN << ssoutfile.str() << NN << endl;
	debugprint << YW << "Main Arguments:" << LCY << argUsed << NN << "\n" <<
	              YW << "Default Arguments:" << LGN << argDefaultUsed << NN << "\n" << 
	              YW << "Unused Arguments:" << LPR << argDefaultNotUsed << NN << "\n\n";
	              

	return myargs["outfile"];
}

int parse_argument(int argc, char* argv[], Simulation& sim, Rloop_model& model, bool &dynamic_flag, bool &orig_flag, bool &sandbox) {
	const char *LPR = "\033[1;95m";
	const char *LBU = "\033[1;94m";
	const char *LRD = "\033[1;31m";
	const char *LCY = "\033[1;36m";
	const char *YW = "\033[1;33m";
	const char *LGN = "\033[1;32m";
	const char *NN = "\033[1;0m";
	
	
	
	std::stringstream debugprint;
	std::string inputFile; // = "inputFile";
	std::string outfile; // = "outfile";
	std::string outdir; // = "./workflow/examples/"
	std::string myscriptname = argv[0];
	std::string myscript = whichFullPath(myscriptname); chomp(myscript);
	myscriptname = getFilename(myscript);
	
	std::stringstream ssoutfile;
	std::map<int, std::string> myorder;        // <line_number, key>
	std::map<std::string, std::string> myargs; // <key, value>

    parseSettingsFile("settings.ini", myorder, myargs, debugprint);

	std::stringstream longarg;
	std::string runScript;
	getDefaultOutfile(argc, argv, myargs, myorder, runScript, longarg);
	
	//cout << "fasta:" << myargs["fasta"] << ", atob:" << atob(myargs["fasta"].c_str()) << "\n";
	//cout << "\n" << YW << myscript;
	// std::string argUsed;
	// std::string argDefaultNotUsed;

	if (atob(myargs["help"].c_str())) {
		printUsage(myscriptname);
		if (atob(myargs["verbose"].c_str())) {
			cout << "Runscript:\n" << runScript << "\n" << debugprint.str() << "\n\n";
		}
		exit(1);
	}
	if (myargs.count("H") > 0) {
	//if (atob(myargs["H"].c_str())) {
		printUsageLong();
		exit(1);
	}

	sim.set_infile(myargs["fasta"].c_str());
	sim.set_outdir(myargs["outdir"].c_str());
	sim.set_outfile(myargs["outfile"].c_str());
	std::string mycmd = "mkdir -p " + myargs["outdir"];
	exec(mycmd.c_str());

	// More details about available options and their descriptions
	sim.setGeneBeg(atoi(myargs["beg"].c_str()));
	sim.setGeneEnd(atoi(myargs["end"].c_str()));
	sim.set_seed(atoi(myargs["seed"].c_str()));

	model.seta(atof(myargs["a"].c_str()));
	model.setnick(atoi(myargs["nick"].c_str()));
	model.setnicklen(atoi(myargs["nicklen"].c_str()));
	model.setselffoldlen(atoi(myargs["selffoldlen"].c_str()));
	model.setN(atoi(myargs["N"].c_str()));
	sim.set_auto_domain_size(atoi(myargs["N"].c_str()));
	model.set_superhelicity(atof(myargs["sigma"].c_str()));
	model.dynamic_setTx_Amb_Sigma(atof(myargs["sigma"].c_str()));
	sim.set_minlength(atoi(myargs["minlength"].c_str()));
	model.dynamic_setWindow_size(atoi(myargs["window"].c_str()));
	sim.dynamic_setWindow_size(atoi(myargs["window"].c_str()));
	model.set_bp_energy_override(atof(myargs["homopolymer"].c_str()));
	sim.reverse_input(atob(myargs["reverse"].c_str()));
	sim.complement_input(atob(myargs["complement"].c_str()));
	sim.reverse_input(atob(myargs["revcomp"].c_str()));
	sim.complement_input(atob(myargs["revcomp"].c_str()));
	
	model.set_unconstrained(atob(myargs["unconstrained"].c_str()));
	sim.set_verbose(atob(myargs["verbose"].c_str()));
	sandbox = atob(myargs["sandbox"].c_str());

	// sim.set_average_g(true);

	sim.set_printbedfile(atob(myargs["bedfile"].c_str()));
	sim.set_circular(atob(myargs["circular"].c_str()));
	sim.set_residuals(atob(myargs["residuals"].c_str()));
	sim.set_power_threshold(atoi(myargs["sensitivity"].c_str()));
	sim.set_top(atoi(myargs["top"].c_str()));
	sim.set_dump(atob(myargs["dump"].c_str()));
	sim.set_import_flag(myargs["import"].c_str());
	sim.set_circular(atob(myargs["circular"].c_str()));
	
	// original rlooper 
	sim.set_average_g(atob(myargs["localaverageenergy"].c_str()));
	
	
	sim.set_naive(atob(myargs["naive"].c_str()));
	sim.set_orig(atob(myargs["orig"].c_str()));
	sim.set_dynamic(myargs["dynamic"].c_str());
	model.dynamic_setTotalSim(myargs["dynamic"].c_str());
	
	cout << NN << "\n";

	sim.add_model(model);

	printf("\n----------------------\n");
	cout << ">Runscript:" << "\n";
	cout << YW << myscript << LCY << runScript << NN << "\n";
	cout << " - inputFile: " << LCY << myargs["fasta"].c_str() << NN << endl;
	cout << " - outFile  : " << LCY << myargs["outfile"] << NN << endl;
	printf("\n----------------------\n");

	if (atob(myargs["verbose"].c_str())) {
		cout << LCY << "#################\nDEBUG PRINT\n#################\n" << endl;
		cout << LGN << debugprint.str() << endl;
		cout << LGN << longarg.str() << endl;
		cout << LCY << "#################\nDEBUG PRINT DONE\n#################\n" << endl;
	}
	// exit(1);
	return 0;
}

void printUsage(std::string& myscript) {

	const char *LPR = "\033[1;95m";
	const char *LBU = "\033[1;94m";
	const char *LRD = "\033[1;31m";
	const char *LCY = "\033[1;36m";
	const char *YW = "\033[1;33m";
	const char *LGN = "\033[1;32m";
	const char *NN = "\033[1;0m";

	cout << endl << "Usage: " << YW << myscript << NN << " --i " << LCY << "<fastaFile.fa>" << NN << " --outdir " << LGN << "<outdir>" << " --o " << LGN << "<outprefix>" << NN << endl;
	
	cout << "\n----------------\n";
	cout << "Options\n" << endl;
	cout << "--i/--fasta     [file.fa]  : fasta input.fa." << LRD << " IMPORTANT: HEADER HAS TO BE IN THIS FORMAT (Example from 1-3908):\n\n" << LGN << endl;
	cout << ">mydef range=mygene:1-3908 5'pad=0 3'pad=0 strand=- repeatMasking=none" << endl;
	cout << "ACGTCGT..." << NN << endl;
	cout << "-o/--out        [string]   : output prefix name" << endl;
	cout << "--outdir        [string]   : output directory" << endl;
	cout << "--N             [int]      : specifies length of superhelical region in bp [1500]" << endl;
	cout << "--a [int]                  : junction energy in kcal/mol, for B DNA it's 10-11 [10]" << endl;
	cout << "--sigma         [float]    : superhelical density of the region [-0.07]" << endl;
	cout << "--minlength     [int]      : minimum rloop peak length in bp [50]" << endl;
	cout << "--nick          [int]      : specifies the ssDNA nick region [N/A]" << endl;
	cout << "--nicklen       [int]      : specifies the length of ssDNA nick region [N/A]" << endl;
	cout << "--beg           [int]      : get seq only from this bp from the input fasta [1]" << endl;
	cout << "--end           [int]      : get seq only until this bp from the input fasta [all]" << endl;
	cout << "--help/--h                 : print this message" << endl;
	cout << "--seed                     : Seed for RNG [42]" << endl;

	// cout << LGN;
	// cout << "--help/--h                 : print this message" << endl;
	// cout << "--H                        : print a LONG help (same as in README.md)" << endl;
	// cout << "--D/--dynamic              : Dynamic" << endl;
	// cout << "--naive                    : Equilibriuim model: naive method" << endl;
	// cout << "--orig                     : Equilibrium model: fixed original method" << endl;
	// cout << "--seed                     : Seed for RNG [42]" << endl;
	// cout << "--beg           [int]      : get seq only from this bp from the input fasta [1]" << endl;
	// cout << "--end           [int]      : get seq only until this bp from the input fasta [all]" << endl;
	
	// cout << NN;
	// cout << "\n----------------\n";
	// cout << "Biophysical Parameter\n" << endl;
	// cout << LGN;
	
	// cout << "--a [int]                  : junction energy in kcal/mol, for B DNA it's 10-11 [10]" << endl;
	
	// cout << "--A             [float]    : turns/bp [1/10.4]" << endl;
	// cout << "--C             [float]    : torsional stiffness for winding (1.8 for ssDNA or 3.6 for dsDNA) [1.8]" << endl;
	// cout << "--T             [float]    : temperature in Kelvin [310]" << endl;
	// cout << "--N             [int/auto] : size of superhelical domain in bp. 'auto' = automatically size to provided sequence. [1500]" << endl;
	// cout << "--sigma         [float]    : superhelical density of the region [-0.07]" << endl;
	// cout << "--txsigma       [float]    : superhelical density behind the RNA polymerase [-0.07] (Kouzine et al: -7%)" << endl;
	// cout << "--minlength     [int]      : minimum rloop peak length in bp [50]" << endl;
	// cout << "--w/--window    [int]      : window_size for dynamic rlooper [15]" << endl;
	// cout << "--unconstrained [F]        : toggle to remove the superhelicity term from the energy calculation [FALSE]" << endl;
	
	// cout << NN;
	// cout << "\n----------------\n";
	// cout << "Sequence Handling Overrides\n" << endl;
	// cout << LGN;
	
	// cout << "--homopolymer   [float]    : override input sequence basepairing energetics and use this as energy, e.g. 0.8 for homopolymer T [FALSE]" << endl;
	// cout << "--reverse       [F]        : toggle to reverse the sequence [FALSE]" << endl;
	// cout << "--revcomp       [F]        : toggle to complement the sequence [FALSE]" << endl;
	// cout << "--complement    [F]        : toggle to complement the sequence [FALSE]" << endl;
	// cout << "--reverse       [F]        : toggle to both reverse and complement the sequence [FALSE]" << endl;
	// cout << "--circular      [F]        : toggle to continously use the beginning sequence when calculation is at the end of the sequence" << endl;
	
	// cout << NN;
	// cout << "\n----------------\n";
	// cout << "Analysis Options\n" << endl;
	// cout << LGN;
	
	// cout << "--residuals     [F]        : toggle output residuals to STDOUT [FALSE]" << endl;
	// cout << "--dump          [F]        : toggle output every possible single R-loop structur to a file [FALSE]" << endl;
	// cout << "--top           [int]      : print out top N most favorable R-loop coordinates to STDOUT [0]" << endl;
	// cout << "--printG        [F]        : toggle print out the local average energy to " << YW << "<output_prefix>_avgG.wig" << LGN << " [FALSE]" << endl;
	// cout << "--printbed      [F]        : print output bed file <output_prefix>.bed [FALSE]" << endl;
	// cout << "--sensitivity   [int]      : --bedfile, peak calling power threshold [10]" << endl;
	
	// cout << NN;
	// cout << "\n----------------\n";
	// cout << "Misc options\n" << endl;
	// cout << LGN;
	// cout << "--import        [F]        : toggle to import peak file [FALSE]" << endl;
	// cout << "--sandbox       [F]        : development sandbox [FALSE]" << endl;
	// cout << NN;
	cout << "\n----------------\n\n\n";
}

void printUsageLong() {
	
	const char *LPR = "\033[1;95m";
	const char *LBU = "\033[1;94m";
	const char *LRD = "\033[1;31m";
	const char *LCY = "\033[1;36m";
	const char *YW = "\033[1;33m";
	const char *LGN = "\033[1;32m";
	const char *NN = "\033[1;0m";

	char buffer[FILENAME_MAX];
	if (getcwd(buffer, sizeof(buffer)) != nullptr) {
		std::string currentDirectory = buffer;
		std::string filePath = currentDirectory + "/src/readme_long.txt";
		ifstream f(filePath);
		if (f.is_open()) {
			cout << "\n\n" << LGN << f.rdbuf() << NN << endl;
		}
		else {
			cout << filePath << " here" << endl;
		}
	} else {
		cerr << "Error getting current directory" << endl;
	}
}

std::string turnSigmaIntoNiceStringForOutput(double sigma, stringstream& debugprint) { // example: -0.035 becomes min3p5 and 0.070 becomes pos7p0
	const char *LPR = "\033[1;95m";
	const char *LBU = "\033[1;94m";
	const char *LRD = "\033[1;31m";
	const char *LCY = "\033[1;36m";
	const char *YW = "\033[1;33m";
	const char *LGN = "\033[1;32m";
	const char *NN = "\033[0m";

	double number = sigma * 100;															//-0.035 -> -3.5
	int integerPart = static_cast<int>(number);											    // -3.5
	float floatfractionalPart = round(std::abs(number) * 10 - std::abs(integerPart) * 10);  // (3.5 - 3) * 10
	int fractionalPart = static_cast<int>(floatfractionalPart);
	std::ostringstream oss;
	if (integerPart < 0)
	{
		oss << "min" << std::abs(integerPart) << "p" << std::abs(fractionalPart);
	}
	else
	{
		oss << "pos" << std::abs(integerPart) << "p" << std::abs(fractionalPart);
	}
	
	debugprint << "rlooper.cpp::turnSigmaIntoNiceStringForOutput:" << NN << "\n" << LGN << 
	"- sigma = " << sigma << "\n" <<
	"- double number = sigma * 100 = " << sigma << "*" << 100 << " = " << number << "\n" <<
	"- int integerPart = int(number) = int(" << number << ") = " << integerPart << "\n" <<
	"- float floatfractionalPart = round(abs(number) * 10 - abs(integerPart * 10) = round(" << std::abs(number) << " * 10 - " << std::abs(integerPart) << " * 10) = " << floatfractionalPart << "\n" <<
	"- int fractionalPart = int(floatfractionalPart) = int(" << floatfractionalPart << ") = " << fractionalPart << "\n" << 
	"- finalIntegerPart = integerPart . p . fractionalPart = " << integerPart << " . p . " << fractionalPart << " = " << oss.str() << "\n" << 
	NN << "\n";
	
	return oss.str();
}

int countCharInString(String s, char delim){ // count number of char in string (like perl tr)
	int count = 0;
	String::size_type pos = s.find_first_of(delim);
	while ((pos = s.find_first_of(delim, pos)) != String::npos){
		count++;
		pos++;
	}
	return count;
}

std::string exec(const char* cmd) { // execute a command in terminal
    char buffer[128];
    std::string result = "";

    // Create a pipe to capture the command's output
    int pipe_fd[2];
    if (pipe(pipe_fd) == -1) {
        throw std::runtime_error("pipe() failed");
    }

    // Fork a child process
    pid_t pid = fork();
    if (pid == -1) {
        throw std::runtime_error("fork() failed");
    }

    if (pid == 0) { // Child process
        close(pipe_fd[0]); // Close the read end of the pipe

        // Redirect stdout to the write end of the pipe
        dup2(pipe_fd[1], STDOUT_FILENO);
        close(pipe_fd[1]);

        // Execute the command in the child process
        if (execl("/bin/sh", "sh", "-c", cmd, NULL) == -1) {
            throw std::runtime_error("execl() failed");
        }
    } else { // Parent process
        close(pipe_fd[1]); // Close the write end of the pipe

        // Read the command's output from the pipe
        ssize_t n;
        while ((n = read(pipe_fd[0], buffer, sizeof(buffer))) > 0) {
            result.append(buffer, n);
        }

        close(pipe_fd[0]); // Close the read end of the pipe

        // Wait for the child process to complete
        int status;
        waitpid(pid, &status, 0);

        if (status != 0) {
            throw std::runtime_error("Command execution failed");
        }
    }

    return result;
}

std::string chomp(string& text) { // remove trailing endline (perl chomp)
    if (!text.empty() && text[text.length() - 1] == '\n') {
        // Remove the trailing newline character
        text.erase(text.length() - 1);
    }
    return text;
}

void parseSettingsFile(const std::string& filename, 
					   std::map<int, std::string>& myorder, std::map<std::string, std::string>& myargs, 
					   std::stringstream& debugprint) { // parse settings.ini file to get default parametes
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    std::string line;
    int lineNumber = 0;

    // Read the file line by line
    while (std::getline(file, line)) {
        ++lineNumber;
        std::istringstream lineStream(line);
        chomp(line);
        std::string key, value;

        // Skip empty lines and lines that do not contain '='
        if (line.empty() || line.find('=') == std::string::npos) {
            continue;
        }

        // Extract key and value
        if (std::getline(lineStream, key, '=') && std::getline(lineStream, value)) {

			trimSpace(key);
			// key = "--" + key;
            // Add to maps
            myorder[lineNumber] = trimSpace(key);        // Store line number and key
            myargs[key]         = trimSpace(value);              // Store key and value
        }
    }

    // Print the maps to verify
    debugprint << "myorder (line_number -> key):" << std::endl;
    for (const auto& entry : myorder) {
        debugprint << entry.first << " -> " << entry.second << std::endl;
    }

    debugprint << "\nmyargs (key -> value):" << std::endl;
    for (const auto& entry : myargs) {
        debugprint << entry.first << " -> " << entry.second << std::endl;
    }

    file.close();
}

std::string getFilename(std::string& filePath) { // get Filename (c++ 11 so no std::filesystem)
    // Find the last slash in the file path
   	std::string filePath2 = filePath;
	filePath2 = getFullPath(filePath2);
    size_t lastSlashPos = filePath.find_last_of("/\\");
    
    // If no slash is found, the entire filePath is a filename
    if (lastSlashPos == std::string::npos) {
        return filePath;
    }
    
    // Extract the filename after the last slash
    return filePath.substr(lastSlashPos + 1);
}

std::string getFolder(std::string& filePath) { // c++ 11
    // Find the last slash in the file path
   	std::string filePath2 = filePath.c_str();
	filePath2 = getFullPath(filePath2);
    size_t lastSlashPos = filePath2.find_last_of("/\\");

    // If no slash is found, the entire filePath is a filename
    if (lastSlashPos == std::string::npos) {
        return "./";
    }
    // Extract the filename after the last slash
    return filePath2.substr(0,lastSlashPos);
}

void getFolderfull(std::string& filePath, std::string& folder, std::string& fileName) { // getFolder and getFilename
	size_t lastSlash = filePath.find_last_of("/\\");
	if (lastSlash != std::string::npos) {
		folder = filePath.substr(0, lastSlash);
		fileName = filePath.substr(lastSlash + 1);
	} else {
		folder = "";
		fileName = filePath;
	}
}

std::string getExtension(std::string& filePath) { // c++ 11
    // Find the last slash in the file path
   	std::string filePath2 = filePath.c_str();
	filePath2 = getFullPath(filePath2);
    size_t lastDotPos = filePath2.find_last_of(".");

    // If no slash is found, the entire filePath is a filename
    if (lastDotPos == std::string::npos) {
        return "";
    }
    // Extract the filename after the last slash
    return filePath2.substr(lastDotPos);
}

string whichFullPath(string& filePath) { // `which file` then get fullPath of it
	std::string mycmd = "realpath `which " + filePath + "`";
	return exec(mycmd.c_str());
}

std::string getFullPath(string& filePath) { // get fullPath of file
	std::string mycmd = "realpath " + filePath;
	return exec(mycmd.c_str());
}

std::string trimSpace(std::string& mystring) { // trim trailing spaces
	mystring.erase(0, mystring.find_first_not_of(" \t\n\r"));
	mystring.erase(mystring.find_last_not_of(" \t\n\r") + 1);
	return(mystring);
}


bool atob(const string& str) {
	string str2 = str;
    std::transform(str2.begin(), str2.end(), str2.begin(), ::tolower);
    std::istringstream is(str2);
    if (str == "False") {
        return(false);
    }
    else if (str == "True") {
        return(true);
    }
    else if (str == "") {
        return(false);
    }
    else {
        return(true);
    }
    bool b;
    is >> std::boolalpha >> b;
    return b;
}


inline bool fileexists (const std::string& file) {
  struct stat buffer;   
  return (stat (file.c_str(), &buffer) == 0); 
}

