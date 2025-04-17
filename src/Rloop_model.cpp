#include "Rloop_model.h"
#include <chrono>
#include "mytime.h"
#include <vector>
#include <cmath>

//
//
// Created by Robert Stolz on 6/27/17.
// Modified and maintanied by Stella Hartono 2023-2025

typedef std::chrono::high_resolution_clock::time_point TimeVar;


Rloop_model::Rloop_model() {
	nick = 0;
	nicklen = 1;
	selffoldlen=0;
	N = 1500; //1500bp is the experimentally determined length of the (-) sc domain after the transcription machinery
	A = 1/10.4; // turns/bp
	C = 1.8; //tortional stiffness of ssDNA winding. (Could be 3.6 for ds or 1.8 for ss winding)
	T = 310;
	k = (2200 * 0.0019858775 * T) / N; //Hooke's law coefficient: (2200*ideal_gas_constant in kcal/mol*absolute_temp_in_kelvin)/N
	a = 10; //Neucleation Free Energy in Kcals (~3-10.2kCals) 5000
	// -- sigma
	currSim = 0;
	iter = 0;
	sigma = -0.07; //measurement of energy upstream of replication domain
	alpha = N*sigma*A; //linking difference: topological parameter
	ambient_sigma = sigma;
	ambient_alpha = N * A * ambient_sigma;
	tx_ambient_sigma = -0.07; // transcriptional sigma
	tx_ambient_alpha = 0;
	tx_alpha = 0;
	tx_sigma = 0;
	sigma_total = ambient_sigma;
	alpha_total = ambient_alpha;


	//energy difference for the described RNA/DNA dinucleotide pairs (RNA/DNA - DNA/DNA) in kcal/mol as described in Huppert 2008
	//note, energies are for the RNA strand in the 5' to 3' direction
	rGG_dCC = -0.36;
	rGC_dCG = -0.16;
	rGA_dCT = -0.1;
	rGT_dCA = -0.06;
	rCG_dGC = 0.97;
	rCC_dGG = 0.34;
	rCA_dGT = 0.45;
	rCT_dGA = 0.38;
	rAG_dTC = -0.12;
	rAC_dTG = -0.16;
	rAA_dTT = 0.6;
	rAT_dTA = -0.12;
	rTG_dAC = .45;
	rTC_dAG = .5;
	rTA_dAT = .28;
	rTT_dAA = .8;

	rGG_rCC = -2.2;
	rGC_rCG = -2.4;
	rGA_rCT = -1.4;
	rGT_rCA = -1.5;
	rCG_rGC = -1.2;
	rCC_rGG = -1.5;
	rCA_rGT = -1.0;
	rCT_rGA = -0.9;
	rAG_rTC = -1.4;
	rAC_rTG = -1.6;
	rAA_rTT = -0.4;
	rAT_rTA = -1.0;
	rTG_rAC = -1.0;
	rTC_rAG = -0.8;
	rTA_rAT = -0.3;
	rTT_rAA = -0.2;

	homopolymer_override = false;
	unconstrained = false;
	override_energy = 0.0;
	
	//dynamic model specific parameters
	dynamic_window_size = 15;
	initiation_step_size = 1;
	elongation_step_size = 5;

}

int Rloop_model::getnick2() {
	return nick;
}

int Rloop_model::getnicklen2() {
	return nicklen;
}
int Rloop_model::getselffoldlen2() {
	return selffoldlen;
}

int Rloop_model::getnick() const {
	return nick;
}

int Rloop_model::getN() const {
	return N;
}

double Rloop_model::getA() const {
	return A;
}

double Rloop_model::getC() const {
	return C;
}

double Rloop_model::getK() const {
	return k;
}

double Rloop_model::geta() const {
	return a;
}

double Rloop_model::getSigma() const {
	return sigma;
}


double Rloop_model::getK2() {
	return k;
}

double Rloop_model::geta2() {
	return a;
}

double Rloop_model::getSigma2() {
	return sigma;
}

double Rloop_model::getAlpha() const {
	return alpha;
}

void Rloop_model::setN(int N) {
	Rloop_model::N = N;
	setAlpha(N*sigma*A);
	k = (2200 * 0.0019858775 * T) / N;
}

void Rloop_model::setnick(int n) {
	Rloop_model::nick = n;
}

void Rloop_model::setnicklen(int n) {
	Rloop_model::nicklen = n;
}

void Rloop_model::setselffoldlen(int n) {
	Rloop_model::selffoldlen = n;
}

void Rloop_model::setA(double A) {
	Rloop_model::A = A;
	setAlpha(N*sigma*A);
}

void Rloop_model::setC(double C) {
	Rloop_model::C = C;
}

void Rloop_model::setK(double k) {
	Rloop_model::k = k;
}

void Rloop_model::seta(double a) {
	Rloop_model::a = a;
}

void Rloop_model::set_superhelicity(double sigma) {
	Rloop_model::sigma = sigma;
	setAlpha(N*sigma*A);
}

void Rloop_model::set_unconstrained(bool value) {
	unconstrained = value;
}

void Rloop_model::setAlpha(double alpha) {
	Rloop_model::alpha = alpha;
}

double Rloop_model::getT() const {
	return T;
}

void Rloop_model::setT(double T) {
	Rloop_model::T = T;
	k = (2200 * 0.0019858775 * T) / N;
}
//member functions


int Rloop_model::dynamic_getIter() const {
	return iter;
}

void Rloop_model::dynamic_setIter(int value) {
	Rloop_model::iter = value;
}

double Rloop_model::dynamic_getSigma() const {
	return sigma;
}

void Rloop_model::dynamic_setSigma(double value) {
	Rloop_model::sigma = value;
}

double Rloop_model::dynamic_getSigmaTotal() const {
	return sigma_total;
}

void Rloop_model::dynamic_setSigmaTotal(double value) {
	Rloop_model::sigma_total = value;
}

double Rloop_model::dynamic_getAlpha() const {
	return alpha;
}

void Rloop_model::dynamic_setAlpha(double value) {
	Rloop_model::alpha = value;
}

double Rloop_model::dynamic_getAlphaTotal() const {
	return alpha_total;
}

void Rloop_model::dynamic_setAlphaTotal(double value) {
	Rloop_model::alpha_total = value;
}

int Rloop_model::dynamic_getCurrentPos() const {
	return current_pos;
}

int Rloop_model::dynamic_getCurrSim() const {
	return currSim;
}

void Rloop_model::dynamic_setCurrSim(int value){
	currSim = value;
}

int Rloop_model::dynamic_getTotalSim() const {
	return totalSim;
}

void Rloop_model::dynamic_setTotalSim(string value){
	if (value == "False") {
		totalSim = 0;
	}
	else {
		totalSim = atoi(value.c_str());
	}
}

void Rloop_model::dynamic_reset_model() {
	N = 1500; //1500bp is the experimentally determined length of the (-) sc domain after the transcription machinery
	A = 1/10.4; // turns/bp
	C = 1.8; //tortional stiffness of ssDNA winding. (Could be 3.6 for ds or 1.8 for ss winding)
	T = 310;
	k = (2200 * 0.0019858775 * T) / N; //Hooke's law coefficient: (2200*ideal_gas_constant in kcal/mol*absolute_temp_in_kelvin)/N
	a = 10; //Neucleation Free Energy in Kcals (~3-10.2kCals) 5000
	// -- sigma
	currSim = 0;
	iter = 0;
	//sigma = -0.07; //measurement of energy upstream of replication domain
	alpha = N*sigma*A; //linking difference: topological parameter
	ambient_sigma = sigma;
	ambient_alpha = N * A * ambient_sigma;
	tx_ambient_sigma = -0.07; // transcriptional sigma
	tx_ambient_alpha = current_pos * A * tx_ambient_sigma; // 0 
	tx_alpha = 0;
	tx_sigma = 0;
	sigma_total = ambient_sigma + tx_sigma; // --sigma + 0
	alpha_total = ambient_alpha + tx_alpha; // N*A*--sigma + 0

	in_rloop = false;
	current_pos = 0;
	current_bp_energy = 0;
	current_junction_energy = 0;
	current_superhelical_energy = 0;
	proposed_bp_energy = 0;
	proposed_junction_energy = 0;
	proposed_superhelical_energy = 0;
	partition_function = 0;
	n_rloops = 0;
	n_rloop_bases = 0; //dynamic_window_size;
	current_pos = dynamic_window_size-1;
	//ambient_alpha = (1+n_rloop_bases)*A*sigma; // alpha
	// double current_alpha_total = ambient_alpha +
		// (dynamic_getCurrentPos()*tx_ambient_sigma*getA());
	
	// Gsigma with new alpha total, but why 4*instead of 2?
	// double current_Gsigma = (
		// (4*pow(pi,2)*C) * (current_alpha_total+n_rloop_bases*A) /
		// (4*pow(pi,2)*C+k*n_rloop_bases)
		// ); 
}

void Rloop_model::dynamic_print_topological_state(std::stringstream &ss_dynamic_log){}

void Rloop_model::dynamic_step_forward_initiation(vector<vector<double>> &df_bfs, vector<vector<double>> &df_bftotals, vector<vector<double>> &df_probs, bool &verbose, int &seed, int &begpos, double &myurn, double &myprob, std::vector<std::string> &myvalues) //, std::stringstream &ss_dynamic_log)
{

	const char *LPR = "\033[95m";
	const char *LBU = "\033[94m";
	const char *LRD = "\033[31m";
	const char *LCY = "\033[36m";
	const char *YW = "\033[33m";
	const char *LGN = "\033[32m";
	const char *NN = "\033[0m";
	TimeVar t1 = timeNow();

	proposed_superhelical_energy = 0;
	proposed_junction_energy = 0;
	partition_function = 0;
	long double Gbps = 0;
	long double Gsigma = 0;
	long double Gtot = 0;
	long double Gjunc = 0;
	long double Gtot_bf = 0;
	proposed_bp_energy = 0;
	current_bp_energy = 0;

	for (int i=current_pos-dynamic_window_size+1; i < current_pos; i++){
		int m=0;
		int my_n = i;
		int my_m = 0;
		proposed_bp_energy = 0;

		for (int j=i+1; j <= current_pos; j++) { 
		
			char nuc1 = seqD[j - 1];
			char nuc2 = seqD[j];
			m = j-i+1+n_rloop_bases;
			my_m = m;
			if (m != 0)
			{
				my_m = m - 1;
			}

			double Gbp = compute_bps_interval(seqD[j-1], seqD[j]);
			proposed_bp_energy += Gbp;
			proposed_superhelical_energy = (2 * pow(pi, 2) * C * k * pow((alpha_total + m * A), 2)) / (4 * pow(pi, 2) * C + k * m);
			proposed_junction_energy = current_junction_energy + .5 * a;
			Gbps += Gbp;
			Gsigma += proposed_superhelical_energy;
			Gjunc += proposed_junction_energy;
			Gtot += proposed_bp_energy + proposed_superhelical_energy + proposed_junction_energy;
			Gbps += Gbp;
			Gsigma += proposed_superhelical_energy;
			Gjunc += proposed_junction_energy;
			Gtot += proposed_bp_energy + proposed_superhelical_energy + proposed_junction_energy;
			Gtot_bf += compute_boltzmann_factor(proposed_bp_energy + proposed_superhelical_energy + proposed_junction_energy,T);
			partition_function += compute_boltzmann_factor(proposed_bp_energy + proposed_superhelical_energy + proposed_junction_energy,T);
		}
	}
	long double G0 = 0.5 * k * pow(alpha_total,2);
	long double G0_bf = compute_boltzmann_factor(.5 * k * pow(alpha_total,2),T);
	long double Gtot_bf0 = partition_function;
	partition_function += G0_bf;
	double myrand = (double)rand();
	
	myrand = (double)rand();
	
	double myrandmax = (double)RAND_MAX;
	double RNG = myrand / myrandmax;
	long double prob_all_rloop_init = 1-(G0_bf/partition_function); 
	myurn = RNG;
	myprob = prob_all_rloop_init;
	stringstream Gjunc_str;
	Gjunc_str << Gjunc;
	stringstream Gbps_str;
	Gbps_str << Gbps;
	stringstream Gsigma_str;
	Gsigma_str << Gsigma;
	stringstream Gtot_str;
	Gtot_str << Gtot;
	stringstream G0_str;
	G0_str << G0;
	stringstream Gtot_bf0_str;
	Gtot_bf0_str << Gtot_bf0;
	stringstream RNG_str;
	RNG_str << RNG;
	stringstream myprob_str;
	myprob_str << myprob;
	stringstream Gtot_bf_str;
	Gtot_bf_str << Gtot_bf;
	stringstream G0_bf_str;
	G0_bf_str << G0_bf;
	stringstream pf;
	pf << partition_function;
	cout << "init rand=" << myrand << ", randmax=" << myrandmax << ", RNG=" << RNG << ", prob=" << prob_all_rloop_init << "\n";
	if (RNG >= prob_all_rloop_init) {
		vector<string> valuesToAdd = {"0_none",Gjunc_str.str(),Gbps_str.str(),Gsigma_str.str(),Gtot_str.str(),G0_str.str(),Gtot_bf_str.str(),G0_bf_str.str(),pf.str(),RNG_str.str(),myprob_str.str()};
		myvalues.insert(myvalues.end(), valuesToAdd.begin(), valuesToAdd.end());

		current_pos += initiation_step_size;
	}
	else {
		vector<string> valuesToAdd = {"1_init",Gjunc_str.str(),Gbps_str.str(),Gsigma_str.str(),Gtot_str.str(),G0_str.str(),Gtot_bf_str.str(),G0_bf_str.str(),pf.str(),RNG_str.str(),myprob_str.str()};
		myvalues.insert(myvalues.end(), valuesToAdd.begin(), valuesToAdd.end());

		in_rloop = true;

		n_rloops++;
		n_rloop_bases += dynamic_window_size;
		for (int i=current_pos-dynamic_window_size+2; i < current_pos; i++){
			current_bp_energy += compute_bps_interval(seqD[i], seqD[i+1]);
		}
		current_junction_energy = .5 * a;
		int m = dynamic_window_size;
		current_superhelical_energy =
			(2 * pow(pi, 2) * C * k * pow((alpha + m * A), 2)) / (4 * pow(pi, 2) * C + k * m);
		begpos = current_pos + 2 - dynamic_window_size;
		current_loci.start_pos = current_pos + 2 - dynamic_window_size;
	}
}

bool Rloop_model::dynamic_step_forward_elongation(vector<vector<double>> &df_bfs, vector<vector<double>> &df_bftotals, vector<vector<double>> &df_probs, bool &verbose, int &seed, int &endpos, double &myurn, double &myprob, vector<string> &myvalues) 
{
	long double prob_all_rloop_init = myprob;
	const char *LPR = "\033[95m";
	const char *LBU = "\033[94m";
	const char *LRD = "\033[31m";
	const char *LCY = "\033[36m";
	const char *YW = "\033[33m";
	const char *LGN = "\033[32m";
	const char *NN = "\033[0m";

	long double Gjunc = 0;
	long double Gbps = current_bp_energy;
	long double Gsigma = 0;
	long double Gtot = 0;
	long double Gtot_bf = 0;

	proposed_bp_energy = current_bp_energy;
	current_pos += elongation_step_size;
	if (current_pos >= seqD.size()-1){ 
		return false;
	}
	for (int i=current_pos-elongation_step_size; i < current_pos; i++){
		double Gbp = compute_bps_interval(seqD[i],seqD[i+1]);
		proposed_bp_energy += compute_bps_interval(seqD[i],seqD[i+1]);
		Gbps += Gbp;
	}

	int m = n_rloop_bases + elongation_step_size;
	proposed_superhelical_energy =
		(2 * pow(pi, 2) * C * k * pow((alpha_total + m * A), 2)) / (4 * pow(pi, 2) * C + k * m); 
	Gsigma += proposed_superhelical_energy;

	current_superhelical_energy =
		(2 * pow(pi, 2) * C * k * pow((alpha_total + n_rloop_bases * A), 2)) / (4 * pow(pi, 2) * C + k * n_rloop_bases);
	proposed_junction_energy = current_junction_energy+0.5*a; 
	
	Gjunc = current_junction_energy; 

	long double G0 = 0.5 * k * pow(alpha_total,2);
	long double G0_bf = compute_boltzmann_factor(.5 * k * pow(alpha_total,2),T);
	Gtot = Gbps + Gsigma + Gjunc; 
	Gtot_bf = compute_boltzmann_factor(Gtot,T);
	long double proposed_total_energy_bf = compute_boltzmann_factor(proposed_bp_energy+proposed_superhelical_energy+current_junction_energy,T);
	partition_function += proposed_total_energy_bf;
	double myrand = (double)rand();
	myrand = (double)rand();
	double myrandmax = (double)RAND_MAX;
	double RNG = myrand / myrandmax * prob_all_rloop_init;
	long double prob_all_rloop_elon = (proposed_total_energy_bf/partition_function); 
	myurn = RNG;
	myprob = prob_all_rloop_elon;
	stringstream Gjunc_str;
	Gjunc_str << Gjunc;
	stringstream Gbps_str;
	Gbps_str << Gbps;
	stringstream Gsigma_str;
	Gsigma_str << Gsigma;
	stringstream Gtot_str;
	Gtot_str << Gtot;
	stringstream G0_str;
	G0_str << G0;
	stringstream RNG_str;
	RNG_str << RNG;
	stringstream myprob_str;
	myprob_str << myprob;
	stringstream Gtot_bf_str;
	Gtot_bf_str << Gtot_bf;
	stringstream G0_bf_str;
	G0_bf_str << G0_bf;
	stringstream pf;
	pf << partition_function;

	cout << "elon rand=" << myrand << ", randmax=" << myrandmax << ", RNG=" << RNG << ", prob=" << prob_all_rloop_elon << "\n";
	if (RNG < prob_all_rloop_elon){ 
		vector<string> valuesToAdd = {"2_elon",Gjunc_str.str(),Gbps_str.str(),Gsigma_str.str(),Gtot_str.str(),G0_str.str(),Gtot_bf_str.str(),G0_bf_str.str(),pf.str(),RNG_str.str(),myprob_str.str()};
		myvalues.insert(myvalues.end(), valuesToAdd.begin(), valuesToAdd.end());
		n_rloop_bases += elongation_step_size;
		current_superhelical_energy = proposed_superhelical_energy;
		current_bp_energy = proposed_bp_energy;
		}
	else{ 
		in_rloop = false;
		current_junction_energy = a;
		current_loci.end_pos = current_pos+1;
		endpos = current_pos + 1;
		rloop_structures.push_back(current_loci);
		vector<string> valuesToAdd = {"3_end",Gjunc_str.str(),Gbps_str.str(),Gsigma_str.str(),Gtot_str.str(),G0_str.str(),Gtot_bf_str.str(),G0_bf_str.str(),pf.str(),RNG_str.str(),myprob_str.str()};
		myvalues.insert(myvalues.end(), valuesToAdd.begin(), valuesToAdd.end());
		current_pos = current_pos + dynamic_window_size; 
		return false;
	}
	return true;
}

string Rloop_model::join(char &delim, const vector<char> &lst) {
	string join_this_vector_char;
	int i = 0;
	for(const auto &s : lst) {
		if (i == 0) {
			join_this_vector_char += delim;
			i = 1;
		}
		join_this_vector_char += s;
	}
	return join_this_vector_char;
}

string Rloop_model::join(const vector<char> &lst) {
	string join_this_vector_char;
	for(const auto &s : lst) {
		join_this_vector_char += s;
	}
	return join_this_vector_char;
}

void Rloop_model::step_forward(vector<char> &sequence, std::vector<char>::iterator &b_0, std::vector<char>::iterator &b_1) {
	if (b_0 == sequence.end()) { 
		b_0 = sequence.begin();   
		b_1 = b_0 + 1;            
	}
	else if (b_0 == sequence.end() - 1) {
		b_1 = sequence.begin();
	}
	else{
		b_1 = b_0 + 1;
	}
}

int Rloop_model::find_distance(vector<char>& sequence, const std::vector<char>::iterator &start, const std::vector<char>::iterator &stop, Structure& structure) {
	int n=1;
	if (stop < start) {
		n += sequence.size() - distance(stop,start);
	}
	else {
		n += distance(start,stop);
	}
	return n;
}

void set_bp_energy_override();

/* Get basepairing energy by dinucleotide. 
 * Changed from original rlooper's, which complements the sequence, then complement the check dinucleotides.
 * e.g. previously, if seq = ACGT -> complement becomes TGCA. Checks below are 'first == T' second == 'G' return rAC_dTG
 * Instead just use the top/non-template sequence like usual. Checks below are 'first == A' second == 'C' return rAC_dTG
 */
double Rloop_model::compute_bps_interval(const char &first, const char &second) {
	if (homopolymer_override == true) {
		return override_energy;
	}
	
	if (first == 'G') { //G
		if (second == 'G') //GG
			return rGG_dCC;
		else if (second == 'C') //GC
			return rGC_dCG;
		else if (second == 'A') //GA
			return rGA_dCT;
		else if (second == 'T'|| second == 'U')
			return rGT_dCA;
	}
	else if (first == 'C') { //G
		if (second == 'G') //GC
			return rCG_dGC;
		else if (second == 'C') //GG
			return rCC_dGG;
		else if (second == 'A') //GT
			return rCA_dGT;
		else if (second == 'T'|| second == 'U')
			return rCT_dGA;
	}
	else if (first == 'A') { //A
		if (second == 'G') //AC
			return rAG_dTC;
		else if (second == 'C') //AG
			return rAC_dTG;
		else if (second == 'A') //AA
			return rAA_dTT;
		else if (second == 'T'|| second == 'U')
			return rAT_dTA;
	}
	else if (first == 'T') { //T
		if (second == 'G') //TG
			return rTG_dAC;
		else if (second == 'C') //TC
			return rTC_dAG;
		else if (second == 'A') //TA
			return rTA_dAT;
		else if (second == 'T'|| second == 'U')
			return rTT_dAA;
	}
	
	// if nothing return at this point, then die coz there's unexpected sequence that's non-A/C/G/T/U
	throw std::invalid_argument("There is a non-A/C/G/T/U in the sequence!\n");
	exit(1);
}

void Rloop_model::set_bp_energy_override(double energy) {
	homopolymer_override = true;
	override_energy = energy;
}

vector<char> Rloop_model::copy(const vector<char> &sequence) {
	vector<char> sequence_copy;
	sequence_copy.insert(sequence_copy.begin(), sequence.begin(), sequence.end());
	return(sequence_copy);
}

void Rloop_model::complement_sequence(vector<char> &sequence) {
	for (int i = 0; i < sequence.size(); i++) {
		if (sequence[i] == 'A')
			sequence[i] = 'T';
		else if (sequence[i] == 'T')
			sequence[i] = 'A';
		else if (sequence[i] == 'C')
			sequence[i] = 'G';
		else //sequence_data[i] == 'G'
			sequence[i] = 'C';
	}
}

void Rloop_model::reverse_sequence(vector<char> &sequence) {
	char temp;
	for (int i = 0; i < (sequence.size()/2); i++) {
		temp = sequence[i];
		sequence[i] = sequence[sequence.size() - 1 - i];
		sequence[sequence.size() - 1 - i] = temp;
	}
}

// intentionally put here, I know it's the same as above and just re-writing
void Rloop_model::reverse_complement_sequence(vector<char> &sequence) {
	complement_sequence(sequence);
	reverse_sequence(sequence);
}

void Rloop_model::compute_structure(vector<char> &sequence, 
                                    const std::vector<char>::iterator &start, 
                                    const std::vector<char>::iterator &stop,
                                    Structure& structure,
                                    float &bp_energy,
                                    float &time1,
                                    float &time2,
                                    float &time3,
                                    bool &verbose) {
	const char* LRD = "\033[31m";
	const char* LCY = "\033[36m";
	const char* YW = "\033[33m";
	const char* LGN = "\033[32m";
	const char* NN  = "\033[0m";
	bool myverbose = true; 
	typedef std::chrono::high_resolution_clock::time_point TimeVar; 
	
	// TIME starts
	TimeVar t1 = timeNow();
	
	long int m = find_distance(sequence,start,stop,structure); 
	long int n = structure.position.start_pos;
	
	time1 += duration(timeNow()-t1)/1e9;
	
	if (!unconstrained) {
		structure.gsigma = (2 * pow(pi, 2) * C * k * pow((alpha + m * A), 2)) / (4 * pow(pi, 2) * C + k * m);
	}
	structure.free_energy += structure.gsigma; 
	
	time2 += duration(timeNow()-t1)/1e9;
	
	vector<char>::iterator b_0 = start + m - 2;
	vector<char>::iterator b_1 = b_0 + 1;
	
	step_forward(sequence, b_0, b_1);
	
	bp_energy += compute_bps_interval(*b_0, *b_1);
	if (bp_energy < 0 && -1 * bp_energy < 5e-9) {bp_energy = 0;} 
	if (bp_energy > 0 &&  1 * bp_energy < 5e-9) {bp_energy = 0;} 
	structure.bp_energy = bp_energy;
	structure.free_energy += structure.bp_energy;
	
	time3 += duration(timeNow()-t1)/1e9;
	
	
	structure.free_energy += a;
	
	structure.boltzmann_factor = compute_boltzmann_factor(structure.free_energy,T);
}

void Rloop_model::compute_external_structure(Structure& structure, Structure& rloop, Peak& external, bool &verbose) {
	std::vector<char>::iterator b_0,b_1;
	long int m = rloop.position.get_length();
	if (!unconstrained) {
		structure.free_energy = (2 * pow(pi, 2) * C * k * pow((alpha + m * A), 2)) / (4 * pow(pi, 2) * C + k * (m-external.position.get_length())) + external.intensity;
	}
	structure.free_energy += rloop.bp_energy;
	structure.boltzmann_factor = compute_boltzmann_factor(structure.free_energy,T);
	
}

void Rloop_model::compute_residuals(Structure &structure, bool &verbose) {
	structure.current_Gsigma = ((4*pow(pi,2)*C) / (4*pow(pi,2)*C+k*structure.position.get_length()-structure.external_length)) * (alpha+structure.position.get_length()*A);
	structure.residual_twist = ((2*pi*k) / (4*pow(pi,2)*C+k*structure.position.get_length()-structure.external_length)) * (alpha+structure.position.get_length()*A);
}

void Rloop_model::ground_state_residuals(double &twist, double &writhe) {
	writhe = alpha;
	twist = k/(2*pi)*C*alpha;
}

long double Rloop_model::ground_state_factor() {
	return compute_boltzmann_factor(((k*pow(alpha, 2)) / 2),T); 
}

long double Rloop_model::ground_state_energy() {
	return ((k*pow(alpha, 2)) / 2); 
}

double Rloop_model::dynamic_compute_res_lk(){
	return ((4*pow(pi,2)*C) / (4*pow(pi,2)*C+k*n_rloop_bases)) * (alpha_total+n_rloop_bases*A);
}

int Rloop_model::dynamic_getWindow_size() const {
	return dynamic_window_size;
}

void Rloop_model::dynamic_setWindow_size(int Dynamic_window_size) {
	Rloop_model::dynamic_window_size = Dynamic_window_size;
}

int Rloop_model::dynamic_getN_rloops() const {
	return n_rloops;
}

void Rloop_model::dynamic_setN_rloops(int n_rloops) {
	Rloop_model::n_rloops = n_rloops;
}

int Rloop_model::dynamic_getInitiation_step_size() const {
	return initiation_step_size;
}

void Rloop_model::dynamic_setInitiation_step_size(int initiation_step_size) {
	Rloop_model::initiation_step_size = initiation_step_size;
}

int Rloop_model::dynamic_getElongation_step_size() const {
	return elongation_step_size;
}

void Rloop_model::dynamic_setElongation_step_size(int elongation_step_size) {
	Rloop_model::elongation_step_size = elongation_step_size;
}

double Rloop_model::dynamic_getTx_Amb_Sigma() const {
	return tx_ambient_sigma;
}

void Rloop_model::dynamic_setTx_Amb_Sigma(double value) {
	Rloop_model::tx_ambient_sigma = value;
}

double Rloop_model::dynamic_getTx_Amb_Alpha() const {
	return tx_ambient_alpha;
}

void Rloop_model::dynamic_setTx_Amb_Alpha(double value) {
	Rloop_model::tx_ambient_alpha = value;
}

void Rloop_model::dynamic_setTx_Sigma(double value) {
	Rloop_model::tx_sigma = value;
}

double Rloop_model::dynamic_getTx_Sigma() const {
	return tx_sigma;
}

void Rloop_model::dynamic_setTx_Alpha(double value) {
	Rloop_model::tx_alpha = value;
}

double Rloop_model::dynamic_getTx_Alpha() const {
	return tx_alpha;
}

int Rloop_model::dynamic_getN_rloop_bases() const {
	return n_rloop_bases;
}

void Rloop_model::dynamic_setN_rloop_bases(int n_rloop_bases) {
	Rloop_model::n_rloop_bases = n_rloop_bases;
}
