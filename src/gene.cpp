//
// Created by Robert Stolz on 6/27/17.
// Modified & maintained by Stella Hartono 2023-2024
//

#include "gene.h"
#include <chrono>
#include "mytime.h"
//constructors
using namespace std;


void Gene::parse_header(){
	unsigned long pos = 0;
	string name, remaining;
	// cout << "name orig : " << header << "\n";
	pos = header.find('=');
	name = header.substr(1,pos);
	// cout << "name =    : " << name << "\n";
	pos = name.find(' ');
	name = name.substr(0,pos);
	// cout << "name space: " << name << "\n";
	gene_name = name;
	position.chromosome = name;
	std::transform(position.chromosome.begin(), position.chromosome.end(), position.chromosome.begin(), ::tolower);
	std::transform(position.chromosome.end()-1, position.chromosome.end(), position.chromosome.end()-1, ::toupper);
	position.start_pos = 1;
	position.end_pos = sequence.size();
	position.strand = "+";
	return;
}
//constructors and destructors
Gene::Gene(){

}

Gene::~Gene(){
	clear_structures();
}

//getters and setters
string Gene::getName(){
	return gene_name;
}

const string &Gene::getHeader() const {
	return header;
}

void Gene::setHeader(const string &header) {
	Gene::header = header;
}

const Loci &Gene::getPosition() const {
	return position;
}

void Gene::setPosition(const Loci &position) {
	Gene::position = position;
}

const vector<char, allocator<char>> &Gene::getSequence() const {
	return sequence;
}

vector<Structure>& Gene::getRloopStructures(){
	return rloop_structures;
}

float Gene::compute_GC_skew(){
	
	if (sequence.size()){
		throw EmptyGeneException();
	}
	float Gs = 0.f;
	float Cs = 0.f;
	for (int i=0; i < sequence.size(); i++){
		if (sequence[i] == 'G'){
			Gs += 1.f;
		}
		else if (sequence[i] == 'C'){
			Cs += 1.f;
		}
	}
	return (Gs-Cs)/(Gs+Cs);
}

bool Gene::read_gene(ifstream& fastafile, int genebeg, int geneend) { //need to test

	const char *LPR = "\033[1;95m";
	const char *LBU = "\033[1;94m";
	const char *LRD = "\033[1;31m";
	const char *LCY = "\033[1;36m";
	const char *YW = "\033[1;33m";
	const char *LGN = "\033[1;32m";
	const char *NN = "\033[1;0m";

	std::vector<char> origseq;
	//initialize variables
	char c,p;
	int currpos = 0;
	//check that fstream is open and not eof, throw exceptions otherwise
	if (!fastafile.is_open()) {
		throw UnexpectedClosedFileException("Gene::read_gene");
	} else if (fastafile.eof()) {
		throw UnexpectedEOFException();
	}
	while (fastafile.get(c)) {
		//read the next character
		c = toupper(c);
		//if the character is the start of a header line
		if (c == '>') {
			//read until the end of the header line
			while (c != '\n' && c != '\r') {
				header.push_back(c);
				c = toupper(fastafile.get());
			}
			parse_header();
		}
		else if (c == 'A' || c == 'T' || c == 'C' || c == 'G' || c == 'N') {
			currpos = currpos + 1;
			//save to sequence vector
			origseq.push_back(c);
			if (genebeg != -1 && currpos <= genebeg) {
				continue;
			}
			if (geneend != -1 && currpos > geneend) {
				if (currpos > 20000) {
					break;
				}
				continue;
			}
			sequence.push_back(c);
		}
		//else if encountered some white space
		else if (c == '\n' || c == ' ' || c == '\t' || c== '\r'){
			p = fastafile.peek();
			if (p == '>'){
				break;
				return false;
			}
			/*else if (p == EOF) {
			 return true;
			}*/
			continue;
		}
		//else the charicter is unrecognized
		else {
			throw InvalidSequenceDataException(c);
		}
	}

	//cout << "gene.cpp: 1. position.start_pos     - position.end_pos (original): " << position.start_pos << "-" << position.end_pos << "\n";


	if (genebeg == -1) {
		fabeg = 0;
	}
	else {
		fabeg = genebeg;
	}
	if (geneend == -1) {
		faend = sequence.size();
	}
	else {
		faend = geneend;
	}

	position.start_pos = 1;
	position.end_pos = sequence.size();// - 1;
	// cout << "gene.cpp: position.start_pos (1) - position.end_pos (sequence.size()): " << position.start_pos << "-" << position.end_pos << "\n";
	// cout << string(sequence.begin(),sequence.end()) << "\n\n";
	windower.set_sequence(sequence);

	// cout << genebeg << "-" << geneend << "\n";
	int seqlen = origseq.end() - origseq.begin();
	string seqshort = string(origseq.begin(), origseq.end());
	if (seqlen > 20) {
		seqshort = string(origseq.begin(), origseq.begin() + 10) + "..." + string(origseq.end() - 10, origseq.end());
	}

	seqlen = sequence.end() - sequence.begin();
	seqshort = string(sequence.begin(), sequence.end());
	if (seqlen > 20) {
		seqshort = string(sequence.begin(), sequence.begin() + 10) + "..." + string(sequence.end() - 10, sequence.end());
	}

	cout << LPR << "gene.cpp::read_gene" << NN << ": original seq (" << LGN << 0 << "-" << int(origseq.size()) << NN << ", " << LGN << seqlen << NN << " bp): " << LCY << seqshort << NN << "\n";
	cout << LPR << "gene.cpp::read_gene" << NN << ":     used seq (" << LGN << fabeg << "-" << faend << NN << ", " << LGN << seqlen << NN << " bp): " << LCY << seqshort << NN << "\n";
	
	cout << LPR << "gene.cpp::read_gene:" << NN << ":" << LGN << "SUCCESS!!" << NN << "\n\n";
	return true;
	//unexpected EOF
	//throw InvalidSequenceDataException();
}


void Gene::print_gene(){
	std::cout << header << std::endl;
	for (std::vector<char>::iterator it = sequence.begin(); it<sequence.end(); ++it){
		cout << *it;
	}
	std::cout << std::endl;
}

void Gene::compute_structures(Model &model, bool &verbose){
	// windower.set_sequence(test);
	//iterate through sequence
	typedef std::chrono::high_resolution_clock::time_point TimeVar;
	const char* LRD = "\033[31m";
	const char* LCY = "\033[36m";
	const char* YW = "\033[33m";
	const char* LGN = "\033[32m";
	const char* NN  = "\033[0m";
	float time1 = 0;
	float time2 = 0;
	float time3 = 0;
	// duration(timeNow()-t1)/1e9
	//check for circular sequence conditions
	if (sequence.size() == 0){
		//throw exception
	}
	//initializing the iterators ensures that the intial comparison in next_window_from_all_windows is not problematic
	std::vector<char>::iterator start = sequence.begin(),stop=sequence.begin()+1;
	windower.reset_window();
	int last_pos = -1;
	float bp_energy = 0;
	int iter0 = 0;
	while (windower.has_next_window()){
		iter0 ++;
		
		windower.next_window_from_all_windows(start,stop);
		Structure temp;
		//set the Loci of the structure using the gene's Loci
		temp.position.chromosome = position.chromosome;
		temp.position.strand     = position.strand;
		temp.position.start_pos  = position.start_pos + windower.get_current_start_offset();
		temp.position.end_pos    = position.start_pos + windower.get_current_stop_offset();
		if (last_pos != temp.position.start_pos) {
			bp_energy = 0;
		}
		last_pos = temp.position.start_pos;
		//pass the structure and window boundaries to the model
		model.compute_structure(sequence,start,stop,temp,bp_energy,time1,time2,time3,verbose);
		//push the now computed structure onto these_structures
		rloop_structures.push_back(temp); //need to make sure the default copy constructor is working properly
		// if (verbose) {
		// 	cout << iter0 << "," << temp.position.start_pos << "," << temp.position.end_pos << "\n";
		// }
	}
	if (verbose) {
		// cout << "Total index   = " << iter0 << "\n";
		// cout << "Sequence size = " << sequence.size() << "\n";
		cout << YW << "- rloop_model.cpp:: TOTAL O(n*n) calculate long int m: " << LGN << time1 << NN << "\n";
		cout << YW << "- rloop_model.cpp:: TOTAL O(n*1) Gsigma: " << LGN << time2 << NN << "\n";
		cout << YW << "- rloop_model.cpp:: TOTAL O(n*m) calculate Gbp: " << LGN << time3 << NN << "\n";
	}
	ground_state_energy = model.ground_state_energy();
}

void Gene::compute_structures_external(vector<Peak> &external_structures, Model &model, bool &verbose){
	//for each external structure
	for (int i=0; i < external_structures.size(); i++) {
		//reset windower position on the sequence
		//re-coordinate this external structure to match R-looper's internal structures
		if (position.strand == "-") { //not thoroughly tested
			external_structures[i].position.start_pos = get_length() - external_structures[i].position.start_pos;
			external_structures[i].position.end_pos = get_length() - external_structures[i].position.end_pos;
		}
		//for every R-loop structure
		int n_rloop_structures = rloop_structures.size();
		for (int j=0; j < n_rloop_structures; j++){
			//cout << j << endl;
			if (external_structures[i].position.start_pos > rloop_structures[j].position.start_pos &&
       external_structures[i].position.end_pos < rloop_structures[j].position.end_pos){
				//cout << "hit" << endl;
				Structure temp;
				//set the Loci of the structure using the gene's Loci
				temp.position.chromosome = position.chromosome;
				temp.position.strand = position.strand;
				temp.position.start_pos = rloop_structures[j].position.start_pos;
				temp.position.end_pos = rloop_structures[j].position.end_pos;
				temp.external = true;
				temp.external_length = external_structures[i].position.get_length();
				//pass the structure and window boundaries to the model
				model.compute_external_structure(temp, rloop_structures[j], external_structures[i],verbose);
				//push the now computed structure onto these_structures
				rloop_structures.push_back(temp); //need to make sure the default copy constructor is working properly
			}
		}
	}
}

vector<Structure> Gene::compute_structures_dynamic(Model& model, vector<char> input_sequence, bool &verbose){
	float time1 = 0;
	float time2 = 0;
	float time3 = 0;   vector<Structure> temp_structures;
	Windower temp_windower;
	temp_windower.set_sequence(input_sequence);
	//check for circular sequence conditions
	if (input_sequence.size() == 0){
		//throw exception
	}
	//initializing the iterators ensures that the intial comparison in next_window_from_all_windows is not problematic
	std::vector<char>::iterator start = input_sequence.begin(),stop=input_sequence.begin()+1;
	temp_windower.reset_window();
	while (temp_windower.has_next_window()){
		temp_windower.next_window_from_all_windows(start,stop);
		Structure temp;
		float bp_energy = 0;
		//set the Loci of the structure using the gene's Loci
		
		temp.position.chromosome = position.chromosome;
		temp.position.strand = position.strand;
		temp.position.start_pos = position.start_pos + windower.get_current_start_offset();
		temp.position.end_pos = position.start_pos + windower.get_current_stop_offset();
		//pass the structure and window boundaries to the model
		model.compute_structure(sequence,start,stop,temp,bp_energy,time1,time2,time3,verbose);
		//push the now computed structure onto these_structures
		temp_structures.push_back(temp); //need to make sure the default copy constructor is working properly
	}
	return temp_structures;
}


void Gene::compute_structures_circular(Model &model, bool &verbose){ //not working
	float time1 = 0;
	float time2 = 0;
	float time3 = 0;    //check for circular sequence conditions
	if (sequence.size() == 0){
		//throw exception
		
	}
	//initializing the iterators ensures that the intial comparison in next_window_from_all_windows is not problematic
	std::vector<char>::iterator start = sequence.begin(),stop=sequence.begin()+1;
	windower.reset_window();
	while (windower.has_next_window_circular()){
		windower.next_window_from_all_windows_circular(start,stop);
		Structure temp;
		float bp_energy = 0;
		//set the Loci of the structure using the gene's Loci
		temp.position.chromosome = position.chromosome;
		temp.position.strand = position.strand;
		temp.position.start_pos = position.start_pos + windower.get_current_start_offset();
		temp.position.end_pos = position.start_pos + windower.get_current_stop_offset();
		//pass the structure and window boundaries to the model
		model.compute_structure(sequence,start,stop,temp,bp_energy,time1,time2,time3,verbose);
		//push the now computed structure onto these_structures
		rloop_structures.push_back(temp); //need to make sure the default copy constructor is working properly
	}
	//cout << rloop_structures.size() << endl;
}

void Gene::compute_residuals(Model &model, bool &verbose){
	//verify that the structures have been computed
	//iterate through all the structures
	for (int i=0; i < rloop_structures.size(); i++){ //not iterating through all the structures???
		model.compute_residuals(rloop_structures[i],verbose);
	}
}

void Gene::clear_structures(){
	rloop_structures.clear();
}

void Gene::complement_sequence(){
	for (int i = 0; i < sequence.size(); i++){
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

void Gene::reverse_sequence(){
	char temp;
	for (int i = 0; i < (sequence.size()/2); i++){
		temp = sequence[i];
		sequence[i] = sequence[sequence.size() - 1 - i];
		sequence[sequence.size() - 1 - i] = temp;
	}
}

// intentionally put here, I know it's the same as above and just re-writing
void Gene::reverse_complement_sequence() {
	complement_sequence();
	reverse_sequence();
}

int Gene::get_length(){
	return sequence.size();
}

void Gene::clear_sequence(){
	//delete sequence data
	sequence.clear();
}

void Gene::dump_structures(string outfilename) {
	ofstream dumpfile(outfilename + gene_name + "_dump.txt", ios::out);
	std::stringstream ss;
	ss << "start_position stop_position energy probability\n";
	ss << "0 0 " << ground_state_energy << ' ' << 0 << endl; //add the ground state probabilit
	cout << get_length() << endl;
	if (position.strand == "+") {
		for (int i = 0; i < rloop_structures.size(); i++) {
			ss << rloop_structures[i].position.start_pos << ' ' << rloop_structures[i].position.end_pos << ' ' <<
				rloop_structures[i].free_energy << ' ' << rloop_structures[i].probability << endl;
		}
	} else if (position.strand == "-") {
		for (int i = 0; i < rloop_structures.size(); i++) {
			ss << (this->getPosition().end_pos-rloop_structures[i].position.end_pos+this->getPosition().start_pos) << ' ' << (this->getPosition().end_pos-rloop_structures[i].position.start_pos+this->getPosition().start_pos) << ' ' <<
				rloop_structures[i].free_energy << ' ' << rloop_structures[i].probability << endl;
		}
	}
	else {
		cout << "Dump error. Strand unspecified.";
		exit(1); //replace with exception
	}
	dumpfile << ss.rdbuf();
	dumpfile.close();
}



// void Gene::parse_header(){
// 	unsigned long pos = 0;
// 	string name, remaining;

// 	/*
// 	Extract gene name, Example:
// 	>pFC9 range=pFC9gene:1-3637 5'pad=0 3'pad=0 strand=+ repeatMasking=none
// 	*/

// 	//delim '=' -> "pFC9 range" from "pFC9 range=pFC9gene..."
// 	pos  = header.find('=');       // @pos = split("=", $header);
// 	name = header.substr(0,pos);   // $name = $pos[0];
// 	// cout << "1. name=" << name << " pos=" << pos << endl;
// 	//delim ' ' -> "pFC9" from "pFC9 range"
// 	pos = header.find(' ');        // @pos = split(" ", $header);
// 	name = name.substr(1,pos);     // $name 
// 	// cout << "2. name=" << name << " pos=" << pos << endl;
// 	remaining = header.substr(pos+1, header.length());
// 	// cout << "3. Remaining=" << remaining << "\n";
	
// 	pos = name.find(' ');
// 	name = name.substr(0,pos);
// 	gene_name = name;

// 	//extract chromosome name
// 	pos = remaining.find('=');
// 	remaining = remaining.substr(pos+1, remaining.length());
// 	pos = remaining.find(':');
// 	position.chromosome = remaining.substr(0,pos);
// 	//for (auto & c: position.chromosome) c = toupper(c); //C++11 string toUpper
// 	std::transform(position.chromosome.begin(), position.chromosome.end(), position.chromosome.begin(), ::tolower);
// 	std::transform(position.chromosome.end()-1, position.chromosome.end(), position.chromosome.end()-1, ::toupper);
// 	remaining = remaining.substr(pos+1, remaining.length());
// 	//extract the start and stop locations
// 	pos = remaining.find("-");
// 	position.start_pos = stol(remaining.substr(0,pos));
// 	//position.start_pos = 0;//stol(remaining.substr(0,pos));
// 	remaining = remaining.substr(pos+1, remaining.length());
// 	pos = remaining.find(" ");
// 	position.end_pos = stol(remaining.substr(0,pos));
// 	cout << "SEQUENCE 3\n" << string(sequence.begin(),sequence.end()) << "\n";
// 	cout << "start=" << position.start_pos << "-" << position.end_pos << "\n";
// 	//position.end_pos = sequence.size(); //(remaining.substr(0,pos));
// 	// if (fabeg != 0) {
// 	// 	position.start_pos = fabeg;
// 	// }
// 	// if (faend != 0) {
// 	// 	position.end_pos = faend;
// 	// }
// 	remaining = remaining.substr(pos+1, remaining.length());
// 	pos = remaining.find("STRAND=");
// 	remaining = remaining.substr(pos+7, remaining.length());
// 	pos = remaining.find(" ");
// 	position.strand = remaining.substr(0,pos);
	
// 	return;
// }



// bool Gene::read_gene(ifstream& fastafile, int genebeg, int geneend) { //need to test

// 	fabeg = -1;
// 	faend = -1;
// 	//initialize variables
// 	char c,p;
// 	int currpos = 0;
// 	string tempseq = "";
// 	//check that fstream is open and not eof, throw exceptions otherwise
// 	if (!fastafile.is_open()) {
// 		throw UnexpectedClosedFileException("Gene::read_gene");
// 	} else if (fastafile.eof()) {
// 		throw UnexpectedEOFException();
// 	}
// 	while (fastafile.get(c)) {
// 		//read the next character
// 		c = toupper(c);
// 		//if the character is the start of a header line
// 		if (c == '>') {
// 			//read until the end of the header line
// 			while (c != '\n' && c != '\r') {
// 				header.push_back(c);
// 				c = toupper(fastafile.get());
// 			}
// 			parse_header();
// 		}
// 		else if (c == 'A' || c == 'T' || c == 'C' || c == 'G') {
// 			//save to sequence vector
// 			tempseq.push_back(c);
// 		}

// 		//else if encountered some white space
// 		else if (c == '\n' || c == ' ' || c == '\t' || c== '\r'){
// 			p = fastafile.peek();
// 			if (p == '>'){
// 				int tempbeg = fabeg;
// 				int tempend = faend;
				
// 				if (tempbeg != -1 && tempend != -1 && tempbeg >= tempend) {
// 					tempend = tempseq.length() + tempend;
// 					tempseq = tempseq + tempseq;
// 					cout << "tempbeg >= tempend so tempseq is doubled, tempbeg-tempend=" << tempbeg << "-" << tempend << ", length tempseq=" << tempseq.length() << "\n\n";
// 					for (size_t i = tempbeg; i <= tempend; i++) {
// 						char d = tempseq[i];
// 						sequence.push_back(d);
// 					}
// 				}
// 				else {
// 					if (fabeg == -1) {
// 						tempbeg = 0;
// 					}
// 					if (faend == -1) {
// 						tempend = tempseq.length()-1;
// 					}
// 					cout << "tempbeg < tempend so all good, tempbeg-tempend=" << tempbeg << "-" << tempend << ", length tempseq=" << tempseq.length() << "\n\n";
// 					for (size_t i = tempbeg; i <= tempend; i++) {
// 						char d = tempseq[i];
// 						sequence.push_back(d);

// 						// if ((fabeg != 0 && faend != 0 && fabeg >= faend) || (fabeg != 0 && faend == 0) || (fabeg == 0 && faend != 0) || (fabeg == 0 && faend == 0)) {
// 						// 	currpos ++;
// 						// 	if (fabeg != 0 && faend != 0) {
// 						// 		if (currpos >= fabeg && currpos <= faend) {
// 						// 			sequence.push_back(d);
// 						// 		}
// 						// 	}
// 						// 	else if (fabeg != 0 && currpos >= fabeg) {
// 						// 		sequence.push_back(d);
// 						// 	}
// 						// 	else if (faend != 0 && currpos <= faend) {
// 						// 		sequence.push_back(d);
// 						// 	}
// 						// 	else {
// 						// 		sequence.push_back(d);
// 						// 	}
// 						// }
// 					//std::cout << d << " ";
// 					}
// 				}

// 				windower.set_sequence(sequence);
// 				return false;
// 			}
// 			/*else if (p == EOF) {
// 			 return true;
// 			}*/
// 			continue;
// 		}
// 		//else the charicter is unrecognized
// 		else {
// 			throw InvalidSequenceDataException(c);
// 		}
// 	}

// 	int tempbeg = fabeg;
// 	int tempend = faend;
// 	string tempseqorig = tempseq;
	
// 	// cout << "fabeg=" << fabeg << "\n";
// 	// cout << "faend=" << faend << "\n";

// 	if (fabeg < 0) {
// 		tempbeg = 0;
// 	}
// 	if (faend > tempseq.length()) {
// 		tempend = tempseq.length();
// 	}
// 	if (fabeg != -1 && faend != -1 && fabeg >= faend) {
// 		tempend = tempseq.length() + tempend;
// 		tempseq = tempseq + tempseq;
// 	}

// 	for (size_t i = tempbeg; i < tempend; i++) {
// 		char d = tempseq[i];
// 		sequence.push_back(d);
// 	}
// 	// cout << "tempbeg >= tempend so CIRCULAR, tempbeg-tempend=" << tempbeg << "-" << tempend << ", length tempseq2=" << tempseq2.length() << ", lengh sequence=" << sequence.size() << "\n\n";
// 	// }
// 	// else {
// 	// 	for (size_t i = tempbeg; i < tempend; i++) {
// 	// 		char d = tempseq[i];
// 	// 		sequence.push_back(d);
// 	// 	}
// 	// 	// cout << "tempbeg < tempend so NOT CIRCULAR, tempbeg-tempend=" << tempbeg << "-" << tempend << ", length tempseq=" << tempseq.length() << ", lengh sequence=" << sequence.size() << "\n\n";
// 	// }

// 	// exit(1);	
// 	string mydot = std::string(tempbeg,'-');
// 	// cout << "gene.cpp: original sequence:\n" << tempseq << "\n";
// 	// cout << "gene.cpp: beg-end (" << tempbeg << "-" << tempend << "):\n" << mydot << string(sequence.begin(),sequence.end()) << "\n\n";	
// 	windower.set_sequence(sequence);
// 	cout << "gene.cpp: 1. position.start_pos     - position.end_pos (original): " << position.start_pos << "-" << position.end_pos << "\n";
// 	position.start_pos = 1;
// 	position.end_pos = sequence.size();// - 1;
// 	cout << "gene.cpp: 2. position.start_pos (1) - position.end_pos (sequence.size()): " << position.start_pos << "-" << position.end_pos << "\n";

// 	if (genebeg != -1) {
// 		fabeg = genebeg;
// 	}
// 	else {
// 		fabeg = 1;
// 	}
// 	if (geneend != -1) {
// 		faend = geneend;
// 	}
// 	else {
// 		faend = sequence.size();
// 	}


// 	return true;
// 	//unexpected EOF
// 	//throw InvalidSequenceDataException();
// }
