// Created by Robert Stolz on 6/28/17.
// Edited and updated by Stella 2023-2025
//

#include "mytime.h"
#include "simulation.h" //

typedef std::chrono::high_resolution_clock::time_point TimeVar;


void initialize_models(char *model)
{
}

Simulation::Simulation()
{
	// default constructor
	minlength = 2; 
	dynamic_window_size = 15;
	reverse_flag = false;
	complement_flag = false;
	power_threshold = 1;
	circular_flag = false;
	auto_domain_size = false;
	import_flag = false;
	top = 0;
	dump = false;
	average_g = true; // false;
	seed = 0;
	dynamic_flag = false;
	naive_flag = false;
	verbose_flag = false;
	orig_flag = false;
}

Simulation::~Simulation()
{
	for (std::vector<Gene *>::iterator it = genes.begin(); std::distance(genes.begin(), it) < genes.size(); ++it)
	{
		delete *it; 
	}
}

void Simulation::set_verbose(bool verbose)
{
	verbose_flag = true;
}

void Simulation::set_dynamic(string dynamic)
{
	if (dynamic == "False") {
		cout << "NOT dynamic\n";
		dynamic_flag = false;
	}
	else {
		cout << "Dynamic is set!" << "\n";
		dynamic_flag = true;
	}
	
}

void Simulation::set_naive(bool naive)
{
	if (naive) {
		cout << "Naive is set!" << "\n";
	}
	else {
		cout << "NOT naive\n";
	}
	naive_flag = naive;
}

void Simulation::set_orig(bool orig)
{
	if (orig) {
		cout << "Orig is set!" << "\n";
	}
	else {
		cout << "NOT orig\n";
	}
	orig_flag = orig;
}

void Simulation::set_infile(string Infilename)
{
	infilename = Infilename;
	infile.open(infilename, ios::in);
}

void Simulation::set_outfile(string Outfilename)
{
	outfilename = Outfilename;
}

void Simulation::set_outdir(string Outdir)
{
	outdir = Outdir;
}

void Simulation::set_minlength(int Minlength)
{
	minlength = Minlength;
}

void Simulation::dynamic_setWindow_size(int Window_Size)
{
	dynamic_window_size = Window_Size;
}

void Simulation::setGeneBeg(int value)
{
	genebeg = value;
}

void Simulation::setGeneEnd(int value)
{
	geneend = value;
}

void Simulation::set_printbedfile(bool value)
{
	printbedfile = value;
}

void Simulation::set_power_threshold(int Power_threshold)
{
	power_threshold = Power_threshold;
}

void Simulation::set_circular(bool value)
{
	circular_flag = value;
}

void Simulation::set_import_flag(string filename)
{
	if (filename == "False") {
		import_flag = false;
	}
	else {
		import_flag = true;
		importfilename = filename;
	}
}

void Simulation::set_residuals(bool value)
{
	residuals = value;
}

void Simulation::set_auto_domain_size(int value)
{
	auto_domain_size = false;
	if (value == 0) {
		auto_domain_size = true;
	}
}

void Simulation::set_dump(bool value)
{
	dump = value;
}

void Simulation::set_average_g(bool value)
{
	average_g = value;
}

void Simulation::set_seed(int value)
{
	seed = value;
	srand(seed);
}

void Simulation::reverse_input(bool flag)
{
	if (flag) {
		cout << "reverse input\n";
	}
	reverse_flag = flag;
}

void Simulation::complement_input(bool flag)
{
	if (flag) {
		cout << "complement input\n";
	}
	complement_flag = flag;
}

void Simulation::set_top(int n)
{
	top = n;
}

std::vector<Model *> Simulation::get_models()
{
	return models;
}

void Simulation::add_model(Model &model)
{
	models.push_back(&model);
}

vector<Peak> Simulation::import_external_structures(string importfilename, Model &model)
{
	ifstream infile(importfilename, ios::in);
	vector<Peak> temp;
	long int start, stop;
	float energy;
	string chr, sign;
	while (infile >> chr >> start >> stop >> sign >> energy)
	{
		temp.push_back(Peak(Loci(chr, sign, start, stop), energy));
	}
	return temp;
}

void Simulation::compute_signal_extbpprobs(Gene &gene, vector<double> *&signal)
{
	signal = new vector<double>(gene.get_length(), 0.0);
	// compute the r-loop involvement probability for each base
	// for each structure in the gene
	for (std::vector<Structure>::iterator it = gene.getRloopStructures().begin();
		 it < gene.getRloopStructures().end(); ++it)
	{
		if (it->external)
		{
			// for each base in the structure
			for (long int i = it->position.start_pos - gene.getPosition().start_pos;
				 i < it->position.end_pos - gene.getPosition().start_pos; i++)
			{
				(*signal)[i] += it->probability;
			}
		}
	}
	// if strand is -, reverse bp_probabilities
	if (gene.getPosition().strand == "-")
	{
		std::reverse(signal->begin(), signal->end());
	}
}

void Simulation::compute_signal_bpprobs(Gene &gene, vector<double> *&signal)
{
	string LPR = "\033[1;95m";
	string LBU = "\033[1;94m";
	string LRD = "\033[1;31m";
	string LCY = "\033[1;36m";
	string YW = "\033[1;33m";
	string LGN = "\033[1;32m";
	string NN = "\033[1;0m";

	signal = new vector<double>(gene.get_length(), 0.0);

	TimeVar t1 = timeNow();
	int iter0 = 1;
	int n = 1;
	int m = 2;
	int last_j = -1;
	int curr_signal = 0;
	for (std::vector<Structure>::iterator it = gene.getRloopStructures().begin() - 1; it >= gene.getRloopStructures().end(); ++it)
	{
		long int j = it->position.start_pos - gene.getPosition().start_pos;
		long int j1 = it->position.end_pos - gene.getPosition().start_pos;
		for (long int i = it->position.start_pos - gene.getPosition().start_pos; i < it->position.end_pos - gene.getPosition().start_pos; i++)
		{
			(*signal)[i] += it->probability;
			j++;
		}
	}
	cout << "\n";
	cout << "- " << YW << "sim.cpp:: compute signal bpprobs(): TOTAL O(n*m) compute_structures: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
	t1 = timeNow(); 

	if (gene.getPosition().strand == "-")
	{
		std::reverse(signal->begin(), signal->end());
	}
}

void Simulation::compute_signal_average_G(Gene &gene, vector<double> *&signal)
{
	signal = new vector<double>(gene.get_length(), 0.0);
	vector<double> bp_partition_functions(gene.get_length(), 0.0);
	for (std::vector<Structure>::iterator it = gene.getRloopStructures().begin();
		 it < gene.getRloopStructures().end(); ++it)
	{
		for (long int i = it->position.start_pos - gene.getPosition().start_pos;
			 i < it->position.end_pos - gene.getPosition().start_pos; i++)
		{
			bp_partition_functions[i] += it->boltzmann_factor;
		}
	}
	for (std::vector<Structure>::iterator it = gene.getRloopStructures().begin();
		 it < gene.getRloopStructures().end(); ++it)
	{
		for (long int i = it->position.start_pos - gene.getPosition().start_pos;
			 i < it->position.end_pos - gene.getPosition().start_pos; i++)
		{
			(*signal)[i] += (it->boltzmann_factor / bp_partition_functions[i]) * it->free_energy;
		}
	}
	if (gene.getPosition().strand == "-")
	{
		std::reverse(signal->begin(), signal->end());
	}
}

void Simulation::compute_signal_mfe(Gene &gene, vector<double> *&signal)
{
	signal = new vector<double>(gene.get_length(), 0.0);
	double current_min = FLT_MAX;
	Structure mfe;
	for (std::vector<Structure>::iterator it = gene.getRloopStructures().begin();
		 it < gene.getRloopStructures().end(); ++it)
	{
		if (it->free_energy < current_min)
		{
			current_min = it->free_energy;
			mfe = *it;
		}
	}
	for (long int i = mfe.position.start_pos - gene.getPosition().start_pos;
		 i < mfe.position.end_pos - gene.getPosition().start_pos; i++)
	{
		(*signal)[i] = 1.0;
	}
	if (gene.getPosition().strand == "-")
	{
		std::reverse(signal->begin(), signal->end());
	}
}

void Simulation::call_peaks_threshold(Gene &gene, vector<double> &signal, vector<Loci> &peaks)
{
	double minimum = 1;
	bool in_peak = false;
	long peak_start = 0, peak_end = 0;
	double magnitude = 0;
	Structure *temp;
	for (int i = 0; i < signal.size(); i++)
	{
		if (signal[i] < minimum && signal[i] != 0)
		{
			minimum = signal[i];
		}
	}
	for (int i = 0; i < signal.size(); i++)
	{
		if (signal[i] > minimum * pow(10, power_threshold))
		{ // the signal is significant
			if (!in_peak)
			{
				in_peak = true;
				peak_start = gene.getPosition().start_pos + i;
			}
		}
		else
		{ // the signal is not significant
			if (in_peak)
			{
				in_peak = false;
				peak_end = gene.getPosition().start_pos + i;
				peaks.emplace_back(Loci(gene.getPosition().chromosome, gene.getPosition().strand, peak_start, peak_end)); 
			}
		}
	}
}

void Simulation::call_peaks_absolute_threshold(Gene &gene, vector<double> &signal, vector<Loci> &peaks)
{
	double minimum = 1;
	bool in_peak = false;
	long peak_start = 0, peak_end = 0;
	double magnitude = 0;
	Structure *temp;
	for (int i = 0; i < signal.size(); i++)
	{
		if (signal[i] > 1 * pow(10, power_threshold))
		{ // the signal is significant
			if (!in_peak)
			{
				in_peak = true;
				peak_start = gene.getPosition().start_pos + i;
			}
		}
		else
		{ // the signal is not significant
			if (in_peak)
			{
				in_peak = false;
				peak_end = gene.getPosition().start_pos + i;
				peaks.emplace_back(Loci(gene.getPosition().chromosome, gene.getPosition().strand, peak_start, peak_end)); 
			}
		}
	}
}

void Simulation::cluster_k_intervals(vector<Loci> &peaks, vector<Loci> &clustered_peaks)
{
	long long int seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::cout << "rng seed: " << seed << "\n";
	vector<double> costs;
	vector<int> chosen_peaks;
	vector<int> clustering_tally;
	for (int i = 0; i < peaks.size(); i++)
	{
		clustering_tally.push_back(0);
	}
	int k;
	k = 5;
	for (int i = 0; i < 1000; i++)
	{
		lloyds_algorithm(peaks, chosen_peaks, k, seed);
		for (int j = 0; j < chosen_peaks.size(); j++)
		{
			clustering_tally[chosen_peaks[j]]++;
		}
		chosen_peaks.empty();
	}
}

double Simulation::lloyds_algorithm(vector<Loci> &peaks, vector<int> &clustering, int k, unsigned seed)
{
	bool swaps = true;
	vector<int> medoid_indeces;		// maps medoid index to actual element in the matrix
	vector<int> medoid_assignments; // the INDEX of the medoid each peak is assigned to.
	vector<vector<double>> pairwise_distance_matrix;
	double configuration_cost = 0;
	for (int i = 0; i < peaks.size(); i++)
	{ // initialize the pairwise distance matrix
		vector<double> temp;
		for (int j = 0; j < peaks.size(); j++)
		{
			temp.push_back(0);
		}
		pairwise_distance_matrix.push_back(temp);
	}
	// choose k different intervals at random as the initial medoids
	// generate k random indeces
	vector<int> shuffled;
	for (int i = 0; i < peaks.size(); i++)
	{ // unshuffled medoid indeces
		shuffled.push_back(i);
		medoid_assignments.push_back(0); // all peaks are temporarily assigned to the first medoid
	}
	std::shuffle(shuffled.begin(), shuffled.end(), std::default_random_engine(seed)); // not tested, need to connect the seed
	for (int i = 0; i < k; i++)
	{
		medoid_indeces.push_back(shuffled[i]); // save the k randomly selected medoid indeces to a list
	}
	// compute the pairwise distance matrix
	for (int i = 0; i < peaks.size(); i++)
	{ // for each peak
		for (int j = 0; j < peaks.size(); j++)
		{ // for each peak
			pairwise_distance_matrix[i][j] = interval_distance(peaks[i], peaks[j]);
		}
	}
	double current_cost = 0; // cost of the current clustering configuration
	// assign each interval to its closest medoid
	for (int i = 0; i < peaks.size(); i++)
	{ // for each peak
		for (int j = 1; j < k; j++)
		{ // for each medoid index
			if (pairwise_distance_matrix[i][medoid_indeces[j]] < pairwise_distance_matrix[i][medoid_indeces[medoid_assignments[i]]])
			{
				medoid_assignments[i] = j;
			}
		}
	}
	// compute full configuration cost
	for (int i = 0; i < medoid_assignments.size(); i++)
	{
		configuration_cost += pairwise_distance_matrix[i][medoid_indeces[medoid_assignments[i]]];
	}

	while (swaps)
	{ // Veroni descent
		swaps = false;
		// assign each interval number to its closest medoid (already done for the first iteration)
		for (int i = 0; i < peaks.size(); i++)
		{
			for (int j = 1; j < k; j++)
			{
				if (pairwise_distance_matrix[i][medoid_indeces[j]] < pairwise_distance_matrix[i][medoid_assignments[i]])
				{
					// configuration_cost -= pairwise_distance_matrix[i][medoid_indeces[medoid_assignments[i]]]; //update the configuration cost
					medoid_assignments[i] = j; // update the medoid assignment with the index of the new medoid
											   // configuration_cost += pairwise_distance_matrix[i][medoid_indeces[medoid_assignments[i]]];
				}
			}
		}
		// for each cluster
		for (int p = 0; p < k; p++)
		{
			// test each object within the cluster as the new medoid of the cluster
			for (int i = 0; i < peaks.size(); i++)
			{
				if (medoid_assignments[i] == p && i != medoid_indeces[p])
				{ // if the medoid is in the currently considered group, but is not the current medoid
					// determine swap cost
					double costA = 0, costB = 0;
					for (int j = 0; j < peaks.size(); j++)
					{
						if (medoid_assignments[i] == p)
						{																				 // if element is in the currently considered cluster
							costA += pairwise_distance_matrix[medoid_indeces[medoid_assignments[i]]][j]; // current configuration
							costB += pairwise_distance_matrix[i][j];									 // currently considered swap
						}
						if (costB < costA)
						{ // swap would reduce the configuration cost
							// update the configuration cost
							configuration_cost -= costA;
							configuration_cost += costB;
							// update medoid_indeces
							medoid_indeces[p] = i;
						}
					}
				}
			}
		}
	}
	// tally the final clustering
	clustering = medoid_indeces;
	return configuration_cost;
}

double Simulation::compute_configuration_cost(vector<vector<double>> &pairwise_distance_matrix, vector<int> medoid_indeces)
{
	double configuration_cost = 0;
	for (int i = 0; i < pairwise_distance_matrix.size(); i++)
	{
		for (int j = 0; j < medoid_indeces.size(); j++)
		{
			configuration_cost += pairwise_distance_matrix[i][medoid_indeces[j]];
		}
	}
	return configuration_cost;
}

double Simulation::interval_distance(const Loci &A, const Loci &B)
{
	double term1 = pow((A.start_pos + A.end_pos) / 2. - (B.start_pos + B.end_pos) / 2., 2);
	double term2 = pow((A.end_pos - A.start_pos) / 2. - (B.end_pos - B.start_pos) / 2., 2) / 3.;
	return term1 + term2;
}

void Simulation::write_wigfile_header(ofstream &outfile, string trackname)
{
	// open stringstream
	std::stringstream ss;
	// compose .wig header
	// adjust browser position
	ss << "track type=wiggle_0 name=\"" << trackname << "\" visibility=full autoscale=off color=50,150,255 priority=10"
	   << "\n";
	outfile << ss.rdbuf();
}

void Simulation::write_wigfile(ofstream &outfile, Gene *gene, std::vector<double> *signal)
{
	// open stringstream
	std::stringstream ss;
	string wigfile_name = gene->getHeader().c_str();
	// compose .wig header
	string name = gene->getName();
	// adjust browser position
	ss << "browser position " << gene->getPosition().chromosome << ':' << gene->getPosition().start_pos << '-' << gene->getPosition().end_pos << "\n";
	ss << '#' << gene->getName() << "\n";
	ss << "fixedStep chrom=" << gene->getPosition().chromosome << " start=" << gene->getPosition().start_pos << " step=1"
	   << "\n";
	for (int i = 0; i < signal->size(); i++)
	{
		ss << (*signal)[i] << "\n";
	}
	// write stringstream to file
	outfile << ss.rdbuf();
}

void Simulation::read_bedfile(ifstream &bedinput, vector<Loci> &peaks)
{
	Loci temp;
	long int pos;
	char buffer[1000];
	string strbuff;
	if (!bedinput.is_open())
	{
		// throw exception
	}
	while (bedinput.getline(buffer, 1000))
	{
		strbuff = std::string(buffer);
		// need to deal with lines that do not contain a bed entry here

		// parse out chromosome name
		pos = strbuff.find('\t');
		temp.chromosome = strbuff.substr(0, pos); // need to handle non-numeric chromosome names as well
		strbuff = strbuff.substr(pos + 1, strbuff.length());
		// parse out start position of the entry
		pos = strbuff.find('\t');
		temp.start_pos = stol(strbuff.substr(0, pos));
		strbuff = strbuff.substr(pos + 1, strbuff.length());
		// parse out end position of the entry
		pos = strbuff.find('\t');
		temp.end_pos = stol(strbuff.substr(0, pos));
		strbuff = strbuff.substr(pos + 1, strbuff.length());
		// discard the next two columns (may need to be made more flexible in the future)
		pos = strbuff.find('\t');
		strbuff = strbuff.substr(pos + 1, strbuff.length());
		pos = strbuff.find('\t');
		strbuff = strbuff.substr(pos + 1, strbuff.length());
		// parse out the strand
		pos = strbuff.find('\t');
		temp.strand = strbuff.substr(0, pos);
		strbuff = strbuff.substr(pos + 1, strbuff.length());
		// save to the peaks vector
		peaks.push_back(temp);
	}
}

void Simulation::write_bedfile_header(ofstream &outfile, string trackname)
{
	// write bedfile
	stringstream ss;
	ss << "track name=rLooper description=\"" << trackname << "\" useScore=1"
	   << "\n";
	outfile << ss.rdbuf();
}

void Simulation::write_bedfile(ofstream &outfile, Gene *gene, vector<Loci> &peaks)
{
	// write bedfile
	stringstream ss;
	string strand_name;
	int start_pos = 0, end_pos = 0;
	if (gene->getPosition().strand == "+")
	{
		strand_name = "POS";
	}
	else
	{
		strand_name = "NEG";
	}
	ss << "browser position " << gene->getPosition().chromosome << ':' << gene->getPosition().start_pos << '-' << gene->getPosition().end_pos << "\n";
	ss << '#' << gene->getName() << "\n";
	// print BED header here
	// print the peaks in BED format
	for (int i = 0; i < peaks.size(); i++)
	{
		ss << peaks[i].chromosome << '\t' << (peaks)[i].start_pos << '\t' << peaks[i].end_pos
		   << '\t' << strand_name << i << '\t' << '0' << '\t' << peaks[i].strand << "\n";
	}
	// write stringstream to file
	outfile << ss.rdbuf();
}

string join(char &delim, const vector<char> &lst)
{
	string join_this_vector_char;
	int i = 0;
	for (const auto &s : lst)
	{
		if (i == 0)
		{
			join_this_vector_char += delim;
			i = 1;
		}
		join_this_vector_char += s;
	}
	return join_this_vector_char;
}

string join(const vector<char> &lst)
{
	string join_this_vector_char;
	for (const auto &s : lst)
	{
		join_this_vector_char += s;
	}
	return join_this_vector_char;
}

/*naive rlooper*/
//-------------------------{
// same as Rloop_model::compute_bps_interval

double getB(const char &first, const char &second)
{

	if (first == 'G')
	{					   // G
		if (second == 'G') // GG
			return -0.36;
		else if (second == 'C') // GC
			return -0.16;
		else if (second == 'A') // GA
			return -0.1;
		else if (second == 'T' || second == 'U') // GT/GU
			return -0.06;
		else if (second == 'N') // GN
			return 0;
	}
	else if (first == 'C')
	{					   // G
		if (second == 'G') // CG
			return 0.97;
		else if (second == 'C') // CC
			return 0.34;
		else if (second == 'A') // CA
			return 0.45;
		else if (second == 'T' || second == 'U') // CT/CU
			return 0.38;
		else if (second == 'N') // CN
			return 0;
	}
	else if (first == 'A')
	{					   // A
		if (second == 'G') // AG
			return -0.12;
		else if (second == 'C') // AC
			return -0.16;
		else if (second == 'A') // AA
			return 0.6;
		else if (second == 'T' || second == 'U') // AG/AU
			return -0.12;
		else if (second == 'N') // AN
			return 0;
	}
	else if (first == 'T')
	{					   // T
		if (second == 'G') // TG
			return 0.45;
		else if (second == 'C') // TC
			return 0.5;
		else if (second == 'A') // TA
			return 0.28;
		else if (second == 'T' || second == 'U') // TT/TU
			return 0.8;
		else if (second == 'N') // TN
			return 0;
	}
	else if (first == 'N')
	{					   // N
		if (second == 'G') // NG
			return 0;
		else if (second == 'C') // NC
			return 0;
		else if (second == 'A') // NA
			return 0;
		else if (second == 'T' || second == 'U') // NT/NU
			return 0;
		else if (second == 'N') // NN
			return 0;
	}
	// if nothing return at this point, then die coz there's unexpected sequence that's non-A/C/G/T/U
	throw std::invalid_argument("There is a non-A/C/G/T/U/N in the sequence!\n");
	exit(1);
}


double getRB(const char &first, const char &second)
{

	if (first == 'G')
	{					   // G
		if (second == 'G') // GG
			return -2.2;
		else if (second == 'C') // GC
			return -2.4;
		else if (second == 'A') // GA
			return -1.4;
		else if (second == 'T' || second == 'U') // GT/GU
			return -1.5;
		else if (second == 'N') // GN
			return 0;
	}
	else if (first == 'C')
	{					   // G
		if (second == 'G') // CG
			return -1.2;
		else if (second == 'C') // CC
			return -1.5;
		else if (second == 'A') // CA
			return -1.0;
		else if (second == 'T' || second == 'U') // CT/CU
			return -0.9;
		else if (second == 'N') // CN
			return 0;
	}
	else if (first == 'A')
	{					   // A
		if (second == 'G') // AG
			return -1.4;
		else if (second == 'C') // AC
			return -1.6;
		else if (second == 'A') // AA
			return -0.4;
		else if (second == 'T' || second == 'U') // AG/AU
			return -1.0;
		else if (second == 'N') // AN
			return 0;
	}
	else if (first == 'T')
	{					   // T
		if (second == 'G') // TG
			return -1.0;
		else if (second == 'C') // TC
			return -0.8;
		else if (second == 'A') // TA
			return -0.3;
		else if (second == 'T' || second == 'U') // TT/TU
			return -0.2;
		else if (second == 'N') // TN
			return 0;
	}
	else if (first == 'N')
	{					   // N
		if (second == 'G') // NG
			return 0;
		else if (second == 'C') // NC
			return 0;
		else if (second == 'A') // NA
			return 0;
		else if (second == 'T' || second == 'U') // NT/NU
			return 0;
		else if (second == 'N') // NN
			return 0;
	}
	// if nothing return at this point, then die coz there's unexpected sequence that's non-A/C/G/T/U
	throw std::invalid_argument("There is a non-A/C/G/T/U/N in the sequence!\n");
	exit(1);
}


double getB2(const char &first, const char &second)
{

	if (first == 'G')
	{					   // G
		if (second == 'G') // GG
			return -2.2 - (-1.84/2);
		else if (second == 'C') // GC
			return -2.4 - (-2.24/2);
		else if (second == 'A') // GA
			return -1.4 - (0.5*-1.3);
		else if (second == 'T' || second == 'U') // GT/GU
			return -1.5 - (0.5*-1.44);
		else if (second == 'N') // GN
			return 0 - (0.5*-0);
	}
	else if (first == 'C')
	{					   // G
		if (second == 'G') // CG
			return -1.2 - (0.5*-2.17);
		else if (second == 'C') // CC
			return -1.5 - (0.5*-1.84);
		else if (second == 'A') // CA
			return -1.0 - (0.5*-1.45);
		else if (second == 'T' || second == 'U') // CT/CU
			return -0.9 - (0.5*-1.28);
		else if (second == 'N') // CN
			return 0 - (0.5*-0);
	}
	else if (first == 'A')
	{					   // A
		if (second == 'G') // AG
			return -1.4 - (0.5*-1.28);
		else if (second == 'C') // AC
			return -1.6 - (0.5*-1.44);
		else if (second == 'A') // AA
			return -0.4 - (0.5*-1.00);
		else if (second == 'T' || second == 'U') // AG/AU
			return -1.0 - (0.5*-0.88);
		else if (second == 'N') // AN
			return 0 - (0.5*-0);
	}
	else if (first == 'T')
	{					   // T
		if (second == 'G') // TG
			return -1.0 - (0.5*-1.45);
		else if (second == 'C') // TC
			return -0.8 - (0.5*-1.30);
		else if (second == 'A') // TA
			return -0.3 - (0.5*-0.58);
		else if (second == 'T' || second == 'U') // TT/TU
			return -0.2 - (0.5*-1.00);
		else if (second == 'N') // TN
			return 0;
	}
	else if (first == 'N')
	{					   // N
		if (second == 'G') // NG
			return 0;
		else if (second == 'C') // NC
			return 0;
		else if (second == 'A') // NA
			return 0;
		else if (second == 'T' || second == 'U') // NT/NU
			return 0;
		else if (second == 'N') // NN
			return 0;
	}
	// if nothing return at this point, then die coz there's unexpected sequence that's non-A/C/G/T/U
	throw std::invalid_argument("There is a non-A/C/G/T/U/N in the sequence!\n");
	exit(1);
}

int stupid_number_converter_part1(const double &number)
{
	if (number <= DBL_MIN)
	{
		return (0);
	}
	int ten_power_this = log(number) / log(10) - 3;
	return (ten_power_this);
}

int stupid_number_converter_part2(const double &number, const int &ten_power_this)
{
	if (number <= DBL_MIN)
	{
		return (0);
	}
	int front_number = (number) / pow(10, ten_power_this) + 0.5;
	return (front_number);
}

void printer(stringstream &ss7, const double &number)
{
	if (number <= DBL_MIN)
	{
		ss7 << "\t" << 0;
	}
	else
	{
		int num2 = stupid_number_converter_part1(number);
		int num1 = stupid_number_converter_part2(number, num2);
		ss7 << "\t" << num1 << "e" << num2;
	}
}

void printerfinal(ofstream &ss, int fabeg, const bool &printSlowButCorrect, const bool &printQuickDirty, const int &n, const int &len, const char &n1, const char &n2, const float &Gbpsindex, const double &Gsigmasindex, const double &Gsindex, const double &bfsindex, const double &bftotal)
{
	std::stringstream ss7;
	int len2 = len;
	if (len2 != 0)
	{
		len2++;
	}
	ss7 << fabeg + n + 1 << "\t" << len2 << "\t" << n1 << "\t" << n2;
	if (printSlowButCorrect)
	{
		ss7 << "\t" << Gbpsindex << "\t" << Gsigmasindex << "\t" << Gsindex << "\t" << bfsindex << "\t" << bfsindex / bftotal << "\n";
	}
	else if (printQuickDirty)
	{
		int Gsindex2 = 10 * Gsindex;
		int Gsigmasindex2 = 10 * Gsigmasindex;
		int Gbpsindex2 = 10 * Gbpsindex;
		ss7 << "\t" << Gbpsindex2 << "e-1";
		ss7 << "\t" << Gsigmasindex2 << "e-1";
		ss7 << "\t" << Gsindex2 << "e-1";
		if (bfsindex <= DBL_MIN)
		{
			ss7 << "\t0\t0";
		}
		else
		{
			printer(ss7, bfsindex);
			printer(ss7, bfsindex / bftotal);
		}
		ss7 << "\n";
	}
	ss << ss7.rdbuf();
}

void naive_forloop_rlooper(const vector<char> &myseq, int fabeg, ofstream &outfile_naive_bigtable, ofstream &outfile_naive_smalltable, ofstream &outfile_naive_sim_reads, ofstream &outfile_naive_sim_count, const double &my_K, const double &my_a, const double &my_sigma, const int &simN, const int &nick0, const int &nicklen, const int &selffoldlen, bool &verbose, const int &minlength)
{
	const char *LRD = "\033[31m";
	const char *LCY = "\033[36m";
	const char *YW = "\033[33m";
	const char *LGN = "\033[32m";
	const char *NN = "\033[0m";
	
	int nick = nick0;
	if (fabeg != -1) {
		nick = nick - fabeg;
	}
			
	
	// 0. printing options (cout)
	bool printSlowButCorrect = true;
	bool printSlowButCorrect2 = false;
	bool printQuickDirty = false;
	TimeVar t1 = timeNow();

	// use prob or Gs
	bool useprobs = true;
	// 1. Precalculate Gsigma and make a table. This sped up by 5-10% (0.2s in 3s runtime)
	/////////////////////////////////////////////////////////////
	//             2 * (pi^2) * C * K                          //
	// Gsigma = ------------------------ * (alpha + m * A)^2   //
	//             4 * (pi^2) * C + K * m                      //
	/////////////////////////////////////////////////////////////
	vector<float> Gsigma;
	int max_length = 2000;
	const float my_alpha = (1500 * my_sigma * 1 / 10.4);
	for (int m = 0; m < myseq.size(); m++)
	{
		int mact = m;
		if (m != 0)
		{
			mact = m + 1; // because using m at least 2bp. Not using m == 1 bp
		}
		Gsigma.push_back(

			(2 * pow(pi, 2) * 1.8 * my_K) * pow((my_alpha + mact * 1 / 10.4), 2) /
			(4 * pow(pi, 2) * 1.8 + my_K * mact)

		);
	}
	if (verbose)
	{
		cout << YW << "  - naive_forloop_rlooper:: O(m) for loop, make Gsigma: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
		t1 = timeNow();
	}

	// 2. For loop O(n(n+1)/2) ~ O(n^2)
	// - Gs: Energy = G.a + sum(G.bp,[1..m]) + G.sigma
	// - bfs: boltzmann factor (exp(-1*Gs))
	// Notes: Storing these as float or int or double doesn't change runtime
	vector<double> Gsigmas;
	vector<float> Gbps;
	vector<double> Gs;	// probably only need float but just in case
	vector<double> bfs; // definitely need double unless we don't care about 1e-40s and below
	vector<char> n1;
	vector<char> n2;
	int index = 0;
	double bftotal = 0;
	double Gs_m0 = (0 + 0 + Gsigma[0]);
	double bf_m0 = exp(-1 * Gs_m0 / (0.0019858775 * 310));
	// total row: N(N+1)/2
	
	for (int n = 0; n < myseq.size(); n++)
	{
		float Gbp = 0;

		// no-rloop (total number of 'no-Rloop' state = 'myseq.size' = N)
		if (n == 0) {
			n1.push_back('-');
			n2.push_back(myseq[n]);
			Gbps.push_back(0);
			Gsigmas.push_back(Gs_m0);
			Gs.push_back(Gs_m0);
			bfs.push_back(bf_m0);
			bftotal += bf_m0;
			index++;
		}

		for (int m = 0; n + m < myseq.size() - 1; m++)
		{
			if (m > max_length - 1)
			{
				break;
			}

			double curr_a = my_a;

			if (n >= nick && n+m >= nick && n+m < nick + selffoldlen)
			{
				Gbp += getB(myseq[n + m], myseq[n + m + 1]);
				// Gbp += getRB(myseq[n + m], myseq[n + m + 1]); // for RNA DNA?
			}
			else {
				Gbp += getB(myseq[n + m], myseq[n + m + 1]);
			}

			if (n >= nick && n < nick + nicklen)
			{
				curr_a = 0;
			}
			n1.push_back(myseq[n + m]);
			n2.push_back(myseq[n + m + 1]);
			
			
			Gbps.push_back(Gbp);
			Gsigmas.push_back(Gsigma[m + 1]);
			Gs.push_back((curr_a + Gbp + Gsigma[m + 1]));
			bfs.push_back(exp(-1 * ((curr_a + Gbp + Gsigma[m + 1]) / (0.0019858775 * 310))));
			bftotal += exp(-1 * ((curr_a + Gbp + Gsigma[m + 1]) / (0.0019858775 * 310)));
			
			index++;
		}
	}

	double probtotal = 0;
	if (verbose)
	{
		cout << YW << "  - naive_forloop_rlooper:: O(n*m/2): sum Gbp, sum to G (a+sum(Gbp)+Gsigma), bfs, sum part_func (bftotal) : " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
		t1 = timeNow();
	} // TIME

	index = 0; // reset index to 0
	vector<double> probs;
	vector<int> ms;
	vector<int> ns;
	// tC8 valBuff[20]; // for printSlowButCorrect2
	for (int n = 0; n < myseq.size(); n++)
	{
		if (n == 0) {
			ns.push_back(n + 1);
			ms.push_back(0);
			probs.push_back(bfs[index] / bftotal);
			probtotal += bfs[index] / bftotal;
			index++;
		}

		for (int m = 0; n + m < myseq.size() - 1; m++)
		{
			if (m > max_length - 1)
			{
				break;
			}
			ns.push_back(n + 1);
			ms.push_back(m + 1);
			probs.push_back(bfs[index] / bftotal);
			probtotal += bfs[index] / bftotal;
			index++;
		}
	}

	if (verbose)
	{
		cout << YW << "  - naive_forloop_rlooper:: O(n*m/2) calculate prob: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
		t1 = timeNow();
	} // TIME

	index = 0;
	for (int n = 0; n < myseq.size(); n++)
	{
		if (n == 0) {
			printerfinal(outfile_naive_smalltable, fabeg, printSlowButCorrect, printQuickDirty, n, 0, n1[index], n2[index], Gbps[index], Gsigmas[index], Gs[index], bfs[index], bftotal);
			index++;
		}
		
		for (int m = 0; n + m < myseq.size() - 1; m++)
		{
			if (m > max_length - 1)
			{
				break;
			}
			
			if (m+2 == 1 || m+2 == 5 || m+2 == 10 || m+2 == 20 || m+2 == 50 || m+2 == 100 || m+2 == 150 || m+2 == 200 || m+2 == 250 || m+2 == 300 || m+2 == 500 || m+2 == 750 || m+2 == 1000 || m+2 == 2000) {
				printerfinal(outfile_naive_smalltable,  fabeg, printSlowButCorrect, printQuickDirty, n, m + 1, n1[index], n2[index], Gbps[index], Gsigmas[index], Gs[index], bfs[index], bftotal);
			}
			index++;
		}
	}
	if (printQuickDirty || printSlowButCorrect || printSlowButCorrect2)
	{
		outfile_naive_smalltable.close();
	}

	if (verbose)
	{
		cout << YW << "  - naive_forloop_rlooper:: O(n*m/2) PRINT _smalltable.tsv: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
		t1 = timeNow();
	} // TIME

	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<int> dd{probs.begin(), probs.end()};

	std::map<int, int> m;
	for (int n = 0; n < simN; ++n)
	{
		++m[dd(gen)];
	}

	if (verbose)
	{
		cout << YW << "  - naive_forloop_rlooper:: O(simN): calculate " << LCY << simN << YW << " simulation " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
		t1 = timeNow();
	} // TIME

	for (auto p : m)
	{
		int pindex = p.first + 1;
		int ptimes = p.second;
		for (int i = 0; i < p.second; i++)
		{
			if (ms[p.first + 1] >= minlength) {
				std::stringstream ss9; // in case we want to print
				ss9 << pindex;
				ss9 << "\t" << ns[p.first + 1] + fabeg;
				ss9 << "\t" << ms[p.first + 1];
				ss9 << "\t" << Gs[p.first + 1];
				ss9 << "\t" << bfs[p.first + 1];
				ss9 << "\t" << probs[p.first + 1];
				ss9 << "\n";
				outfile_naive_sim_reads << ss9.rdbuf();
			}
		}
	}
	outfile_naive_sim_reads.close();
	cout << "MS = " << ms[66] << "\n";
	if (verbose)
	{
		cout << YW << "  - naive_forloop_rlooper:: O(simN): PRINT sim_reads.tsv: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
		t1 = timeNow();
	} // TIME
	cout << myseq.size() << "\n";
	
	vector<int> mycount(myseq.size());
	for (auto p : m)
	{
		int pindex = p.first + 1;
		int ptimes = p.second;
		for (int i = 0; i < p.second; i++)
		{
			int beg = ns[p.first + 1];
			int len = ms[p.first + 1];
			int end = beg + len;
			if (ms[p.first + 1] >= minlength) {
	
				for (int j = beg; j < end; j++)
				{
					mycount[j]++;
				}
			}
		}
	}
	outfile_naive_sim_count << "#position\ttotal\ttotal_perc\n";
	for (int i = 0; i < mycount.size(); i++)
	{
		std::stringstream ss10;
		float myperc = 100 * mycount[i] / simN;
		int npos = i + fabeg;
		ss10 << npos << "\t" << mycount[i] << "\t" << setprecision(2) << myperc << "\n";
		outfile_naive_sim_count << ss10.rdbuf();
	}
	outfile_naive_sim_count.close();
	if (verbose)
	{
		cout << YW << "  - naive_forloop_rlooper:: O(n): PRINT sim_count.tsv " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
		t1 = timeNow();
	} // TIME

	int total_expected = myseq.size() * (myseq.size() + 1) / 2 - myseq.size() + 1;
	if (verbose)
	{
		cout << "------------------------\n";
		cout << "- Done!\n";
		cout << "  - naive_forloop_rlooper total rows = (n*m/2)\n";
		cout << "  - total    rows, where n = " << myseq.size() << ", m = min(seq_size = " << myseq.size() + 1 << ",max_length = " << max_length + 1 << "): " << index << "\n";
		cout << "  - expected rows, where n = " << myseq.size() << ", m = seq_size = " << myseq.size() + 1 << ": " << total_expected << "\n";
		cout << "  - difference: " << total_expected << " - " << index << " = " << total_expected - index << "\n";
		cout << "  - probtotal = " << probtotal << "\n";
		cout << "------------------------\n";
	}
}
//-------------------------}

void Simulation::simulation_A()
{ // some of this code might be migrated into new objects and functions in the future
	// initialize variables
	const char *LPR = "\033[1;95m";
	const char *LBU = "\033[1;94m";
	const char *LRD = "\033[1;31m";
	const char *LCY = "\033[1;36m";
	const char *YW = "\033[1;33m";
	const char *LGN = "\033[1;32m";
	const char *NN = "\033[1;0m";

	cout << LPR << "\nsimulation.cpp::" << "simulation_A" << NN << "\n";

	bool verbose = verbose_flag;
	int max_length = 2000;
	vector<Peak> external_structures;


	if (!infile.is_open())
	{
		throw UnexpectedClosedFileException("Simulation::simulation_A_D");
	}

	std::string outfilenames = "Outfiles:\n";
	if (dynamic_flag)
	{
		outfilenames += "\n- Dynamic outputs:\n" + string(LCY) +
						outfilename + "_dynamic_log.txt\n" +string(NN);
	}
	if (naive_flag)
	{
		outfilenames += "\n- Naive eq. outputs:\n" + string(LCY) +
						outfilename + "_naive_smalltable.tsv\n" +
						outfilename + "_naive_sim_reads.tsv\n" +
						outfilename + "_naive_sim_count.tsv\n" + string(NN);
	}
	if (orig_flag)
	{
		outfilenames += "\n- Orig eq. outputs:\n" + string(LCY) +
						outfilename + "_orig_bigtable.tsv\n" + string(NN);
	}

	ofstream outfile_naive_smalltable(outfilename + "_naive_smalltable.tsv", ios::out);
	ofstream outfile_naive_sim_reads(outfilename + "_naive_sim_reads.tsv", ios::out);
	ofstream outfile_naive_sim_count(outfilename + "_naive_sim_count.tsv", ios::out);
	ofstream outfile_dynamic_log(outfilename + "_dynamic_log.txt", ios::out);
	ofstream outfile_orig_bigtable(outfilename + "_orig_bigtable.tsv", ios::out);
	ofstream outfile_naive_bigtable(outfilename + "_naive_bigtable.tsv", ios::out);
	ofstream outfile_dynamic_bigtable(outfilename + "_dynamic_bigtable.tsv", ios::out);

	if (import_flag)
	{
	}

	bool eof = false;
	// Rloop_model* dynamic_model = static_cast<Rloop_model *>(models[0]);

	if (models.size() < 1)
	{
		// throw exception
	}
	if (verbose)
	{
		cout << "\n";
	}

	if (verbose)
	{
		cout << "-> is verbose!"
			 << "\n";
	}
	else
	{
		cout << "-> NOT verbose!"
			 << "\n";
	}
	if (dynamic_flag)
	{
		cout << "-> is dynamic!"
			 << "\n";
	}
	else
	{
		cout << "-> NOT dynamic!"
			 << "\n";
	}
	if (orig_flag)
	{
		cout << "-> is orig!"
			 << "\n";
	}
	else
	{
		cout << "-> NOT orig!"
			 << "\n";
	}

	while (!eof)
	{
		Gene *this_gene = new Gene();


		TimeVar t1 = timeNow();

		if (verbose)
		{
			cout << YW << "sim.cpp::time: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
			t1 = timeNow();
		}

		// 1. set min window size as minlength
		this_gene->windower.set_min_window_size(minlength);

		// 2. read gene
		eof = this_gene->read_gene(infile, genebeg, geneend);

		cout << "processing gene: " << this_gene->getName() << "..."
			 << "\n\n";

		// VERBOSE: PRINT SEQUENCE
		// if (verbose)
		// {
			// cout << "Sequence:\n"
				//  << LCY << join(this_gene->getSequence()) << NN << "\n";
		// }
		// VERBOSE: TIME
		if (verbose)
		{
			cout << YW << "sim.cpp::set gene io " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
			t1 = timeNow();
		} // TIME

		// compute structures using models
		if (auto_domain_size)
		{
			static_cast<Rloop_model *>(models[0])->setN(this_gene->get_length()); // need to compute this from the actual sequence.
		}
		if (this_gene->getPosition().strand == "-")
		{
			this_gene->reverse_complement_sequence();
		}
		if (complement_flag)
		{
			this_gene->complement_sequence();
		}
		if (reverse_flag)
		{
			this_gene->reverse_sequence();
		}

		// VERBOSE: TIME
		if (verbose)
		{
			cout << YW << "sim.cpp::set flags: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
			t1 = timeNow();
		} // TIME

		// ORIGINAL EQ RLOOPER
		// NAIVE RLOOPER
		// naive/just simple forloop no constructor etc. 2310bp = 1s without printing, 2-3s with printing

		if (naive_flag)
		{
			if (verbose)
			{
				cout << LBU << "\n##############################\n";
				cout << LBU << "sim.cpp:: Running naive_forloop_rlooper()!" << NN << "\n";
			}
			double my_K = models[0]->getK2();		  // 2200 * 0.0019858775 * 310 / 1500; //models[0]->getK();
			double my_sigma = models[0]->getSigma2(); //-0.07; //models[0]->getSigma();
			double my_a = models[0]->geta2();		  // 10; //models[0]->geta();
			const int simN = 10000;
			const int nick = models[0]->getnick2();
			const int nicklen = models[0]->getnicklen2();
			const int selffoldlen = models[0]->getselffoldlen2();
			cout << "K           : " << models[0]->getK2() << "\n";
			cout << "sigma       : " << models[0]->getSigma2() << "\n";
			cout << "a           : " << models[0]->geta2() << "\n";
			cout << "nick        : " << models[0]->getnick2() << "\n";
			cout << "nicklen     : " << models[0]->getnicklen2() << "\n";
			cout << "selffoldlen : " << models[0]->getselffoldlen2() << "\n";
			naive_forloop_rlooper(this_gene->getSequence(), this_gene->fabeg, outfile_naive_bigtable, outfile_naive_smalltable, outfile_naive_sim_reads, outfile_naive_sim_count, my_K, my_a, my_sigma, simN, nick, nicklen, selffoldlen, verbose, minlength);

			if (verbose)
			{
				cout << LBU << "##############################\n";
			}
			if (verbose)
			{
				cout << LBU << "sim.cpp:: TOTAL O(m*n) naive_forloop_rlooper: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n\n";
				t1 = timeNow();
			} // TIME
		}


		// DYNAMIC RLOOPER
		// dynamic model
		// WORK IN PROGRESS! Does not affect anything above this
		if (dynamic_flag)
		{
			srand(seed);
						// TIME
			if (verbose)
			{
				cout << LPR << "\n##############################" << NN << "\n";
				cout << LPR << "sim.cpp:: Running dynamic_rlooper()!" << NN << "\n";
				t1 = timeNow();
			}
			bool myprint = false;
			int nick = 0;
			vector<vector<double>> df_bfs;
			vector<vector<double>> df_bftotals;
			vector<vector<double>> df_probs;

			if (verbose)
			{
				cout << "- " << YW << "sim.cpp:: TOTAL O(n*m) precalculating bfs/probs/bftoatls: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
				t1 = timeNow();
			}
			cout << ">>" << YW << "SEED = " << seed << NN << "\n";


			// TIME
			if (verbose)
			{
				cout << LPR << "\n##############################" << NN << "\n";
				cout << LPR << "sim.cpp:: Running dynamic_rlooper()!" << NN << "\n";
				t1 = timeNow();
			}
			Rloop_model *Dmodel = static_cast<Rloop_model *>(models[0]);
			Dmodel->seqD = this_gene->getSequence();
			cout << "\n\n";
			cout << LCY << join(Dmodel->seqD) << NN << "\n";
			cout << "- amb_link_diff; " << Dmodel->ambient_alpha << "\n";
			cout << "- window size: " << LGN << dynamic_window_size << NN << "\n";
			cout << "- seq size: " << LGN << this_gene->getSequence().size() << NN << "\n";
			cout << "- A:  " << LGN << Dmodel->getA() << NN << "\n\n";

			Dmodel->dynamic_setWindow_size(dynamic_window_size);
			vector<int> maxlengths;
	
			for (int currSim = 0; currSim < Dmodel->dynamic_getTotalSim(); currSim++)
			{
				Dmodel->dynamic_setCurrSim(currSim);
				if (currSim % 20 == 0)
				{
					cout << "- Done " << LGN << currSim << NN << " / " << LCY << Dmodel->dynamic_getTotalSim() << NN << "\n";
				}
				stringstream ss_dynamic_log;
				Dmodel->dynamic_reset_model();

				// calculate the rloop independent superhelical conditions at the current polymerase position in the current rloop state.
				int iter = -1;
				int begpos = -1;
				int endpos = -1;
				double urninit = 0;
				double urnelon = 0;
				double probinit = 0;
				double probelon = 0;
				bool mybreak = false;

				// 1. Ambient sigma, --sigma (e.g. 14% below) THESE NEVER CHANGE
				double amb_sigma = Dmodel->ambient_sigma;
				double amb_alpha = Dmodel->ambient_alpha; // ambient alpha NEVER changes! same as Dmodel->getN() * Dmodel->getA() * Dmodel->dynamic_getSigma(); // 1500 * 1/10.5 * -14%  = -20 (ambient alpha or delta Lk) // dynamic_compute_res_lk();//->getAlpha();

				// for (my $iter = 0; $iter < $plasmid_length; $iter ++)
				while (Dmodel->dynamic_getCurrentPos() < this_gene->getSequence().size())
				{ // until end of sequence
					iter++;
					Dmodel->dynamic_setIter(iter);

				
					// 2. Tx sigma is always 7% over 1500bp
					double tx_window = min(Dmodel->dynamic_getCurrentPos(), Dmodel->getN());           // currpos 15: min(15,1500) = 15; currpos 1501: min(1501,1500) = 1500
					double tx_sigma = tx_window / Dmodel->getN() * Dmodel->dynamic_getTx_Amb_Sigma(); // 15bp / 1500bp * -7% = -0.07% (-7% over 15bp "diluted" in 1500bp)
					double tx_alpha = Dmodel->getN() * Dmodel->getA() * tx_sigma;					  // 1500 * 1/10.5 * -0.07% = -0.1

					Dmodel->dynamic_setTx_Sigma(tx_sigma);
					Dmodel->dynamic_setTx_Alpha(tx_alpha);

					// 3. Rloop
					double rloop_alpha_local = (4*pow(pi,2)*Dmodel->getC()) / (4*pow(pi,2)*Dmodel->getC()+Dmodel->getK()*Dmodel->dynamic_getN_rloop_bases()) * 
					(Dmodel->dynamic_getAlphaTotal()+Dmodel->dynamic_getN_rloop_bases()*Dmodel->getA());
					double rloop_sigma = rloop_alpha_local / (Dmodel->getN() * Dmodel->getA());
					double rloop_alpha = Dmodel->getN() * Dmodel->getA() * rloop_sigma;
					// return ((4*pow(pi,2)*C) / (4*pow(pi,2)*C+k*n_rloop_bases)) * (alpha_total+n_rloop_bases*A);


					// 4. Alpha total
					double sigma_total = amb_sigma + tx_sigma; // plasmid sigma 14% + tx sigma -0.07% = 14.07% per 15bp
					double alpha_total = amb_alpha + tx_alpha;	// plasmid alpha -20 + tx alpha -0.1 = -20.1

					Dmodel->dynamic_setSigmaTotal(sigma_total);
					Dmodel->dynamic_setAlphaTotal(alpha_total);
					Dmodel->dynamic_setAlpha(Dmodel->dynamic_getAlphaTotal());

					int currpos = Dmodel->dynamic_getCurrentPos(); // for printing
					if (this_gene->fabeg != -1) {
						currpos += this_gene->fabeg;
					}

				    std::map<std::string, std::vector<std::string>> data;

					// vector<string> header = {"currSim","iter","currpos",};
					vector<string> myheader = {
						"CurrSim",
						"Iter",
						"N     ",
						"CurrentPos  ",
						"SigmaAmb    ",
						"SigmaTxAmb  ",
						"SigmaTx     ",
						"RloopSigma  ",
						"SigmaTotal  ",
						"AlphaAmb    ",
						// "AlphaTxAmbInit_" + to_string(Dmodel->dynamic_getWindow_size()),
						// "AlphaTxAmbElon_" + to_string(Dmodel->dynamic_getElongation_step_size()),
						"AlphaTx     ",
						"AlphaTotal  ",
						// "RloopAlphaL ",
						"RloopAlpha  ",
						"AlphaCurr   ",
						"Rloop_bases ",
						"Rloop_state ",
						"Gjunc       ",
						"Gbp         ",
						"Gsigma      ",
						"Gtot        ",
						"G0          ",
						"Gtot_bf     ",
						"G0_bf       ",
						"part_func   ",
						"RNG         ",
						"Prob        "
						};
//		vector<string> valuesToAdd = {"2_elon",to_string(Gjunc),to_string(Gbps),to_string(Gsigma),to_string(Gtot),to_string(G0),to_string(Gtot_bf),to_string(G0_bf),to_string(partition_function),to_string(urn),to_string(prob_all_rloop_elon)};

					if (Dmodel->ambient_sigma != Dmodel->getSigma()) {
						cout << "ERROR! Dmodel->ambient_sigma " << Dmodel->ambient_sigma << " != Dmodel->getSigma() or sigma " << Dmodel->getSigma() << "!\n";
						exit(1);
					}
					// cout <<	to_string(Dmodel->dynamic_getTx_Sigma()) << "=" << to_string(tx_window) << "/" << to_string(Dmodel->getN()) << "*" << to_string(Dmodel->dynamic_getTx_Amb_Sigma()) << " -> ";
					// cout <<	Dmodel->dynamic_getTx_Sigma() << "=" << tx_window << "/" << Dmodel->getN() << "*" << Dmodel->dynamic_getTx_Amb_Sigma() << " -> ";
					// cout << tx_window / Dmodel->getN() *Dmodel->dynamic_getTx_Amb_Sigma() << "\n";
					
					vector<string> myvalues = {
						to_string(Dmodel->dynamic_getCurrSim()),
						to_string(Dmodel->dynamic_getIter()),
						to_string(Dmodel->getN()),
						to_string(Dmodel->dynamic_getCurrentPos()),
						to_string(Dmodel->ambient_sigma),
						to_string(Dmodel->dynamic_getTx_Amb_Sigma()),
						to_string(Dmodel->dynamic_getTx_Sigma()), // + "=" + to_string(tx_window) + "/" + to_string(Dmodel->getN()) + "*" + to_string(Dmodel->dynamic_getTx_Amb_Sigma()), // 15bp / 1500bp * -7% = -0.07% (-7% over 15bp "diluted" in 1500bp)
						to_string(rloop_sigma),
						to_string(Dmodel->dynamic_getSigmaTotal()),
						to_string(Dmodel->ambient_alpha),
						// to_string(Dmodel->dynamic_getTx_Amb_Alpha()),
						// to_string(Dmodel->dynamic_getTx_Amb_Alpha()),
						to_string(Dmodel->dynamic_getTx_Alpha()),
						to_string(Dmodel->dynamic_getAlphaTotal()),
						// to_string(rloop_alpha_local),
						to_string(rloop_alpha),
						to_string(Dmodel->dynamic_getAlpha()),
						to_string(Dmodel->dynamic_getN_rloop_bases())
					};
					// header
					// if (currSim == 0 && iter == 0)
					// {
					// // for (const auto& pair : data) {
					// // 	ss_dynamic_log2 << pair.first << ",";
					// // }
					// 	// ss_dynamic_log2 << "currSim,iter,currpos,N,bg_sigma,added_sigma,getAlpha,ambient_alpha,added_alpha,total_alpha,res_lk,n_rloop,rloop_state,G0,Gtot,G0_bf,Gtot_bf,random,exp\n";
					// }
					// for (const auto& pair : data) {
					// 	for (const std::string& value : pair.second) {
					// 		ss_dynamic_log2 << value << ",";
					// 	}
					// 	// ss_dynamic_log2 << "\n";
					// }

					// if nrloop is more than 0, then set alpha to be N*A*bg_sigma
					// if (Dmodel->dynamic_getN_rloop_bases() > 0) { // to avoid /0 errors
					// 	Dmodel->setAlpha(Dmodel->dynamic_compute_res_lk()); // returns residual linking difference ?????
					// 	//Dmodel->setAlpha(Dmodel->getN()*Dmodel->getA()*Dmodel->getSigma());
					// }
					// else {// otherwise set alpha as alpha total
					// 	Dmodel->setAlpha(Dmodel->dynamic_getAlphaTotal());
					// 	// setAlpha(alpha_total)
					// }

					// Dmodel->dynamic_print_topological_state();
					// ss_dynamic_log2 << std::setprecision(3) << currSim << ","
					// 				<< iter << ","
					// 				<< Dmodel->dynamic_getCurrentPos() << ","
					// 				<< Dmodel->getN() << ","
					// 				<< Dmodel->getSigma() << ","
					// 				<< Dmodel->dynamic_getTx_Amb_Sigma() << ","
					// 				<< Dmodel->getAlpha() << ","
					// 				<< amb_alpha << ","
					// 				<< tx_alpha << ","
					// 				<< Dmodel->dynamic_getAlphaTotal() << ","
					// 				<< Dmodel->dynamic_compute_res_lk() << ","
					// 				//<< Dmodel->ground_state_factor() << ","
					// 				//<< gs << ","
					// 				<< Dmodel->dynamic_getN_rloop_bases();

					if (Dmodel->in_rloop == false) { // if in the initiation regime
						begpos = -1;
						Dmodel->dynamic_step_forward_initiation(df_bfs, df_bftotals, df_probs, verbose, seed, begpos, urninit, probinit, myvalues);//, ss_dynamic_log2);
					}
					else {
						bool mycontinue = Dmodel->dynamic_step_forward_elongation(df_bfs, df_bftotals, df_probs, verbose, seed, endpos, urnelon, probinit, myvalues);//, ss_dynamic_log2);

						if (mycontinue == false) {
							int offset = 0;
							if (this_gene->fabeg != -1) {
								offset = this_gene->fabeg;
							}
							std::stringstream ssdyn; // in case we want to print
							if (begpos == -1) {
								ssdyn << infilename << "\t" << offset + begpos << "\t" << offset + endpos << "\t"
									  << "NOPK_" << currSim << "\t" << endpos - begpos << "\t" << this_gene->getPosition().strand << "\t" << urninit << "\t" << probinit << "\t" << urnelon << "\t" << probelon << "\n";
							}
							else {
								ssdyn << infilename << "\t" << offset + begpos << "\t" << offset + endpos << "\t"
									  << "PEAK_" << currSim << "\t" << endpos - begpos << "\t" << this_gene->getPosition().strand << "\t" << urninit << "\t" << probinit << "\t" << urnelon << "\t" << probelon << "\n";
							}
							outfile_dynamic_bigtable << ssdyn.rdbuf();
							ssdyn.str("");
							begpos = -1;
							endpos = -1;

							// mybreak = true; //only allows one R-loop per simulation
							//continue; // allows any number of R-loops per simulation
						}
					}

					// for (size_t i = 0; i < header.size(); ++i) {
					// 	data[header[i]].push_back(myvalues[i]);
					// }

					std::stringstream ss_dynamic_log2;
					int maxvallen = 15;
					int tempiter = 0;
					if (currSim == 0 && iter == 0)
					{
						tempiter = 0;
						for (int headerInd = 0; headerInd < myheader.size(); headerInd ++) {
						// for (const auto& pair : data) {
							string currheader = myheader[headerInd];
							if (tempiter != 0) {
								ss_dynamic_log2 << " ";
							}
							tempiter ++;
							int currvallen = currheader.size(); //pair.first.size();
							int addlen = currvallen; //max(maxvallen - currvallen + 0, 0);
							maxlengths.push_back(addlen);
							ss_dynamic_log2 << currheader; //pair.first;
						}
						ss_dynamic_log2 << "\n";
						// ss_dynamic_log2 << "currSim,iter,currpos,N,bg_sigma,added_sigma,getAlpha,ambient_alpha,added_alpha,total_alpha,res_lk,n_rloop,rloop_state,G0,Gtot,G0_bf,Gtot_bf,random,exp\n";
					}
					// for (const auto& pair : data) {
					// 	for (const std::string& currvalue : pair.second) {
					tempiter = 0;
					for (int myvaluesInd = 0; myvaluesInd < myheader.size(); myvaluesInd ++) {
						string currvalue = myvalues[myvaluesInd];
						int currmaxlengths = maxlengths[myvaluesInd];
						if (tempiter != 0) {
							ss_dynamic_log2 << " ";
						}
						tempiter ++;
						int currvallen = currvalue.size();
						int addlen = 0;
						if (currmaxlengths > currvallen) {
							addlen = currmaxlengths - currvallen;
						 } //max(maxvallen - currvallen + 0, 0);
						ss_dynamic_log2 << currvalue;
						for (int templen = 0; templen < addlen; templen++) {
							ss_dynamic_log2 << " ";
						}
						// }
					}
					ss_dynamic_log2 << "\n";

					outfile_dynamic_log << ss_dynamic_log2.rdbuf();
					

					// std::map<std::string, int> column_widths;
					// for (const auto& pair : data) {
					// 	int max_width = pair.first.size();
					// 	for (const std::string& value : pair.second) {
					// 		max_width = std::max(max_width, static_cast<int>(value.size()));
					// 	}
					// 	column_widths[pair.first] = max_width;
					// }

					// if (Dmodel->dynamic_getCurrSim() == 0 && Dmodel->dynamic_getIter() == 0)
					// {
					// 	for (size_t i = 0; i < header.size(); ++i) {
					// 		ss_dynamic_log2 << std::left << std::setw(column_widths[header[i]]) << header[i] << ",";
					// 	}
					// }
					// for (size_t i = 0; i < header.size(); ++i) {
					// 	for (const std::string& value : data[header[i]]) {
					// 		ss_dynamic_log2 << value << ",";
					// 	}
					// 	// ss_dynamic_log2 << std::endl;
					// }							

					// outfile_dynamic_log << ss_dynamic_log2.rdbuf();
						
					if (mybreak == true) {
						break;
					}
				}

				std::stringstream ssdyn; // in case we want to print
				if (begpos == -1) {
					ssdyn << infilename << "\t" << begpos << "\t" << endpos << "\t"
						  << "NOPK_" << currSim << "\t" << endpos - begpos << "\t" << this_gene->getPosition().strand << "\t" << urninit << "\t" << probinit << "\t" << urnelon << "\t" << probelon << "\n";
				}
				else {
					ssdyn << infilename << "\t" << begpos << "\t" << endpos << "\t"
						  << "PEAK_" << currSim << "\t" << endpos - begpos << "\t" << this_gene->getPosition().strand << "\t" << urninit << "\t" << probinit << "\t" << urnelon << "\t" << probelon << "\n";
				}
				outfile_dynamic_bigtable << ssdyn.rdbuf();
				ssdyn.str("");
				begpos = -1;
				endpos = -1;

				// cout << Dmodel->dynamic_write_buffer.rdbuf();
			}

			// TIME
			if (verbose)
			{
				float mydur = ((duration(timeNow() - t1) / 1e9) / Dmodel->dynamic_getTotalSim());
				cout << LPR << "sim.cpp:: TOTAL O(simN*m*n) dynamic_rlooper: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s (" << mydur << "s/simN)"
					 << "\n\n";
				cout << LPR << "\n##############################" << NN << "\n";
				t1 = timeNow();
			} // TIME
		}

		// modified original equilibrium model, modified original (fixed the 3 nested for loops when calculating bp_energy -> N^3 into N^2). 2310bp = 1s without printing
		if (orig_flag)
		{
			if (verbose)
			{
				cout << LCY << "\n##############################\n";
				cout << LCY << "sim.cpp:: Running original_rlooper()!" << NN << "\n";
			}
			if (circular_flag)
			{
				this_gene->compute_structures_circular(*models[0], verbose);
			}
			else
			{
				this_gene->compute_structures(*models[0], verbose);
			}
			if (verbose)
			{
				cout << "- " << YW << "sim.cpp:: TOTAL O(n*m) compute_structures: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
				t1 = timeNow();
			} // TIME

			if (import_flag)
			{ // need to eventually match imported structures with their associated genes / plasmids
				cout << "- "
					 << "importing external structures from " << importfilename << "..."
					 << "\n";
				external_structures = import_external_structures(importfilename, *models[0]);
				this_gene->compute_structures_external(external_structures, *models[0], verbose);
				cout << "- "
					 << "complete!\n";
			}
			cout << "\n";
			cout << "- Total index  : " << this_gene->getRloopStructures().size() << "\n";
			cout << "- Sequence size: " << this_gene->getSequence().size() << "\n";

			// ensemble analysis, free energies and boltzmann factors have already been computed in compute_structures
			// compute partition function
			long double partition_function = 0; // partition function is the denominator of the p(G(s)) = boltzman(G(s)) / sum(boltzman(G(s)))
			long double sanity_check = 0;
			int iter0 = 0;
			int n = 1;
			int m = 2;
			for (vector<Structure>::iterator it = this_gene->getRloopStructures().begin(); it < this_gene->getRloopStructures().end(); ++it)
			{
				partition_function += it->boltzmann_factor; // sum all the boltzman(G(s))
				if (m == 2)
				{
					partition_function += models[0]->ground_state_factor();
					iter0++;
				}
				iter0++;
				if (n + m == 1 + this_gene->getSequence().size())
				{
					n++;
					m = 2;
				}
				else
				{
					m++;
				}
			}
			partition_function += models[0]->ground_state_factor(); // state N (ACTG-> state at G, where there will be 100% no rloop/ground state

			cout << "- "
				 << "ground_state: " << models[0]->ground_state_factor() << "\n"; // cout the m = 0 ground state boltzman too
			// partition_function += models[0]->ground_state_factor(); //
			// compute boltzmann weights and store in the structures
			iter0 = 1;
			n = 1;
			m = 2;
			for (vector<Structure>::iterator it = this_gene->getRloopStructures().begin();
				 it < this_gene->getRloopStructures().end(); ++it)
			{
				stringstream ssorig;
				it->probability = it->boltzmann_factor / partition_function;
				sanity_check += it->boltzmann_factor / partition_function;
				// int mprint = m;
				// if (m == 1) {
				// 	mprint = 0;
				// }
				if (m == 2)
				{
					ssorig << n << '\t' << 0 << '\t' << 0 << '\t' << models[0]->ground_state_energy() << '\t' << models[0]->ground_state_energy() << '\t' << models[0]->ground_state_factor() << '\t' << models[0]->ground_state_factor() / partition_function << "\n";
					iter0++;
				}
				ssorig << n << '\t' << m << '\t' << it->bp_energy << '\t' << it->gsigma << '\t' << it->free_energy << '\t' << it->boltzmann_factor << '\t' << it->probability << "\n";
				iter0++;
				if (n + m == 1 + this_gene->getSequence().size())
				{
					n++;
					m = 2;
				}
				else
				{
					m++;
				}
				outfile_orig_bigtable << ssorig.rdbuf();
				ssorig.str(std::string());
			}

			stringstream ssorig;
			ssorig << n << '\t' << 0 << '\t' << 0 << '\t' << models[0]->ground_state_energy() << '\t' << models[0]->ground_state_energy() << '\t' << models[0]->ground_state_factor() << '\t' << models[0]->ground_state_factor() / partition_function << "\n";
			outfile_orig_bigtable << ssorig.rdbuf();
			ssorig.str(std::string());

			sanity_check += models[0]->ground_state_factor() / partition_function;
			cout << "- "
				 << "P(ground state): " << models[0]->ground_state_factor() / partition_function << "\n";
			// if (fabs(1-sanity_check) > .00001){
			// throw SimulationException("Ensemble probability sum != 1"); //this throw is uncaught
			// }
			if (verbose)
			{
				cout << "- " << YW << "sim.cpp:: Main calculation done! " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
				t1 = timeNow();
			} // TIME
			// std::sort(this_gene->getRloopStructures().begin(), this_gene->getRloopStructures().end());
			// compute signals and output .wig tracks
			vector<double> *signal = NULL, *signal2 = NULL, *signal3 = NULL, *signal4 = NULL;
			vector<Loci> peaks;

			if (verbose)
			{
				cout << "- " << YW << "sim.cpp:: sorting this_gene: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
				t1 = timeNow();
			} // TIME

			t1 = timeNow();
			compute_signal_bpprobs(*this_gene, signal);
			if (verbose)
			{
				cout << "- " << YW << "- sim.cpp:: TOTAL O(n*m) compute_signal_bpprobs: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
				t1 = timeNow();
			} // TIME
			if (import_flag)
			{
				compute_signal_extbpprobs(*this_gene, signal4);
				if (verbose)
				{
					cout << "- " << YW << "- sim.cpp:: compute_signal_extbpporbs: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
					t1 = timeNow();
				} // TIME
			}
			if (average_g)
			{
				compute_signal_average_G(*this_gene, signal2);
				if (verbose)
				{
					cout << "- " << YW << "- sim.cpp:: compute_signal_avgG: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
					t1 = timeNow();
				} // TIME
			}
			compute_signal_mfe(*this_gene, signal3);
			if (verbose)
			{
				cout << "- " << YW << "- sim.cpp:: compute_signal_mfe: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
				t1 = timeNow();
			} // TIME
			// write signals
			// write_wigfile(outfile1, this_gene, signal);
			if (verbose)
			{
				cout << "- " << YW << "- sim.cpp:: write_wigfile bpprobs: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
				t1 = timeNow();
			} // TIME

			if (import_flag)
			{
				// write_wigfile(outfile6, this_gene, signal4);
				if (verbose)
				{
					cout << "- " << YW << "- sim.cpp:: write_wigfile extbpprobs: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
					t1 = timeNow();
				} // TIME
			}
			if (average_g)
			{
				// write_wigfile(outfile2, this_gene, signal2);
				if (verbose)
				{
					cout << "- " << YW << "- sim.cpp:: write_wigfile avgG.wig: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
					t1 = timeNow();
				} // TIME
			}
			// write_wigfile(outfile3, this_gene, signal3);
			if (verbose)
			{
				cout << "- " << YW << "- sim.cpp:: write_wigfile mfe.wig: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
				t1 = timeNow();
			} // TIME
			// call peaks and write results to .bed files
			if (printbedfile)
			{

				call_peaks_absolute_threshold(*this_gene, *signal, peaks); // possible null pointer exception generated here
				if (verbose)
				{
					cout << "- " << YW << "- sim.cpp:: call_peaks_abs_thres with bpprobs.wig: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
					t1 = timeNow();
				} // TIME
				// write to bedfile
				// write_bedfile(outfile4, this_gene, peaks);
				if (verbose)
				{
					cout << "- " << YW << "- sim.cpp:: write_bedfile bpprobs.bed: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
					t1 = timeNow();
				} // TIME
				peaks.clear();
				call_peaks_absolute_threshold(*this_gene, *signal3, peaks); // possible null pointer exception generated here
				// write to bedfile
				if (verbose)
				{
					cout << "- " << YW << "- sim.cpp:: call_peaks_abs_thres with mfe.wig: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
					t1 = timeNow();
				} // TIME
				// write_bedfile(outfile5, this_gene, peaks);
				if (verbose)
				{
					cout << "- " << YW << "- sim.cpp:: write_bedfile mfe.bed: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
					t1 = timeNow();
				} // TIME
			}
			cout << "- complete!\n";

			if (verbose)
			{
				cout << "- " << YW << "- sim.cpp:: Done with main orig_rlooper: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
				t1 = timeNow();
			} // TIME
			break;
			// output residuals if the option is selected
			if (residuals)
			{
				int gq_length = 0;
				double ensemble_residual_twist = 0, ensemble_current_Gsigma = 0,
					   ensemble_wrapping_absorption = 0, ensemble_strand_separation_absorption = 0;
				Rloop_model *temp = (Rloop_model *)models[0];
				this_gene->compute_residuals(*models[0], verbose);
				if (import_flag)
				{
					gq_length = external_structures[0].position.get_length(); // only works for one imported structure for now
				}
				for (vector<Structure>::iterator it = this_gene->getRloopStructures().begin(); it < this_gene->getRloopStructures().end(); ++it)
				{
					ensemble_residual_twist += it->residual_twist * it->probability;
					ensemble_current_Gsigma += it->current_Gsigma * it->probability;
					ensemble_wrapping_absorption += (temp->getAlpha() - it->current_Gsigma + (it->position.get_length() - gq_length) * temp->getA()) * it->probability;
					ensemble_strand_separation_absorption -= ((it->position.get_length() - gq_length) * temp->getA()) * it->probability;
				}
				// consider the ground state as well
				double twist = 0, writhe = 0;
				models[0]->ground_state_residuals(twist, writhe);
				ensemble_residual_twist += twist * (models[0]->ground_state_factor() / partition_function);
				ensemble_current_Gsigma += writhe * (models[0]->ground_state_factor() / partition_function);
				cout << "\t- "
					 << "ensemble_residual_twist: " << ensemble_residual_twist << "\n";
				cout << "\t- "
					 << "ensemble_current_Gsigma: " << ensemble_current_Gsigma << "\n";
				// convert linking difference to superhelicity
				cout << "\t- "
					 << "ensemble_residual_superhelicity: " << ensemble_current_Gsigma / (temp->getN() * temp->getA()) << "\n";
				cout << "\t- "
					 << "ensemble_wrapping_absorption: " << ensemble_wrapping_absorption << "\n";
				cout << "\t- "
					 << "ensemble_strand_separation_absorption: " << ensemble_strand_separation_absorption << "\n";

				if (verbose)
				{
					cout << "- " << YW << "- sim.cpp:: residuals: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
					t1 = timeNow();
				} // TIME
			}
			if (top > 0)
			{
				// sort top N structures into a new vector
				std::sort(this_gene->getRloopStructures().begin(), this_gene->getRloopStructures().end());
				Rloop_model *temp = (Rloop_model *)models[0];
				// output structures to .bed file
				for (int i = 0; i < top; i++)
				{
					// if the sequence has been reversed, output the reversed coordinates for the top structures
					if (this_gene->getPosition().strand == "-")
					{
						cout << "\t- " << this_gene->getSequence().size() - this_gene->getRloopStructures()[i].position.start_pos << ' '
							 << this_gene->getSequence().size() -
									this_gene->getRloopStructures()[i].position.end_pos
							 << ' ';
					}
					else
					{ // gene is on + strand
						cout << "\t- " << this_gene->getRloopStructures()[i].position.start_pos << ' '
							 << this_gene->getRloopStructures()[i].position.end_pos << ' ';
					}
					cout << "\t- " << this_gene->getRloopStructures()[i].free_energy << ' '
						 << this_gene->getRloopStructures()[i].probability << ' '
						 << this_gene->getRloopStructures()[i].residual_twist << ' '
						 << this_gene->getRloopStructures()[i].current_Gsigma << ' '
						 << this_gene->getRloopStructures()[i].current_Gsigma / (temp->getN() * temp->getA()) << "\n";
				}

				if (verbose)
				{
					cout << "- " << YW << "- sim.cpp:: top: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
					t1 = timeNow();
				} // TIME
			}
			if (dump)
			{
				this_gene->dump_structures(outfilename);
				if (verbose)
				{
					cout << "- " << YW << "- sim.cpp:: dump: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n";
					t1 = timeNow();
				} // TIME
			}
			delete signal;

			// clear_sequence the sequence data from the gene to save memory
			this_gene->clear_sequence();
			this_gene->clear_structures();
			// store the gene in the genes vector
			genes.push_back(this_gene);

			if (verbose)
			{
				cout << LCY << "##############################\n";
			}
			if (verbose)
			{
				cout << LCY << "sim.cpp:: TOTAL O(m*n) original_rlooper: " << LGN << duration(timeNow() - t1) / 1e9 << NN << " s\n\n";
				t1 = timeNow();
			} // TIME
		}
	}
	outfile_orig_bigtable.close();
	outfile_dynamic_bigtable.close();
	outfile_naive_bigtable.close();
	outfile_naive_smalltable.close();
	outfile_naive_sim_count.close();
	outfile_naive_sim_reads.close();

	cout << "\n"
		 << outfilenames << NN << "\n\n";
}

// computes P(R-Loop) for the provided supercoiling value
void Simulation::simulation_B(float superhelicity, ofstream &outfile)
{
	bool verbose = verbose_flag;
	vector<Peak> external_structures;
	if (!infile.is_open())
	{
		throw UnexpectedClosedFileException("Simulation::simulation_B");
	}
	if (models.size() < 1)
	{
		// throw exception
	}
	float p_rloop = 0;
	Gene *this_gene;
	if (!genes.size())
	{
		this_gene = new Gene();
		this_gene->read_gene(infile, genebeg, geneend);
		this_gene->windower.set_min_window_size(minlength);
		if (this_gene->getPosition().strand == "-")
		{
			this_gene->reverse_complement_sequence();
		}
		if (complement_flag)
		{
			this_gene->complement_sequence();
		}
		if (reverse_flag)
		{
			this_gene->reverse_sequence();
		}
		// this_gene->complement_sequence();
		// this_gene->reverse_sequence();
		genes.push_back(this_gene);
	}
	else
	{
		this_gene = genes[0];
		this_gene->clear_structures();
	}
	models[0]->set_superhelicity(superhelicity);		  // set the superhelicity in the model to the provided value
	this_gene->clear_structures();						  // saves memory
	this_gene->compute_structures(*(models[0]), verbose); // from gene.h::compute_structures()
	if (import_flag)
	{ // need to eventually match imported structures with their associated genes / plasmids
		cout << "\t- "
			 << "importing external structures from " << importfilename << "..."
			 << "\n";
		external_structures = import_external_structures(importfilename, *models[0]);
		this_gene->compute_structures_external(external_structures, *models[0], verbose);
		cout << "\t- "
			 << "complete!"
			 << "\n";
	}
	// determine P(ground state)
	long double partition_function = 0;
	long double ground_state_factor = 0;
	int index = this_gene->getRloopStructures().size();
	int count = 0;
	for (vector<Structure>::iterator it = this_gene->getRloopStructures().begin();
		 it < this_gene->getRloopStructures().end(); ++it)
	{
		partition_function += it->boltzmann_factor;
		count++;
	}
	ground_state_factor = models[0]->ground_state_factor();
	partition_function += ground_state_factor;
	// determine P(R-Loop) as 1-P(ground state)
	p_rloop = 1 - (ground_state_factor / partition_function);
	// display result
	outfile << superhelicity << ' ' << p_rloop << "\n";
}

void Simulation::simulation_C(float superhelicity, ofstream &outfile)
{
	bool verbose = verbose_flag;
	if (!infile.is_open())
	{
		throw UnexpectedClosedFileException("Simulation::simulation_C");
	}
	if (models.size() < 1)
	{
		// throw exception
	}
	Gene *this_gene;
	if (!genes.size())
	{
		this_gene = new Gene();
		this_gene->read_gene(infile, genebeg, geneend);
		this_gene->windower.set_min_window_size(minlength);
		genes.push_back(this_gene);
	}
	else
	{
		this_gene = genes[0];
		this_gene->clear_structures();
	}
	if (this_gene->getPosition().strand == "-")
	{
		this_gene->reverse_complement_sequence();
	}
	if (complement_flag)
	{
		this_gene->complement_sequence();
	}
	if (reverse_flag)
	{
		this_gene->reverse_sequence();
	}
	models[0]->set_superhelicity(superhelicity); // set the superhelicity in the model to the provided value
	if (circular_flag)
	{
		this_gene->compute_structures_circular(*models[0], verbose);
	}
	else
	{
		this_gene->compute_structures(*models[0], verbose);
	}

	// determine P(ground state)
	long double partition_function = 0;
	long double ground_state_factor = 0;
	unsigned long index = this_gene->getRloopStructures().size();
	int count = 0;
	for (vector<Structure>::iterator it = this_gene->getRloopStructures().begin();
		 it < this_gene->getRloopStructures().end(); ++it)
	{

		partition_function += it->boltzmann_factor;
		count++;
	}
	ground_state_factor = models[0]->ground_state_factor();
	partition_function += ground_state_factor;
	// determine expected length at the given superhelicity value
	double expected_length = 0, n = 0;
	for (vector<Structure>::iterator it = this_gene->getRloopStructures().begin();
		 it < this_gene->getRloopStructures().end(); ++it)
	{
		if (top > 0 && it->position.get_length() > top)
		{
			expected_length += (it->boltzmann_factor / partition_function) * it->position.get_length();
		}
		else if (top == 0)
		{
			expected_length += (it->boltzmann_factor / partition_function) * it->position.get_length();
		}
	}
	double var = 0;
	n = 0;
	// compute and report weighted variance
	for (vector<Structure>::iterator it = this_gene->getRloopStructures().begin();
		 it < this_gene->getRloopStructures().end(); ++it)
	{
		if (top > 0)
		{
			// do nothing
		}
		else if (top == 0)
		{
			// build weighted variance for the whole ensemble
			var += (it->boltzmann_factor / partition_function) * pow(it->position.get_length() - expected_length, 2);
			n += (it->boltzmann_factor / partition_function);
		}
	}
	var /= n;
	// display result
	outfile << superhelicity << ' ' << expected_length << ' ' << var << "\n";
}

void Simulation::sandbox()
{ // test/debug environment
	bool verbose = verbose_flag;
	// sum sequence favorability

	// R-loop length histogram
	ofstream outfile(outfilename, ios::out);
	Gene *this_gene;
	this_gene = new Gene();
	this_gene->read_gene(infile, genebeg, geneend);
	this_gene->windower.set_min_window_size(minlength);

	if (this_gene->getPosition().strand == "-")
	{
		this_gene->reverse_complement_sequence();
	}

	if (complement_flag)
	{
		this_gene->complement_sequence();
	}
	if (reverse_flag)
	{
		this_gene->reverse_sequence();
	}
	genes.push_back(this_gene);
	if (circular_flag)
	{
		this_gene->compute_structures_circular(*models[0], verbose);
	}
	else
	{
		this_gene->compute_structures(*models[0], verbose);
	}
	// determine P(ground state)
	long double partition_function = 0;
	long double ground_state_factor = 0;
	long double sanity_check = 0;
	for (vector<Structure>::iterator it = this_gene->getRloopStructures().begin();
		 it < this_gene->getRloopStructures().end(); ++it)
	{
		partition_function += it->boltzmann_factor;
	}
	ground_state_factor = models[0]->ground_state_factor();

	partition_function += ground_state_factor;
	// sanity check code
	for (vector<Structure>::iterator it = this_gene->getRloopStructures().begin();
		 it < this_gene->getRloopStructures().end(); ++it)
	{
		sanity_check += it->boltzmann_factor / partition_function;
	}
	sanity_check += models[0]->ground_state_factor() / partition_function;
	cout << "\t- "
		 << "Sanity check: " << sanity_check << "\n";
	vector<long double> values;
	values.assign(this_gene->getSequence().size() + 1, 0); // fill vector with 0s
	values[0] = ground_state_factor / partition_function;
	// iterate through structures and record each probability to the appropriate place in the values array
	for (int i = 1; i < this_gene->getRloopStructures().size(); i++)
	{
		values[this_gene->getRloopStructures()[i].position.get_length()] +=
			this_gene->getRloopStructures()[i].boltzmann_factor / partition_function;
	}
	for (int i = 0; i < values.size(); i++)
	{
		outfile << i << ' ' << values[i] << "\n";
	}
	outfile.close();
}
