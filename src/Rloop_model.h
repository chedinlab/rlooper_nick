//
// Created by Robert Stolz on 6/27/17.
//
#ifndef RLOOPER2_RLOOP_EQUILIBRIUM_MODEL_H
#define RLOOPER2_RLOOP_EQUILIBRIUM_MODEL_H

#include <vector>
#include "model.h"
#include <string>
// #include "colors.h"
// #include "mytime.h"


class Rloop_model: public Model{
protected:
	//model parameters
	//equilibrium energetics parameters
	int nick; //nick start location
	int nicklen; // nick length bp (end = nick + nicklen, default is 0)
	int selffoldlen; // self fold length from nick position, default
	int N; //experimentally determined length of the (-) sc domain adter the transcription machinery
	double A; //turns/bp of the B-form double helix
	double C; //tortional stiffness of ssDNA winding. (Could be 3.6 for ds or 1.8 for ss winding)
	double k; //Hooke's law coefficient: (2200*ideal_gas_constant in kcal/mol*absolute_temp_in_kelvin)/N
	double T; //temparature
	double a; //Neucleation Free Energy (junction energy) in Kcals (~3-10.2kCals) 5000
	double sigma; //measurement of energy upstream of replication domain. Moved to up the Model object.
	double alpha; //linking difference: topological parameter. refers to the percent of twists difference.
	//energy difference for the described RNA/DNA dinucleotide pairs (RNA/DNA - DNA/DNA) in kcal/mol as described in Huppert et. al. 2008
	//note, energies are for the RNA (non-template/sense DNA) strand in the 5' to 3' direction
	double rGG_dCC;
	double rGC_dCG;
	double rGA_dCT;
	double rGT_dCA;
	double rCG_dGC;
	double rCC_dGG;
	double rCA_dGT;
	double rCT_dGA;
	double rAG_dTC;
	double rAC_dTG;
	double rAA_dTT;
	double rAT_dTA;
	double rTG_dAC;
	double rTC_dAG;
	double rTA_dAT;
	double rTT_dAA;
	
	double rGG_rCC;
	double rGC_rCG;
	double rGA_rCT;
	double rGT_rCA;
	double rCG_rGC;
	double rCC_rGG;
	double rCA_rGT;
	double rCT_rGA;
	double rAG_rTC;
	double rAC_rTG;
	double rAA_rTT;
	double rAT_rTA;
	double rTG_rAC;
	double rTC_rAG;
	double rTA_rAT;
	double rTT_rAA;
	
	bool homopolymer_override;
	bool unconstrained;
	double override_energy;
	
	//dynamic model specific parameters
	double sigma_total;
	double alpha_total;
	double current_bp_energy;
	double current_junction_energy;
	double current_superhelical_energy;
	double proposed_bp_energy;
	double proposed_junction_energy;
	double proposed_superhelical_energy;
	double partition_function;
	int n_rloops;
	int n_rloop_bases;
	

	int totalSim;
	int current_pos; //The INDEX value of the current "polymerase" position. current_pos + 1 is the base-pair number.
	int initiation_step_size;
	int elongation_step_size;
	// double tx_sigma;
	// double tx_alpha;

	Loci current_loci; //used to store the coordinates of the current rloop structure.
	
public:
	//constructors
	Rloop_model();
	//need a special constructor that lets your specify some or all of these parameters
	
	int getnick() const;
	int getnicklen() const;
	int getN() const;
	double getA() const;
	double getC() const;
	double getK() const;
	double getT() const;
	double geta() const;
	double getSigma() const;
	double getAlpha() const;
	void setnick(int n);
	void setnicklen(int n);
	void setselffoldlen(int n);
	void setN(int N);
	void setA(double A);
	void setC(double C);
	void setK(double k);
	void setT(double T);
	void seta(double a);
	int getnick2();
	int getnicklen2();
	int getselffoldlen2();
	double geta2();
	double getK2();
	double getSigma2();
	
	void set_superhelicity(double sigma);
	void set_unconstrained(bool value);
	void setAlpha(double alpha);
	void set_bp_energy_override(double energy);
	
	//member functions
	//added for eq model
	string join(char &delim, const vector<char> &lst); // perl/python like join
	string join(const vector<char> &lst); // perl/python like join
	vector<char> copy(const vector<char> &sequence);
	void reverse_sequence(vector<char> &sequence);
	void complement_sequence(vector<char> &sequence);
	void reverse_complement_sequence(vector<char> &sequence);
	void step_forward(vector<char> &sequence, vector<char>::iterator& b_0, vector<char>::iterator& b_1);
	
	// original eq model
	int find_distance(vector<char>& sequence,const vector<char>::iterator& first, const vector<char>::iterator& second, Structure& structure);
	double step_forward_bps(const vector<char>::iterator& first, const vector<char>::iterator& second);
	double compute_bps_interval(const char &first, const char &second);
	void compute_structure(vector<char>& sequence, const std::vector<char>::iterator &start, const std::vector<char>::iterator &stop, Structure& structure, float &bp_energy, float &time1, float &time2, float &time3, bool &verbose);
	void compute_external_structure(Structure& structure, Structure& rloop, Peak& external, bool &verbose);
	void compute_residuals(Structure& structure, bool &verbose);
	void ground_state_residuals(double& twist, double& writhe);
	long double ground_state_factor();
	long double ground_state_energy();
	
	// dynamic model
	//constructors
	//Rloop_dynamic_model();
	
	bool in_rloop;
	vector<char> seqD;
	vector<Loci> rloop_structures;
	double ambient_alpha;
	double ambient_sigma;
	int currSim;
	int iter;
	double tx_sigma;
	double tx_alpha;
	double tx_ambient_sigma;
	double tx_ambient_alpha;

	stringstream dynamic_write_buffer;
	int dynamic_window_size;
	
	int dynamic_getCurrentPos() const;
	
	void dynamic_setIter(int value);
	int dynamic_getIter() const;

	void dynamic_setTotalSim(string value);
	int dynamic_getTotalSim() const;

	void dynamic_setCurrSim(int value);
	int dynamic_getCurrSim() const;
	
	void dynamic_setAlpha(double value);
	double dynamic_getAlpha() const;

	void dynamic_setSigmaTotal(double value);
	double dynamic_getSigmaTotal() const;

	void dynamic_setSigma(double value);
	double dynamic_getSigma() const;

	void dynamic_setTx_Sigma(double value);
	double dynamic_getTx_Sigma() const;

	void dynamic_setTx_Alpha(double value);
	double dynamic_getTx_Alpha() const;

	void dynamic_setAlphaTotal(double value);
	double dynamic_getAlphaTotal() const;

	double dynamic_getTx_Amb_Sigma() const;
	void dynamic_setTx_Amb_Sigma(double value);

	double dynamic_getTx_Amb_Alpha() const;
	void dynamic_setTx_Amb_Alpha(double value);

	void dynamic_reset_model();
	void dynamic_step_forward_initiation(vector<vector<double>> &df_bfs, vector<vector<double>> &df_bftotals, vector<vector<double>> &df_probs, bool &verbose, int &seed, int &begpos, double &myurn, double &myprob, std::vector<std::string> &values);//, std::stringstream &ss_dynamic_log);
	bool dynamic_step_forward_elongation(vector<vector<double>> &df_bfs, vector<vector<double>> &df_bftotals, vector<vector<double>> &df_probs, bool &verbose, int &seed, int &endpos, double &myurn, double &myprob, std::vector<std::string> &values);//, std::stringstream &ss_dynamic_log);
	void dynamic_print_topological_state(std::stringstream &ss_dynamic_log);
	// void Rloop_model::dynamic_update_topology();

	double dynamic_compute_res_lk();
	
	//automatically generated getters and setters vvv (refactor later)
	int dynamic_getWindow_size() const;
	void dynamic_setWindow_size(int Dynamic_window_size);
	int dynamic_getN_rloops() const;
	void dynamic_setN_rloops(int n_rloops);
	int dynamic_getInitiation_step_size() const;
	void dynamic_setInitiation_step_size(int initiation_step_size);
	int dynamic_getElongation_step_size() const;
	void dynamic_setElongation_step_size(int elongation_step_size);
	int dynamic_getN_rloop_bases() const;
	void dynamic_setN_rloop_bases(int n_rloop_bases);


};
#endif //RLOOPER2_MODEL_H
