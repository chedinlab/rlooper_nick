//
// Created by Robert Stolz on 6/27/17.
// Modified & maintained by Stella Hartono 2023-2024
//

#ifndef RLOOPER2_GENE_H
#define RLOOPER2_GENE_H
#include "structure.h"
#include "exception_handling.h"
#include "model.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

void declarestring();
// void declarestring(string &coldotYW, string &coldotLPR,string &coldotLBU,string &coldotLRD,string &coldotLCY,string &coldotLGN, string &coldotNN) {

class Gene{
private:
    //private member functions
    std::string gene_name;
    std::string header;
    Loci position;
    std::vector<char> sequence;
    std::vector<Structure> rloop_structures;
    long double ground_state_energy;
    
private:

    /**
     * parses FASTA headers into the Gene's Loci member.
     */
    void parse_header();

public:
    int fabeg;
    int faend;

    Windower windower;
    //constructors
    Gene();
    //destructors
    ~Gene();
    //Gene(std::fstream* fastafile);
    //getters and setters
    string getName();
    const string &getHeader() const;
    void setHeader(const string &header);
    const Loci &getPosition() const;
    void setPosition(const Loci &position);
    const vector<char, allocator<char>> &getSequence() const;
    //void setSequence(const vector<char, allocator<char>> &sequence);
    vector<Structure>& getRloopStructures();
    //member functions
    /**
     * Reads the next FASTA record from the input file. Calls parse_header.
     * @param   fastafile   An address to the open ifstream with at least one FASTA record in it.
     * @return              returns a boolean indicating if the gene was teh last FASTA record in the file.
     */
    bool read_gene(std::ifstream& fastafile, int genebeg, int geneend);

    /**
     * useful for debugging. Prints the header and the sequence stored in the current gene.
     */
    void print_gene();

    /**
     * Takes a model and uses it to allocate and populate a vector<Structure> of structures. Pushes the result to structures.
     * @param   model   a Model object with an implemented compute_structure method.
     */
    void compute_structures(Model& model, bool &verbose);

    vector<Structure> compute_structures_dynamic(Model& model, vector<char> input_sequence, bool &verbose);

    void compute_structures_circular(Model& model, bool &verbose);

    void compute_structures_external(vector<Peak> &external_structures, Model &model, bool &verbose);

    /**
     * computes residual twist and superhelicity for the ensemble
     */
    void compute_residuals(Model& model, bool &verbose);

    /**
     * Clears all structures associated with the gene, correctly deallocating them in memory.
     */
    void clear_structures();

    /**
     * Computes the GC skew of the provided sequence.
     * @return  returns a float representing the GC skew.
     */
    float compute_GC_skew();

    /**
     * Computes the UY skew of the provided sequence.
     * @return  returns a float representing the GC skew.
     */
    float compute_UY_skew();

    float compute_GC_content();

    /**
     * complements the sequence data (A<->T, G<->C)
     */
    void complement_sequence();

    /**
     * reverses the sequence data (GATTACA <-> ACATTAG)
     */
    void reverse_sequence();

    /**
     * reverse complement the sequence data (GATTACA <-> TGTAATC)
     * intentionally put here per biologist tradition
     * basically just complement_sequence() then reverse_sequence()
     */
    void reverse_complement_sequence();

    /**
     * Returns the length of the gene
     * @return  returns an int computed from the gene's location
     */
    int get_length();

    void clear_sequence();

    void dump_structures(string outfilename);

};


#endif //RLOOPER2_GENE_H
