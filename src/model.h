//
// Created by Robert Stolz on 6/28/17.
//

#ifndef RLOOPER2_MODEL_H
#define RLOOPER2_MODEL_H

#include "structure.h"
#include "windower.h"
#include "biophysics.h"

class Model{ //abstract class
protected:
    std::vector<char>* target_sequence; //keeps track of what gene we are in
public:
    virtual void compute_structure(std::vector<char>& sequence, const std::vector<char>::iterator &start, const std::vector<char>::iterator &stop, Structure& structure, float &bp_energy, float &time1, float &time2, float &time3, bool &verbose) = 0;
    virtual void compute_external_structure(Structure& structure, Structure& rloop, Peak& external, bool &verbose) = 0;
    virtual void compute_residuals(Structure& structure, bool &verbose) = 0;
    virtual void ground_state_residuals(double& twist, double& writhe) = 0;
    virtual long double ground_state_factor() = 0;
    virtual long double ground_state_energy() = 0;
    virtual double geta2() = 0;
    virtual int getnick2() = 0;
    virtual int getnicklen2() = 0;
    virtual int getselffoldlen2() = 0;
    virtual double getK2() = 0;
    virtual double getSigma2() = 0;
    virtual void set_superhelicity(double sigma) = 0;
    virtual void set_unconstrained(bool value) = 0;

};

#endif //RLOOPER2_MODEL_H
