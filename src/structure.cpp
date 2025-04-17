//
// Created by Robert Stolz on 6/27/17.
// a

#include "structure.h"

Loci::Loci(std::string C, std::string c, long int s, long int e): chromosome(C), strand(c), start_pos(s), end_pos(e) {}

int Loci::get_length(){
    return end_pos - start_pos;
}

Peak::Peak(Loci l, int i): position(l), intensity(i) {}

Structure::Structure(Loci l, float g, float e, float f, float b, float p) {
    position = l;
    gsigma = g;
    bp_energy = e;
    free_energy = f;
    boltzmann_factor = b;
    probability = p;
}
