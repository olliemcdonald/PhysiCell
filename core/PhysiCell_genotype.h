/*
Additional file to deal with genotype changes
*/

#ifndef __PhysiCell_genotype_h__
#define __PhysiCell_genotype_h__

#include <vector>
#include <string>

#include "../BioFVM/BioFVM.h"
#include "./PhysiCell_utilities.h"

using namespace BioFVM;

namespace PhysiCell{

class Cell;

class Genotype
{
private:
public:
  std::string genotype;
  double cna_prob;
  double snv_prob;

  int num_events;  // cell-specific counter of number of alterations
  // cna_counter and snv_counter keeps track up the current alteration number to
  //append to new clones shared between instances of Genotype
  static int cna_counter;
  static int snv_counter;

  double punctuated_probability;
  double punctuated_poisson_parameter;

  Genotype();

  void (*alter_genotype)( Cell* pCell, Genotype& genotype);

  void copynumber_change();
  void update_copynumber_gradual();
  void update_copynumber_punctuation();
};

};

#endif
