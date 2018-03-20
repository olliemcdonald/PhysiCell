/*
Additional file to deal with genotype changes
*/
#include "./PhysiCell_genotype.h"

#include "../BioFVM/BioFVM.h"
#include "./PhysiCell_constants.h"

using namespace BioFVM;

namespace PhysiCell{
  int Genotype::cna_counter;
  int Genotype::snv_counter;

Genotype::Genotype()
{
  genotype = "";
  cna_prob = 0.0;
  snv_prob = 0.0;
  int cna_counter = 0;
  int snv_counter = 0;

  punctuated_probability = 0.0;
  punctuated_poisson_parameter = 0;

  alter_genotype = NULL;
}

void Genotype::copynumber_change()
{
  num_events++;
  cna_counter++;
  double runif_punctuation = UniformRandom();
  if (runif_punctuation < punctuated_probability)
  {
    update_copynumber_punctuation();
  }
  else
  {
    update_copynumber_gradual();
  }

}

void Genotype::update_copynumber_gradual()
{
  genotype = genotype + ">" + std::to_string(cna_counter);

  // include fitness adjustments later
}

void Genotype::update_copynumber_punctuation()
{
  int pois_rv = PoissonRandom(punctuated_poisson_parameter) + 1;  // rand_pois;
  genotype = genotype + ">" + std::to_string(cna_counter) + "." +
             std::to_string(pois_rv);

  // include fitness adjustments later
}

};
