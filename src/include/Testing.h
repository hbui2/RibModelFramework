#ifndef Testing_H
#define Testing_H


#include "Utility.h"
#include "SequenceSummary.h"
#include "Gene.h"
#include "Genome.h"
<<<<<<< HEAD


void testSequenceSummary();
void testGene();
void testGenome(std::string testFileDir);
void testCodonTable();
void initMutation(std::vector<double> mutationValues, unsigned mixtureElement, std::string aa);


=======
#include "base/Parameter.h"
#include "CovarianceMatrix.h"
#include "MCMCAlgorithm.h"


int testUtility();
int testSequenceSummary();
int testGene();
void testGenomePAHelper(Genome* genome, bool simulated); //Only used internally for testGenome
int testGenome(std::string testFileDir);
int testParameter();
//int testParameterWithFile(std::string filename); //TODO: Rework or remove
int testCovarianceMatrix();
//int testPAParameter(); //TODO: Rework or remove
int testMCMCAlgorithm();
>>>>>>> upstream/master

//Blank header
#endif // Testing_H
