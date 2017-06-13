#include "include/CodonTable.h"
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

#include <iostream>

std::vector <std::string> CodonTable::defaultVector;

//----------------------------------------------//
//----------Constructors & Destructors----------//
//----------------------------------------------//
CodonTable::CodonTable()
{
    tableId = 0; //standard codon table by NCBI
    splitAA = true;
}


CodonTable::CodonTable(unsigned _tableId, std::string model, bool _splitAA, std::vector <std::string> _groupList) : tableId(_tableId), splitAA(_splitAA)
{
	//We assume the user gives the id of the codon table as described by the NCBI tables, which are 1 indexed, not zero.
	//Some of the tables are skipped so we check for invalid input here.
	if(tableId == 7 || tableId == 8 || tableId == 15 || tableId == 17 || tableId == 18 || tableId == 19 || tableId == 20)
	{
        std::cerr << "Invalid codon table: " << tableId << " using default codon table (NCBI codon table 1)\n";
        tableId = 1; //standard codon table by NCBI
    }
    tableId--; //Make the table id zero indexed from here.

	//TODO: error check the user giving something completely invalid, such as zero.


	//If the user does not pass a groupList...
	if (_groupList == defaultVector)
	{
		if (model == "RFP")
		{
			groupList = { "GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
				"TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
				"ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
				"AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
				"CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
				"ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
				"AGT", "TAA", "TAG", "TGA" };
		}
		else // currently ROC and FONSE 
		{
			groupList = splitAA ? defaultSplitAAGroupListings[tableId] : defaultAAGroupListings[tableId];
		}
	} //End of default group lists.
	else //User has provided a groupList to use (presumably a subset of the entire group list as shown above).
	{
		groupList = _groupList; //TODO: error check given lists
	}
}


CodonTable::~CodonTable()
{
    //dtor
}


CodonTable::CodonTable(const CodonTable& other)
{
    tableId = other.tableId;
    splitAA = other.splitAA;
}


CodonTable& CodonTable::operator=(const CodonTable& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    tableId = rhs.tableId;
    splitAA = rhs.splitAA;
    return *this;
}



//------------------------------------//
//----------Getter Functions----------//
//------------------------------------//
unsigned  CodonTable::getTableId()
{
    return tableId;
}


bool CodonTable::getSplitAA()
{
    return splitAA;
}

std::map<std::string, std::vector<std::string>> CodonTable::getAAToCodonMap()
{
	return AAToCodonMap;
}

std::map<std::string, std::vector<unsigned>> CodonTable::getAAToCodonIndexMap()
{
	return AAToCodonIndexMap;
}

std::map<std::string, std::vector<std::string>> CodonTable::getAAToCodonMapWithoutReference()
{
	return AAToCodonMapWithoutReference;
}

std::map<std::string, std::vector<unsigned>> CodonTable::getAAToCodonIndexMapWithoutReference()
{
	return AAToCodonIndexMapWithoutReference;
}

std::map <std::string, std::string> CodonTable::getCodonToAAMap()
{
    return codonToAAMap;
}

std::vector<std::string> CodonTable::getGroupList()
{
	return groupList;
}


std::map <std::string, unsigned> CodonTable::getAAMap()
{
    return AAMap;
}


std::map <std::string, unsigned> CodonTable::getAAToNumCodonsMap()
{
    return AAToNumCodonsMap;
}

unsigned CodonTable::getNumCodonsForAA(std::string aa, bool withoutReference)
{
	unsigned rv;
    unsigned numCodons = AAToNumCodonsMap.find(aa) -> second;

	if (numCodons == 0)
	{
		rv = 0;
	}
	else
	{
		rv = withoutReference ? numCodons - 1 : numCodons;
	}
	return rv;
}


unsigned CodonTable::getNumCodonsForAAIndex(unsigned aaIndex, bool withoutReference)
{
    std::string aa = indexToAA(aaIndex);
    return getNumCodonsForAA(aa, withoutReference);
}


//--------------------------------------//
//----------Mapping Operations----------//
//--------------------------------------//
unsigned CodonTable::AAToAAIndex(std::string aa)
{
	unsigned rv;
	if (AAMap.find(aa) == AAMap.end())
	{
		std::cerr << "Warning: invalid AA in CodonTable::AAToAAIndex";
		rv = 42; // 42 is the error check since AAs will go from 0 to around 30, but never above 40. We just like 42.
	}
	else
	{
		rv = AAMap.find(aa)->second;
	}
	return rv;
}


std::vector <unsigned> CodonTable::AAIndexToCodonRange(unsigned aaIndex, bool withoutReference)
{
    return withoutReference ? AAToCodonIndexMapWithoutReference.find(aminoAcidArray[aaIndex])->second : AAToCodonIndexMap.find(aminoAcidArray[aaIndex])->second;
}


std::string CodonTable::indexToCodon(unsigned index)
{
    return CodonTable::codonArray[index];
}


std::vector <unsigned> CodonTable::AAToCodonRange(std::string aa, bool withoutReference)
{
    unsigned aaIndex = AAToAAIndex(aa);
    return AAIndexToCodonRange(aaIndex, withoutReference);
}


std::vector<std::string> CodonTable::AAToCodon(std::string aa, bool withoutReference)
{
    std::vector <std::string> codons;
    std::vector <unsigned> aaRange = AAToCodonRange(aa, withoutReference);
	for (unsigned i = 0; i < aaRange.size(); i++)
	{
		codons.push_back(indexToCodon(aaRange[i]));
	}
    return codons;
}


std::string CodonTable::codonToAA(std::string& codon)
{
    return codonToAAMap.find(codon) -> second;
}


unsigned CodonTable::codonToIndex(std::string& codon)
{
	unsigned rv = 64;
	std::map <std::string, unsigned int>::const_iterator mit;
	mit = codonToIndexWithReference.find(codon);
	if (mit != codonToIndexWithReference.end()) rv = mit->second;
    return rv;
}


unsigned CodonTable::codonToAAIndex(std::string& codon)
{
    std::string AA = codonToAA(codon);
    return AAMap.find(AA) -> second;
}


std::string CodonTable::indexToAA(unsigned aaIndex)
{
    return aminoAcidArray[aaIndex];
}



//-----------------------------------//
//----------Other Functions----------//
//-----------------------------------//
void CodonTable::setupCodonTable()
{
	std::vector <unsigned> codonIndices;
	std::vector <std::string> codons;

	std::vector <std::vector <std::vector <std::string> > > ctl = splitAA ? codonTableListing : codonTableListingWithoutSplit;

	for (unsigned groupIndex = 0; groupIndex < aminoAcidArray.size(); groupIndex++) {
		unsigned AANum = fullAAMap.find(aminoAcidArray[groupIndex])->second; //TODO: error check
		AAToNumCodonsMap.insert(std::make_pair(aminoAcidArray[groupIndex], numCodonsPerAAForTable[tableId][AANum]));
	}

	// setup AAToCodon

	for (unsigned groupIndex = 0; groupIndex < aminoAcidArray.size(); groupIndex++)
	{
		unsigned AAIndex = fullAAMap.find(aminoAcidArray[groupIndex])->second;
		AAToCodonMap.insert(std::make_pair(aminoAcidArray[groupIndex], ctl[tableId][AAIndex]));
		for (unsigned codonId = 0; codonId < ctl[tableId][AAIndex].size(); codonId++)
		{
			codonIndices.push_back(codonToIndexWithReference.find(ctl[tableId][AAIndex][codonId])->second);
		}

		AAToCodonIndexMap.insert(std::make_pair(aminoAcidArray[groupIndex], codonIndices));
		codonIndices.resize(0);
	}

	for (unsigned groupIndex = 0; groupIndex < aminoAcidArray.size(); groupIndex++)
	{
		unsigned AAIndex = fullAAMap.find(aminoAcidArray[groupIndex])->second;
		if (ctl[tableId][AAIndex].size() != 0)
		{
			for (unsigned codonId = 0; codonId < ctl[tableId][AAIndex].size() - 1; codonId++)
			{
				codons.push_back(ctl[tableId][AAIndex][codonId]);
				codonIndices.push_back(codonToIndexWithReference.find(ctl[tableId][AAIndex][codonId])->second);
			}
		}
		AAToCodonMapWithoutReference.insert(std::make_pair(aminoAcidArray[groupIndex], codons));
		AAToCodonIndexMapWithoutReference.insert(std::make_pair(aminoAcidArray[groupIndex], codonIndices));
		codonIndices.resize(0);
		codons.resize(0);
	}

	for (unsigned groupIndex = 0; groupIndex < aminoAcidArray.size(); groupIndex++) {
		std::string AA = aminoAcidArray[groupIndex];
		std::vector <std::string> codonList = AAToCodonMap.find(AA)->second;
		for (unsigned codonIndex = 0; codonIndex < codonList.size(); codonIndex++) {
			std::string codon = codonList[codonIndex];
			codonToAAMap.insert(std::make_pair(codon, AA));
		}
	}
}


bool CodonTable::checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound)
{
    bool check = false;
    if (lowerbound <= index && index <= upperbound)
    {
        check = true;
    }
    else
    {
        std::cerr <<"Error with the index\nGIVEN: " << index <<"\n";
        std::cerr <<"MUST BE BETWEEN:	" << lowerbound << " & " << upperbound <<"\n";
    }
    return check;
}



//----------------------------------------------//
//---------Static Variables & Functions---------//
//----------------------------------------------//
CodonTable* CodonTable::codonTable;
const std::string CodonTable::Ser2 = "Z";
const std::string CodonTable::Ser1 = "J";
const std::string CodonTable::Thr4_1 = "T";
const std::string CodonTable::Thr4_2 = "B";
const std::string CodonTable::Leu1 = "U";

const std::vector <std::string> CodonTable::aminoAcidArray = {
        "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", CodonTable::Leu1, "M", "N", "P", "Q", "R",
        CodonTable::Ser1, CodonTable::Ser2, "S", "T", CodonTable::Thr4_1, CodonTable::Thr4_2, "V", "W", "Y", "X"};
const std::vector <std::string> CodonTable::aminoAcidArrayWithoutSplit = {
        "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X"};

const std::map<std::string, unsigned> CodonTable::fullAAMap = { { "A", 0 }, { "C", 1 }, { "D", 2 }, { "E", 3 }, { "F", 4 },
{ "G", 5 }, { "H", 6 }, { "I", 7 }, { "K", 8 }, { "L", 9 }, { CodonTable::Leu1, 10}, { "M", 11 }, { "N", 12 }, { "P", 13 },
{ "Q", 14 }, { "R", 15 }, { CodonTable::Ser1, 16 }, { CodonTable::Ser2, 17 }, { "S", 18 }, { "T", 19 }, { CodonTable::Thr4_1, 20 },
{ CodonTable::Thr4_2, 21 }, { "V", 22 }, { "W", 23 }, { "Y", 24 }, { "X", 25 } };


const std::map<std::string, unsigned> CodonTable::codonToIndexWithReference = {{"GCA", 0}, {"GCC", 1}, {"GCG", 2},
   {"GCT", 3}, {"TGC", 4}, {"TGT", 5}, {"GAC", 6}, {"GAT", 7}, {"GAA", 8}, {"GAG", 9}, {"TTC", 10}, {"TTT", 11},
   {"GGA", 12}, {"GGC", 13}, {"GGG", 14}, {"GGT", 15}, {"CAC", 16}, {"CAT", 17}, {"ATA", 18}, {"ATC", 19}, {"ATT", 20},
   {"AAA", 21}, {"AAG", 22}, {"CTA", 23}, {"CTC", 24}, {"CTG", 25}, {"CTT", 26}, {"TTA", 27}, {"TTG", 28}, {"ATG", 29},
   {"AAC", 30}, {"AAT", 31}, {"CCA", 32}, {"CCC", 33}, {"CCG", 34}, {"CCT", 35}, {"CAA", 36}, {"CAG", 37}, {"AGA", 38},
   {"AGG", 39}, {"CGA", 40}, {"CGC", 41}, {"CGG", 42}, {"CGT", 43}, {"TCA", 44}, {"TCC", 45}, {"TCG", 46}, {"TCT", 47},
   {"ACA", 48}, {"ACC", 49}, {"ACG", 50}, {"ACT", 51}, {"GTA", 52}, {"GTC", 53}, {"GTG", 54}, {"GTT", 55}, {"TGG", 56},
   {"TAC", 57}, {"TAT", 58}, {"AGC", 59}, {"AGT", 60}, {"TAA", 61}, {"TAG", 62}, {"TGA", 63}};


const std::string CodonTable::codonArray[] =
        {"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
         "TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
         "ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
         "AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
         "CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
         "ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
         "AGT", "TAA", "TAG", "TGA"};


// http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
const std::vector <std::string> CodonTable::codonTableDefinition = {"1. The Standard Code",
   "2. The Vertebrate Mitochondrial Code", "3. The Yeast Mitochondrial Code",
   "4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
   "5. The Invertebrate Mitochondrial Code", "6. The Ciliate, Dasycladacean and Hexamita Nuclear Code",
   "7. Invalid Codon Table", "8. Invalid Codon Table", "9. The Echinoderm and Flatworm Mitochondrial Code",
   "10. The Euplotid Nuclear Code", "11. The Bacterial, Archaeal and Plant Plastid Code",
   "12. The Alternative Yeast Nuclear Code", "13. The Ascidian Mitochondrial Code",
   "14. The Alternative Flatworm Mitochondrial Code", "15. Invalid Codon Table", "16. Chlorophycean Mitochondrial Code",
   "17. Invalid Codon Table", "18. Invalid Codon Table", "19. Invalid Codon Table", "20. Invalid Codon Table",
   "21. Trematode Mitochondrial Code", "22. Scenedesmus obliquus Mitochondrial Code",
   "23. Thraustochytrium Mitochondrial Code", "24. Pterobranchia Mitochondrial Code",
   "25. Candidate Division SR1 and Gracilibacteria Code"};

const std::map <std::string, unsigned> CodonTable::AAMap = { { "A", 0 }, { "C", 1 }, { "D", 2 }, { "E", 3 }, { "F", 4 }, { "G", 5 }, { "H", 6 }, { "I", 7 },
{ "K", 8 }, { "L", 9 }, { CodonTable::Leu1, 10 }, { "M", 11 }, { "N", 12 }, { "P", 13 }, { "Q", 14 }, { "R", 15 },
{ CodonTable::Ser1, 16 }, { CodonTable::Ser2, 17 }, { "S", 18 }, { "T", 19 }, { CodonTable::Thr4_1, 20 }, { CodonTable::Thr4_2, 21 },
{ "V", 22 }, { "W", 23 }, { "Y", 24 }, { "X", 25 } };

const unsigned CodonTable::numCodonsPerAAForTable[25][26] = {
        {4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,3}, // 1. The Standard Code
        {4,2,2,2,2,4,2,2,2,6,0,2,2,4,2,4,0,2,4,4,0,0,4,2,2,4}, // 2. The Vertebrate Mitochondrial Code
        {4,2,2,2,2,4,2,2,2,2,0,2,2,4,2,6,0,2,4,0,4,4,4,2,2,2}, // 3. The Yeast Mitochondrial Code
        {4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,2,2,2}, // 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
        {4,2,2,2,2,4,2,2,2,6,0,2,2,4,2,4,0,4,4,4,0,0,4,2,2,2}, // 5. The Invertebrate Mitochondrial Code
        {4,2,2,2,2,4,2,3,2,6,0,1,2,4,4,6,0,2,4,4,0,0,4,1,2,1}, // 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 7. Invalid Codon Table
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 8. Invalid Codon Table
        {4,2,2,2,2,4,2,3,1,6,0,1,3,4,2,4,0,4,4,4,0,0,4,2,2,2}, // 9. The Echinoderm and Flatworm Mitochondrial Code
        {4,3,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2}, // 10. The Euplotid Nuclear Code
        {4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,3}, // 11. The Bacterial, Archaeal and Plant Plastid Code
        {4,2,2,2,2,4,2,3,2,5,0,1,2,4,2,6,1,2,4,4,0,0,4,1,2,3}, // 12. The Alternative Yeast Nuclear Code
        {4,2,2,2,2,6,2,2,2,6,0,2,2,4,2,4,0,2,4,4,0,0,4,2,2,2}, // 13. The Ascidian Mitochondrial Code
        {4,2,2,2,2,4,2,3,1,6,0,1,3,4,2,4,0,4,4,4,0,0,4,2,3,1}, // 14. The Alternative Flatworm Mitochondrial Code
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 15. Invalid Codon Table
        {4,2,2,2,2,4,2,3,2,6,1,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2}, // 16. Chlorophycean Mitochondrial Code
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 17. Invalid Codon Table
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 18. Invalid Codon Table
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 19. Invalid Codon Table
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 20. Invalid Codon Table
        {4,2,2,2,2,4,2,2,1,6,0,2,3,4,2,4,0,4,4,4,0,0,4,2,2,2}, // 21. Trematode Mitochondrial Code
        {4,2,2,2,2,4,3,2,2,6,1,1,2,4,2,6,0,2,3,4,0,0,4,1,2,3}, // 22. Scenedesmus obliquus Mitochondrial Code
        {4,2,2,2,2,4,2,3,2,5,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,4}, // 23. Thraustochytrium Mitochondrial Code
        {4,2,2,2,2,4,2,3,3,6,0,1,2,4,2,4,0,3,4,4,0,0,4,2,2,2}, // 24. Pterobranchia Mitochondrial Code
        {4,2,2,2,2,5,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2}, // 25. Candidate Division SR1 and Gracilibacteria Code
};

const std::vector <std::vector <std::string> > CodonTable::defaultAAGroupListings = {
	// 1. The Standard Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	"S", "T", "V", "Y" },
	// 2. The Vertebrate Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
	"S", "T", "V", "W", "Y" },
	// 3. The Yeast Mitochondrial Cod
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
	"S", "T", "V", "W", "Y" },
	// 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	"S", "T", "V", "W", "Y" },
	// 5. The Invertebrate Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
	"S", "T", "V", "W", "Y" },
	// 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	"S", "T", "V", "Y" },
	// 7. Invalid Codon Table
	{},
	// 8. Invalid Codon Table
	{},
	// 9. The Echinoderm and Flatworm Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "L", "N", "P", "Q", "R",
	"S", "T", "V", "W", "Y" },
	// 10. The Euplotid Nuclear Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	"S", "T", "V", "Y" },
	// 11. The Bacterial, Archaeal and Plant Plastid Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	"S", "T", "V", "Y" },
	// 12. The Alternative Yeast Nuclear Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	"S", "T", "V", "Y" },
	// 13. The Ascidian Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
	"S", "T", "V", "W", "Y" },
	// 14. The Alternative Flatworm Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "L", "N", "P", "Q", "R",
	"S", "T", "V", "W", "Y" },
	// 15. Invalid Codon Table
	{},
	// 16. Chlorophycean Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	"S", "T", "V", "Y" },
	// 17. Invalid Codon Table
	{},
	// 18. Invalid Codon Table
	{},
	// 19. Invalid Codon Table
	{},
	// 20. Invalid Codon Table
	{},
	// 21. Trematode Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "L", "M", "N", "P", "Q", "R",
	"S", "T", "V", "W", "Y" },
	// 22. Scenedesmus obliquus Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	"S", "T", "V", "Y" },
	// 23. Thraustochytrium Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	"S", "T", "V", "Y" },
	// 24. Pterobranchia Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	"S", "T", "V", "W", "Y" },
	// 25. Candidate Division SR1 and Gracilibacteria Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	"S", "T", "V", "Y" },
};

const std::vector <std::vector <std::string> > CodonTable::defaultSplitAAGroupListings = {
	// 1. The Standard Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "Y" },
	// 2. The Vertebrate Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "W", "Y" },
	// 3. The Yeast Mitochondrial Cod
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", CodonTable::Thr4_1, CodonTable::Thr4_2, "V", "W", "Y" },
	// 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "W", "Y" },
	// 5. The Invertebrate Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "W", "Y" },
	// 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "Y" },
	// 7. Invalid Codon Table
	{},
	// 8. Invalid Codon Table
	{},
	// 9. The Echinoderm and Flatworm Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "W", "Y" },
	// 10. The Euplotid Nuclear Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "Y" },
	// 11. The Bacterial, Archaeal and Plant Plastid Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "Y" },
	// 12. The Alternative Yeast Nuclear Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "Y" },
	// 13. The Ascidian Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "W", "Y" },
	// 14. The Alternative Flatworm Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "W", "Y" },
	// 15. Invalid Codon Table
	{},
	// 16. Chlorophycean Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "Y" },
	// 17. Invalid Codon Table
	{},
	// 18. Invalid Codon Table
	{},
	// 19. Invalid Codon Table
	{},
	// 20. Invalid Codon Table
	{},
	// 21. Trematode Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "L", "M", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "W", "Y" },
	// 22. Scenedesmus obliquus Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "Y" },
	// 23. Thraustochytrium Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "Y" },
	// 24. Pterobranchia Mitochondrial Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "W", "Y" },
	// 25. Candidate Division SR1 and Gracilibacteria Code
	{ "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R",
	CodonTable::Ser2, "S", "T", "V", "Y" },
};

const std::vector <std::vector <std::vector <std::string> > > CodonTable::codonTableListing = {
	// 1. The Standard Code
	{
		{"GCT", "GCC", "GCA", "GCG"}, // A
		{"TGT", "TGC" }, // C
		{"GAT", "GAC" }, // D
		{"GAA", "GAG"}, // E
		{"TTT", "TTC"}, // F
		{"GGT", "GGC", "GGA", "GGG" }, // G
		{"CAT", "CAC" }, // H
		{"ATT", "ATC", "ATA"}, // I
		{"AAA", "AAG" }, // K
		{"TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{"ATG"}, // M
		{"AAT", "AAC"}, // N
		{"CCT", "CCC", "CCA", "CCG"}, // P
		{"CAA", "CAG"}, // Q
		{"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}, // R
		{}, // CodonTable::Ser1
		{"AGT", "AGC"}, // CodonTable::Ser2
		{"TCT", "TCC", "TCA", "TCG"}, // S
		{"ACT", "ACC", "ACA", "ACG"}, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{"GTT", "GTC", "GTA", "GTG"}, // V
		{"TGG"}, // W
		{"TAT", "TAC"}, // Y
		{"TAA", "TAG", "TGA"}, // X
	},
	
	// 2. The Vertebrate Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG", "ATA" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG", "AGA", "AGG"}, // X
	},

	// 3. The Yeast Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG", "ATA" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{ "CTT", "CTC", "CTA", "CTG" }, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
	// 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
			// 5. The Invertebrate Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG", "ATA" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC", "AGA", "AGG" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},

			// 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG", "TAA", "TAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TGA" }, // X
	},
	// 7. Invalid Codon Table
	{},
	// 8. Invalid Codon Table
	{},
	// 9. The Echinoderm and Flatworm Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC", "AAA" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG"}, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC", "AGA", "AGG" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
	// 10. The Euplotid Nuclear Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC", "TGA" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG"}, // X
	},
	// 11. The Bacterial, Archaeal and Plant Plastid Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG", "TGA" }, // X
	},
	// 12. The Alternative Yeast Nuclear Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{ "CTG" }, // CodonTable::Ser1
		{ "AGT", "AGC" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG", "TGA" }, // X
	},
	// 13. The Ascidian Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG", "AGA", "AGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG", "ATA" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
	// 14. The Alternative Flatworm Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC", "AAA" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG"}, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC", "AGA", "AGG" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC", "TAA" }, // Y
		{ "TAG" }, // X
	},
	// 15. Invalid Codon Table
	// 16. Chlorophycean Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{ "TAG" }, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TGA" }, // X
	},
	// 17. Invalid Codon Table
	// 18. Invalid Codon Table
	// 19. Invalid Codon Table
	// 20. Invalid Codon Table
	// 21. Trematode Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC" }, // I
		{ "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG", "ATA" }, // M
		{ "AAT", "AAC", "AAA" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG"}, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC", "AGA", "AGG" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG"}, // X
	},
	// 22. Scenedesmus obliquus Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{ "TAG" }, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TGA", "TCA" }, // X
	},
	// 23. Thraustochytrium Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG", "TGA", "TTA" }, // X
	},
	// 24. Pterobranchia Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG", "AGG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG"}, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC", "AGA" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
	// 25. Candidate Division SR1 and Gracilibacteria Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG", "TGA" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
};

const std::vector <std::vector <std::vector <std::string> > > CodonTable::codonTableListingWithoutSplit = {
	// 1. The Standard Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG", "TGA" }, // X
	},

	// 2. The Vertebrate Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG", "ATA" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG", "AGA", "AGG" }, // X
	},

	// 3. The Yeast Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG", "ATA" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" }, // S
		{ "ACT", "ACC", "ACA", "ACG", "CTT", "CTC", "CTA", "CTG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
	// 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
	// 5. The Invertebrate Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG", "ATA" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG" }, // R
		{}, // CodonTable::Ser1
		{ "AGT", "AGC", "AGA", "AGG" }, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "AGA", "AGG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},

	// 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG", "TAA", "TAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TGA" }, // X
	},
	// 7. Invalid Codon Table
	{},
	// 8. Invalid Codon Table
	{},
	// 9. The Echinoderm and Flatworm Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC", "AAA" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "AGA", "AGG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
	// 10. The Euplotid Nuclear Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC", "TGA" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
	// 11. The Bacterial, Archaeal and Plant Plastid Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG", "TGA" }, // X
	},
	// 12. The Alternative Yeast Nuclear Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "CTG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG", "TGA" }, // X
	},
	// 13. The Ascidian Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG", "AGA", "AGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG", "ATA" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
	// 14. The Alternative Flatworm Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC", "AAA" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "AGA", "AGG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC", "TAA" }, // Y
		{ "TAG" }, // X
	},
	// 15. Invalid Codon Table
	// 16. Chlorophycean Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "TAG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TGA" }, // X
	},
	// 17. Invalid Codon Table
	// 18. Invalid Codon Table
	// 19. Invalid Codon Table
	// 20. Invalid Codon Table
	// 21. Trematode Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC" }, // I
		{ "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG", "ATA" }, // M
		{ "AAT", "AAC", "AAA" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "AGA", "AGG" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
	// 22. Scenedesmus obliquus Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "TAG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCG", "AGT", "AGC" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TGA", "TCA" }, // X
	},
	// 23. Thraustochytrium Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG", "TGA", "TTA" }, // X
	},
	// 24. Pterobranchia Mitochondrial Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG", "AGG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "AGA" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG", "TGA" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
	// 25. Candidate Division SR1 and Gracilibacteria Code
	{
		{ "GCT", "GCC", "GCA", "GCG" }, // A
		{ "TGT", "TGC" }, // C
		{ "GAT", "GAC" }, // D
		{ "GAA", "GAG" }, // E
		{ "TTT", "TTC" }, // F
		{ "GGT", "GGC", "GGA", "GGG", "TGA" }, // G
		{ "CAT", "CAC" }, // H
		{ "ATT", "ATC", "ATA" }, // I
		{ "AAA", "AAG" }, // K
		{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" }, // L
		{}, // CodonTable::Leu1
		{ "ATG" }, // M
		{ "AAT", "AAC" }, // N
		{ "CCT", "CCC", "CCA", "CCG" }, // P
		{ "CAA", "CAG" }, // Q
		{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" }, // R
		{}, // CodonTable::Ser1
		{}, // CodonTable::Ser2
		{ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" }, // S
		{ "ACT", "ACC", "ACA", "ACG" }, // T
		{}, // CodonTable::Thr4_1
		{}, // CodonTable::Thr4_2
		{ "GTT", "GTC", "GTA", "GTG" }, // V
		{ "TGG" }, // W
		{ "TAT", "TAC" }, // Y
		{ "TAA", "TAG" }, // X
	},
};

void CodonTable::createCodonTable(unsigned tableId, std::string model, bool split, std::vector <std::string> groupList)
{

    codonTable = new CodonTable(tableId, model, split, groupList);
    codonTable -> setupCodonTable();
}


CodonTable* CodonTable::getInstance()
{
    return codonTable;
}