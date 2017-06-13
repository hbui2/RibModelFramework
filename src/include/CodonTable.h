#ifndef CodonTable_H
#define CodonTable_H

#include <string>
#include <map>
#include <algorithm>
#include <cctype>
#include <vector>
#include <array>
#include "Utility.h"

class CodonTable
{
    private:

		static std::vector <std::string> defaultVector;

        static CodonTable *codonTable;
        unsigned tableId;
        bool splitAA;
        //Stored by AA index then codon index. Excludes the last codon index in every AA grouping.

		std::map <std::string, std::vector <std::string>> AAToCodonMap;
		std::map <std::string, std::vector <unsigned>> AAToCodonIndexMap;
		std::map <std::string, std::vector <std::string>> AAToCodonMapWithoutReference;
		std::map <std::string, std::vector <unsigned>> AAToCodonIndexMapWithoutReference;
        std::map <std::string, std::string> codonToAAMap; //Maps ALL codons for current conditions to AAs.
        const static std::map <std::string, unsigned> AAMap; //Maps currently used AAs to indices.
        std::map <std::string, unsigned> AAToNumCodonsMap;
        //Maps currently used AAs to the number of codons that code for them.

        //Maps codons to indices (not including last in each AA grouping).
		std::vector <std::string> groupList;

		const static std::vector <std::vector <std::string> > defaultAAGroupListings;
		const static std::vector <std::vector <std::string> > defaultSplitAAGroupListings;

		const static std::vector <std::vector <std::vector <std::string> > > codonTableListing;
		const static std::vector <std::vector <std::vector <std::string> > > codonTableListingWithoutSplit;

    public:

        //Constructors & destructors:
        explicit CodonTable(); //Defaults to table 1 and splitting AA
        CodonTable(unsigned _tableId, std::string model, bool _splitAA = true, std::vector <std::string> _groupList = defaultVector);
        virtual ~CodonTable();
        CodonTable(const CodonTable& other); //Todo: Need? If so update the function.
        CodonTable& operator=(const CodonTable& other); //Todo: Need? if so update the function.



        //Getter functions:
        unsigned getTableId();
        bool getSplitAA();
		std::map <std::string, std::vector <std::string>> getAAToCodonMap();
		std::map <std::string, std::vector <unsigned>> getAAToCodonIndexMap();
		std::map <std::string, std::vector <std::string>> getAAToCodonMapWithoutReference();
		std::map <std::string, std::vector <unsigned>> getAAToCodonIndexMapWithoutReference();
        std::map <std::string, std::string> getCodonToAAMap(); //Maps ALL codons for current conditions to AAs.
		std::map <std::string, unsigned> getAAMap(); //Maps currently used AAs to indices.
		std::map <std::string, unsigned> getAAToNumCodonsMap();
		std::vector <std::string> getGroupList();
        
        //Mapping operations:
        unsigned AAToAAIndex(std::string aa);
        std::vector <unsigned> AAIndexToCodonRange(unsigned aaIndex, bool withoutReference = false);
        static std::string indexToCodon(unsigned index);
        std::vector <unsigned> AAToCodonRange(std::string aa, bool withoutReference = false);
        std::vector<std::string> AAToCodon(std::string aa, bool withoutReference = false);
        std::string codonToAA(std::string& codon);
		static unsigned codonToIndex(std::string& codon);
        unsigned codonToAAIndex(std::string& codon);
        std::string indexToAA(unsigned aaIndex);
		unsigned getNumCodonsForAA(std::string aa, bool withoutReference = false);
		unsigned getNumCodonsForAAIndex(unsigned aaIndex, bool withoutReference = false);



        //Other functions:
        void setupCodonTable(); //Sets up the private variables that do all the mappings.
        bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);


        //Static variables & functions:
        static const std::string Ser2;
        static const std::string Ser1; //Necessary for codon table 12
        static const std::string Thr4_1; //Necessary for codon table 3
        static const std::string Thr4_2; //Necessary for codon table 3
        static const std::string Leu1; //Necessary for codon table 16, 22


        static const std::vector <std::string> aminoAcidArray; //Index = AA
        static const std::vector <std::string> aminoAcidArrayWithoutSplit; //Array containing all non-split AAs.
		static const std::map<std::string, unsigned> fullAAMap; 
        static const std::map<std::string, unsigned> codonToIndexWithReference; //Map of indices to all codons.
        static const std::string codonArray[]; //List of codons.
        static const std::vector <std::string> codonTableDefinition;

        //Description title for each codon table according to NCBI.

        static const unsigned numCodonsPerAAForTable[25][26]; //Sized on tableId and AA.


        static void createCodonTable(unsigned tableId, std::string model, bool split = true, std::vector <std::string> groupList = defaultVector); //Used to create the singleton instance.
        static CodonTable* getInstance(); //Get the singleton instance.



        //--------------------R WRAPPERS--------------------//
        //Getter functions:
#ifndef STANDALONE
        unsigned getTableIdR();
        std::vector<std::vector<unsigned>> getCodonIndexListingR();
        std::vector<std::vector<unsigned>> getCodonIndexListingWithoutReferenceR();
        std::map <std::string, unsigned> getAAMapR();
        std::map <std::string, unsigned> getForParamVectorMapR();

        unsigned getNumCodonsForAAIndexR(unsigned aaIndex, bool forParamVector = false);
        std::string getForParamVectorCodonR(unsigned codonIndex);



        //Mapping operations:
        unsigned AAToAAIndexR(std::string aa);
        std::vector <unsigned> AAIndexToCodonRangeR(unsigned aaIndex, bool forParamVector = false);
        std::string indexToCodonR(unsigned index);
        std::vector <unsigned> AAToCodonRangeR(std::string aa, bool forParamVector = false);
        std::vector<std::string> AAToCodonR(std::string aa, bool forParamVector = false);
        std::string codonToAAR(std::string& codon);
        unsigned codonToIndexR(std::string& codon);
        unsigned codonToAAIndexR(std::string& codon);
        std::string indexToAAR(unsigned aaIndex);



        //Static getter functions:
        static std::string getSer2R();
        static std::string getSer1R();
        static std::string getThr4_1R();
        static std::string getThr4_2R();
        static std::string getLeu1R();

        static std::vector <std::string> getAminoAcidArrayR();
        static std::vector <std::string> getAminoAcidArrayWithoutSplitR();
        static std::vector <std::vector <unsigned>> getNumCodonsPerAAForTableR();
        static std::vector <std::string> getCodonTableDefinitionR();

        static std::vector<std::string> getCodonArrayR();

#endif

};

#endif
