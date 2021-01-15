//============================================================================
// Name        : USH.cpp
// Author      : Samuele Girotto
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "FISH/FISH.h"

int main(int argc, char* argv[]) {
	string dir_output = "output/";
	FileParameter param;

	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i], "-si") == 0)
		{
			i++;
			if(!param.init(argv[i], ""))
			{
				cerr<<endl<<"Please enter an input filename single-end: -si <AbsPathFile>"<<flush;
				return 0;
			}
		}
		else if(strcmp(argv[i], "-pi") == 0)
		{
			i++;
			if(!param.init(argv[i], argv[i+1]))
			{
				cerr<<endl<<"Please enter an input filenames paired-end: -pi <AbsPathFile1> <AbsPathFile2>"<<flush;
				return 0;
			}
		}
		else if(strcmp(argv[i], "-q") == 0)
		{
			i++;
			vector<string> lines;
			string pathQmers(argv[i]);
			if(getLines(pathQmers, lines))
			{
				vector<bool> correctQmer(lines.size(), false);
				regex rgx("^1(0|1)*1$");
				for(size_t i = 0; i < lines.size(); i++)
				{
					correctQmer[i] = regex_match(lines[i], rgx);
					if(!correctQmer[i])
					{
						cerr<<endl<<"Error on " << to_string(i+1) << "° spaced seed. Enter q-mer with 1 at begin and end of the string on input files. "
						"Ex. 1**1*11*1. 1 is the simbol considered, any others are not valid simbols"<<flush;
						return 0;
					}
					else
						param.addSpacedQmer(lines[i], lines[i]);
				}
			}
			else
			{
				cerr<<endl<<"Please enter a spaced seeds path as -q <AbsPathFile>. Every file's line must contain a spaced seeds."<<flush;
				return 0;
			}
		}
		else if(strcmp(argv[i], "-dirO") == 0)
		{
			i++;
			if(string(argv[i]) == "")
			{
				cerr<<endl<<"Enter valid path for output directory"<<flush;
				return 0;
			}
			else
				dir_output = argv[i];
		}
	}
	if(param.getVSpaced().empty())//Applied default
	{
		param.addSpacedQmer("CLARK-S", "1111011101110010111001011011111");
		param.addSpacedQmer("CLARK-S", "1111101011100101101110011011111");
		param.addSpacedQmer("CLARK-S", "1111101001110101101100111011111");
		param.addSpacedQmer("rasbhari_minimizing_overlap_complexity", "1111010111010011001110111110111");
		param.addSpacedQmer("rasbhari_minimizing_overlap_complexity", "1110111011101111010010110011111");
		param.addSpacedQmer("rasbhari_minimizing_overlap_complexity", "1111101001011100111110101101111");
		param.addSpacedQmer("rasbhari_maximizing_sensitivity", "1111011110011010111110101011011");
		param.addSpacedQmer("rasbhari_maximizing_sensitivity", "1110101011101100110100111111111");
		param.addSpacedQmer("rasbhari_maximizing_sensitivity", "1111110101101011100111011001111");
		cout << endl << "Applied default spaced seed" << flush;
	}
	//Creo cartella output se non presente
	createDirAndSubDir(dir_output);

	cout << endl << "Applied the following spaced seed..." << flush;
	for(size_t i = 0; i < param.getVSpaced().size(); i++)
		cout << endl << "Type:" << param.getVSpaced()[i].first << ", Spaced seed: " << param.getVSpaced()[i].second.toString() << flush;

	bool single_test_equals = false;
	string dir_output_1 = dir_output + "single/";
	FISH test_single;
	if(test_single.load_sequences(param))
	{
		for(size_t j = 0; j < param.getVSpaced().size(); ++j)
		{
			string dir_output_tmp = dir_output_1 + param.getVSpaced()[j].first +  param.getVSpaced()[j].second.toString() + "/";

			test_single.single_run(param.getVSpaced()[j].second, single_test_equals);
			test_single.single_save(param, dir_output_tmp);
		}
	}

	////////////////////////////////////////////////////////////////7

	bool multi_test_equals = false;
	string dir_output_2 = dir_output + "multi/";
	vector<SpacedQmer> multi_spaced;
	for(size_t i = 0; i < param.getVSpaced().size(); ++i)
		multi_spaced.push_back(param.getVSpaced()[i].second);

	FISH test_multi;
	if(test_multi.load_sequences(param))
	{
		test_multi.multi_run(multi_spaced, multi_test_equals);
		test_multi.multi_save(param, multi_spaced, dir_output_2);
	}
}
