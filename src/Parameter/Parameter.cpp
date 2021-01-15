/*
 * Parameter.cpp
 *
 *  Created on: 05/feb/2016
 *      Author: samuele
 */

#include "Parameter.h"

FileParameter::FileParameter() {

}

bool FileParameter::init(string pathFile1, string pathFile2) {
	return this->input_files.init(pathFile1, pathFile2);
}

void FileParameter::addSpacedQmer(string nametype, string spaced) {
	this->spaced_seed.push_back(make_pair(nametype+"/", SpacedQmer(spaced)));
}

const PairFiles& FileParameter::getInputFiles() const {
	return this->input_files;
}

const vector<pair<string, SpacedQmer>>& FileParameter::getVSpaced() const {
	return this->spaced_seed;
}

