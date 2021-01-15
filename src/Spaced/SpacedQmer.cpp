/*
 * SpacedQmer.cpp
 *
 *  Created on: 20/lug/2016
 *      Author: samuele
 */

#include "SpacedQmer.h"

SpacedQmer::SpacedQmer()
{
	this->reset("");
}

SpacedQmer::SpacedQmer(string spaced_qmer) {
	this->reset(spaced_qmer);
}

void SpacedQmer::reset(string spaced_qmer) {
	this->spaced_q = spaced_qmer;
	this->SaveIndexOne();
	this->GetUnit(this->unit_map);
}

void SpacedQmer::SaveIndexOne() {
	this->pos_one.clear();this->pos_one.shrink_to_fit();
	for(size_t i = 0; i < this->spaced_q.length(); ++i)
		if(this->isOne(i))
			this->pos_one.push_back(i);
}

void SpacedQmer::GetUnit(Unit& v_unit) {
	//reset
	v_unit.n_one.clear();v_unit.n_one.reserve(this->pos_one.size()+1);//Dimensione massima è il weight+1
	v_unit.v_pos.clear();v_unit.v_pos.reserve(this->pos_one.size()+1);//Dimensione massima è il weight+1
	//Calcola
	size_t init_unit = 0;
	size_t size_unit = 0;
	for(size_t i = 0; i < this->pos_one.size(); ++i)
	{
		if(this->pos_one[i] != this->pos_one[init_unit]+size_unit)
		{
			Pos_Ones new_unit;
			new_unit.n_one = size_unit;
			new_unit.pos_start = this->pos_one[init_unit];
			new_unit.n_one_before = init_unit;
			v_unit.v_pos.push_back(move(new_unit));

			v_unit.n_one.push_back(size_unit); //inserisco su vettore, poi elimino quindi perdo posizione

			init_unit = i;
			size_unit = 0;
		}
		++size_unit;
	}
	if(size_unit > 0)
	{
		Pos_Ones new_unit;
		new_unit.n_one = size_unit;
		new_unit.pos_start = this->pos_one[init_unit];
		new_unit.n_one_before = init_unit;
		v_unit.v_pos.push_back(move(new_unit));

		v_unit.n_one.push_back(size_unit); //inserisco su vettore, poi elimino quindi perdo posizione
	}
	//order and remove duplicate
	sort(v_unit.n_one.begin(), v_unit.n_one.end());
	v_unit.n_one.resize(distance(v_unit.n_one.begin(), unique(v_unit.n_one.begin(), v_unit.n_one.end())));

	//assign index
	for(size_t i = 0; i < v_unit.v_pos.size(); ++i)
		v_unit.v_pos[i].index_one = distance(v_unit.n_one.begin(), find(v_unit.n_one.begin(), v_unit.n_one.end(), v_unit.v_pos[i].n_one));
}
