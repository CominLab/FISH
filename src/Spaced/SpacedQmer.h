/*
 * SpacedQmer.h
 *
 *  Created on: 20/lug/2016
 *      Author: samuele
 */

#ifndef SRC_HASH_SPACEDQMER_H_
#define SRC_HASH_SPACEDQMER_H_
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

//For unit
struct Pos_Ones {
	size_t index_one = numeric_limits<size_t>::max();
	size_t n_one = 0;
	size_t pos_start = 0;
	size_t n_one_before = 0;//numero di 1 che precedono il primo 1 dello unit
};
typedef vector<Pos_Ones> V_Pos_Ones;
struct Unit {
	vector<size_t> n_one;
	V_Pos_Ones v_pos;
};

//for previous
typedef vector<size_t> Position;
struct PreviousShift {
	Position one_to_change;
	Position one_to_remove;
	size_t one_exit = 0;
	size_t shift_min = 0;
	inline size_t GetSize() const {
		return this->one_to_change.size() + this->one_to_remove.size() + this->one_exit + (this->one_to_remove.empty() ? 0:2);
	}
};
typedef vector<PreviousShift> V_PreviusShift;

class SpacedQmer {
public:
	SpacedQmer();
	SpacedQmer(string spaced_qmer);
	inline size_t GetWeight() const {
		return this->pos_one.size();
	}
	inline size_t GetQ() const {
		return this->spaced_q.length();
	}
	inline bool isOne(size_t index) const {
		return this->spaced_q[index] == '1';
	}
	inline const Position& GetPosOne() const {
		return this->pos_one;
	}
	inline const Unit& GetUnit() const {
		return this->unit_map;
	}
	inline const string& toString() const {
		return this->spaced_q;
	}
	void reset(string spaced_qmer);
private:
	string spaced_q;
	Position pos_one;
	Unit unit_map;
	V_PreviusShift shift_min_change;

	void SaveIndexOne();
	void GetUnit(Unit& v_unit);
};

#endif /* SRC_HASH_SPACEDQMER_H_ */
