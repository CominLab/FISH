#include "SpacedQmer.h"
#include "../Hash/HashType.h"

//For multi unit
struct MapUnit {
	vector<size_t> n_one;
	vector<V_Pos_Ones> v_v_pos;
};

inline static void GetMapUnit(const vector<SpacedQmer>& v_spaced, MapUnit& map_unit) {
	map_unit.v_v_pos.resize(v_spaced.size());
	for(size_t i = 0; i < v_spaced.size(); ++i)
	{
		const Unit& unit = v_spaced[i].GetUnit();
		for(size_t j = 0; j < unit.n_one.size(); ++j)
			map_unit.n_one.push_back(unit.n_one[j]);
		for(size_t j = 0; j < unit.v_pos.size(); ++j)
			map_unit.v_v_pos[i].push_back(unit.v_pos[j]);
	}
	//order and remove duplicate
	sort(map_unit.n_one.begin(), map_unit.n_one.end());
	map_unit.n_one.resize(distance(map_unit.n_one.begin(), unique(map_unit.n_one.begin(), map_unit.n_one.end())));

	//order and remove duplicate
	//assign index
	#pragma omp parallel for
	for(size_t i = 0; i < map_unit.v_v_pos.size(); ++i)
		#pragma omp parallel for
		for(size_t j = 0; j < map_unit.v_v_pos[i].size(); ++j)
			map_unit.v_v_pos[i][j].index_one = distance(map_unit.n_one.begin(), find(map_unit.n_one.begin(), map_unit.n_one.end(), map_unit.v_v_pos[i][j].n_one));
}

class SpacedQmer_Multi {
public:
	inline SpacedQmer& operator[](size_t i){return this->v_spaced[i];}
	inline const SpacedQmer& operator[](size_t i) const{return this->v_spaced[i];}
    inline size_t size() const {return this->v_spaced.size();}
	inline const MapUnit& getMapUnit() const {return map_unit;}

	inline void init(const vector<SpacedQmer>& v_spaced) {
		this->v_spaced = v_spaced;
		GetMapUnit(this->v_spaced, this->map_unit);
	}

private:
	vector<SpacedQmer> v_spaced;
	MapUnit map_unit;
};
