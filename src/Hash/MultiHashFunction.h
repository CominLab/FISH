/*
 * MultiHashFunction.h
 *
 *  Created on: 10/feb/2017
 *      Author: samuele
 */

#ifndef HASH_MULTIHASHFUNCTION_H_
#define HASH_MULTIHASHFUNCTION_H_

#include "HashFunction.h"

inline static void GetHashes_speedup_multi_unit(const string& s_Str, const SpacedQmer_Multi& spaced_qmers,
		Hash_Err_V_V& vHashes, hash_type (*fConvertion)(char)) {
	//Inizializzo vettore per gli hash e errori
//	if(vHashes.size() != v_spaced.size())
//	{
//		vHashes.resize(v_spaced.size());
//		vHashes.shrink_to_fit();
//	}

	//Get hash v for all unit present
	const MapUnit& map_unit = spaced_qmers.getMapUnit();
	Hash_Err_V_V hash_v(map_unit.n_one.size());
	#pragma omp parallel for
	for(size_t i = 0; i < map_unit.n_one.size(); ++i)//parallel computation
		GetHashes_speedup_previous(s_Str, map_unit.n_one[i], hash_v[i], fConvertion);

	#pragma omp parallel for
	for(size_t s = 0; s < spaced_qmers.size(); ++s)
	{
		const SpacedQmer& spaced = spaced_qmers[s];
		long n_hashes = s_Str.size()-spaced.GetQ()+1;

		Hash_Err_V& hashes = vHashes[s];
		hashes.clear();
		if(n_hashes>0)
		{
			hashes.resize(n_hashes); //Crea vettore

			//Combine different hash
			const V_Pos_Ones& v_pos = map_unit.v_v_pos[s];
			#pragma omp parallel for
			for(size_t i = 0; i < hashes.size(); ++i)
			{
				Hash_Err& curr_hash = hashes[i];
				for(const Pos_Ones& unit_pos : v_pos)
				{
					Hash_Err_V& hash_unit_v = hash_v[unit_pos.index_one];

					Hash_Err& hash_unit = hash_unit_v[i+unit_pos.pos_start];

					curr_hash.hash |= (hash_unit.hash << unit_pos.n_one_before*2);
					if(!hash_unit.isCorrect())
						curr_hash.add_pos_err(unit_pos.n_one_before, hash_unit);//aggiungi errore posizione corretta
				}
			}
		}
	}
}

#endif /* HASH_MULTIHASHFUNCTION_H_ */
