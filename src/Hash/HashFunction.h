/*
 * HashFunction.h
 *
 *  Created on: 20/lug/2016
 *      Author: samuele
 */

#ifndef HASHFUNCTION_H_
#define HASHFUNCTION_H_

#include "HashType.h"
#include "../Spaced/SpacedQmer_Multi.h"
#include <algorithm>
#include <iostream>

inline static hash_type CharToInt(char ch)
{
	if(ch == 'A')
		return 0;
	if(ch == 'C')
		return 1;
	if(ch == 'G')
		return 2;
	if(ch == 'T')
		return 3;
	return 4; //ERROR CODE
}

inline static hash_type CharToIntComplement(char ch)
{
	if(ch == 'A')
		return 3;
	if(ch == 'C')
		return 2;
	if(ch == 'G')
		return 1;
	if(ch == 'T')
		return 0;
	return 4; //ERROR CODE
}

//Hash per tutti 1 su spaced qmer
inline static void GetHash(const string& s_Str, size_t startQmer, size_t length, Hash_Err& hash_err, hash_type (*fConvertion)(char))
{
	hash_err.reset();
//	#pragma omp parallel for ordered
	for(size_t i = startQmer; i < startQmer + length; ++i)
	{
		hash_type ch = (*fConvertion)(s_Str[i]);
//		#pragma omp ordered
		if(ch == 4) //Errore conversione
			hash_err.push_back_error(i);
		else
			hash_err.hash |= ch << ((i - startQmer) * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
	}
}

//Hash per spaced qmer con *
inline static void GetHash(const string& s_Str, size_t startQmer, const SpacedQmer& spaced_qmer,
		Hash_Err& hash_err, hash_type (*fConvertion)(char))
{
	hash_err.reset();
	const Position& pos_one = spaced_qmer.GetPosOne();
	for(size_t j = 0; j < pos_one.size(); ++j)
	{
		hash_type ch = (*fConvertion)(s_Str[startQmer+pos_one[j]]);
		if(ch == 4) //Errore conversione
			hash_err.push_back_error(j);
		else
			hash_err.hash |= ch << (j * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
	}
}

//Hash veloce con spaced qmer tutti 1
inline static void GetHashes_speedup_previous(const string& s_Str, size_t length,
		Hash_Err_V& vHash, hash_type (*fConvertion)(char)) {
	vHash.clear();
	if(s_Str.size() >= length)
	{
		size_t n_hashes = s_Str.size() - length + 1;
		vHash.resize(n_hashes); //Crea vettore

		GetHash(s_Str, 0, length, vHash[0], fConvertion);//primo da computare a parte
		for(size_t pos=1; pos < vHash.size(); ++pos)
		{
			Hash_Err& prev_hash = vHash[pos-1];
			Hash_Err& curr_hash = vHash[pos];

			//copia hash e sottrai una posizione dal precedente
			curr_hash.hash = prev_hash.hash;
			curr_hash.hash >>= 2; //sposta 2 bit, esce una lettera
			curr_hash.sub_pos_err(1, prev_hash);

			hash_type enter = (*fConvertion)(s_Str[pos+length-1]);
			if(enter == 4)
				curr_hash.push_back_error(length-1);
			else
				curr_hash.hash |= enter << ((length - 1) * 2);	//aggiungi ultimo elemento OR possibile perchè prima ho
																//diviso per 4 e la posizione dove scrivo ha sicuramente 0
		}
	}
}

inline static void GetHashes_naive(const string& s_Str, const SpacedQmer& spaced_qmer,
		Hash_Err_V& vHash, hash_type (*fConvertion)(char))
{
//	bool isAllOne = spaced_qmer.GetWeight() == spaced_qmer.GetQ();
//	if(isAllOne)
//		GetHashes_speedup_previous(s_Str, spaced_qmer.GetQ(), vHash, fConvertion);
//	else
//	{
		vHash.clear();
		if(s_Str.size() >= spaced_qmer.GetQ())
		{
			size_t n_hashes = s_Str.size() - spaced_qmer.GetQ() + 1;
			vHash.resize(n_hashes); //Crea vettore
			#pragma omp parallel for
			for(size_t pos=0; pos < vHash.size(); ++pos)
				GetHash(s_Str, pos, spaced_qmer, vHash[pos], fConvertion);
		}
//	}
}

inline static void GetHashes_speedup_unit(const string& s_Str, const SpacedQmer& spaced_qmer,
		Hash_Err_V& vHash, hash_type (*fConvertion)(char)) {
//	bool isAllOne = spaced_qmer.GetWeight() == spaced_qmer.GetQ();
//	if(isAllOne)
//		GetHashes_speedup_previous(s_Str, spaced_qmer.GetQ(), vHash, fConvertion);
//	else
//	{
		vHash.clear();
		if(s_Str.size() >= spaced_qmer.GetQ())
		{
			const Unit& spaced_v = spaced_qmer.GetUnit();

			//Get hash v for all unit present
			//TODO: si può migliorare, dato che ordinate si può prendere quelle più piccole per comporre quelle più grandi
			vector<Hash_Err_V> hash_v(spaced_v.n_one.size());
			#pragma omp parallel for
			for(size_t i = 0; i < spaced_v.n_one.size(); ++i)//parallel computation
				GetHashes_speedup_previous(s_Str, spaced_v.n_one[i], hash_v[i], fConvertion);

			//Combine different hash
			size_t n_hashes = s_Str.size() - spaced_qmer.GetQ() + 1;
			vHash.resize(n_hashes); //Crea vettore

			const V_Pos_Ones& v_pos = spaced_v.v_pos;
			#pragma omp parallel for
			for(size_t i = 0; i < vHash.size(); ++i)
			{
				Hash_Err& curr_hash = vHash[i];
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
//	}
}

#endif /* HASHFUNCTION_H_ */
