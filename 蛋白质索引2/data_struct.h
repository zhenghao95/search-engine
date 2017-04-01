#ifndef _DATA_STRUCT_H_
#define _DATA_STRUCT_H_

#include<utility>
using namespace std;


#define H2O 18
const int  max_mem_size = 1024 * 1024 *512;
const int  max_sq_len = 50000;

typedef pair<size_t, size_t>P;

struct Protein_Meta
{
	size_t offset = 0;
	size_t Pro_Num = 0;
	unsigned char File_Num = 0;
};

struct Protein_Index_Table
{
	size_t AC_Pos = 0;
	size_t ac_Len = 0;
	size_t DE_Pos = 0;
	size_t DE_Len = 0;
	size_t SQ_Pos = 0;
	size_t SQ_Len = 0;
};

struct Peptide_Dictionary
{
	size_t Pos = 0;
	int    Len = 0;
	double  Mass = 0;
	int    MissCleaveSites = 0;
	size_t InvertFilePos = 0;
	double Godel = 0.0;		
	Peptide_Dictionary& operator=(Peptide_Dictionary& value)
	{
		Pos = value.Pos;
		Len = value.Len;
		Mass = value.Mass;
		MissCleaveSites = value.MissCleaveSites;
		InvertFilePos = value.InvertFilePos;
		Godel = value.Godel;
		return *this;
	}

	bool operator>(const Peptide_Dictionary b) const			//按质量排序
	{
		return this->Mass > b.Mass;
	}
};

struct Peptide_Dictionary_FINAL
{
	size_t Pos = 0;
	int    Len = 0;
	double  Mass = 0;
	int    MissCleaveSites = 0;
	fpos_t InvertFilePos = 0;// 用于存放位置
	
};

struct Peptide_Mass_Dictionary

{
	double Mass =0;
	size_t Pep_ID=0;
};

#endif // !_DATA_STRUCT_H_


