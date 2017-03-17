#include<utility>
#define MAX_PRO_NUM 50

using namespace std;

typedef pair<size_t, size_t>P;

struct Protein_Meta
{
	size_t offset;
	size_t Pro_Num;
	unsigned char File_Num;
};

struct Protein_Index_Table
{
	size_t AC_Pos;
	size_t ac_Len;
	size_t DE_Pos;
	size_t DE_Len;
	size_t SQ_Pos;
	size_t SQ_Len;
};

struct Peptide_Dictionary
{
	size_t Pos;
	unsigned char Len;
	size_t Mass;
	unsigned char MissCleaveSites;
	size_t InvertFilePos;
};

struct Peptide_Mass_Dictionary
{
	size_t Mass;
	size_t Pep_ID;
};

struct Peptide_Protein_Invert
{
	size_t Pro_num;
	size_t Pro_ID[MAX_PRO_NUM];
	//记得设立错误判别机制
   //即如果一个pep对应的蛋白质数目超过了MAX_PRO_NUM
	//要报错

};
