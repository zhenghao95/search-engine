#ifndef  _PEPTIDE_INDEX_H_
#define	 _PEPTIDE_INDEX_H_

#include"file_operation.h"
#include"data_struct.h"
using namespace std;

int create_peptide_index();				//创建肽段索引
void peptide_query();						//肽段查询
	
int Digestion();							//对所有蛋白质酶切
int specific_enzyme(const char*str, const bool Mflag, Peptide_Dictionary* pep_dic, int & pep_dic_len, const Peptide_Dictionary*basic_dic, const int basic_len);	//特异酶切
int  create_basic_pepDic(const char enzyme_teminal, const int proID, const char*sq_str, const double aa[], Peptide_Dictionary* basic_dic);			//历遍一遍sq,在每个酶切点酶切，产生basic_dic
int k_merge(int ls[], Peptide_Dictionary b[], int k);	//归并去冗余


bool read_sq(const size_t filenum, size_t&file_id, FILE*&index_fp, FILE*&sq_fp, char*sq_str); //从sq文件中读取一条sq数据
int ans_enzyme(const char*sq_str);													//统计sq中酶切位点的数目
static bool isBelong(const char a);													//判断某个字符是否是酶切点
static int create_prime_array();													//创建素数数组
static double Cal_Godel_value(const char*str, const int len);						//计算哥德尔编码，传入肽段的起始位置和长度
void printf_pepDic(const Peptide_Dictionary *Dic, const int dic_len);				//输出肽段词典中全部肽段信息（Peptide_Dictionary）
void printf_pepDicFinal(const Peptide_Dictionary_FINAL*dic, const int dic_len);		//输出肽段词典中全部肽段信息(Peptide_Dictionary_FINAL)



static void sort_word(Peptide_Dictionary *array, size_t n, int word, Peptide_Dictionary*buf);	//对double的四个字节排序
void fsort_mass(Peptide_Dictionary *array, size_t n);											//基数排序
static void adjust(int ls[], Peptide_Dictionary b[], int k, int s);								//败者树调整
static void create_loser_tree(int ls[], Peptide_Dictionary b[], int k);							//创建败者树


#endif // ! _PEPTIDE_INDEX_H_