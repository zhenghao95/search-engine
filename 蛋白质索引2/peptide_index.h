#ifndef  _PEPTIDE_INDEX_H_
#define	 _PEPTIDE_INDEX_H_

#include"file_operation.h"
#include"data_struct.h"
using namespace std;

int create_peptide_index();				//�����Ķ�����
void peptide_query();						//�Ķβ�ѯ
	
int Digestion();							//�����е�����ø��
int specific_enzyme(const char*str, const bool Mflag, Peptide_Dictionary* pep_dic, int & pep_dic_len, const Peptide_Dictionary*basic_dic, const int basic_len);	//����ø��
int  create_basic_pepDic(const char enzyme_teminal, const int proID, const char*sq_str, const double aa[], Peptide_Dictionary* basic_dic);			//����һ��sq,��ÿ��ø�е�ø�У�����basic_dic
int k_merge(int ls[], Peptide_Dictionary b[], int k);	//�鲢ȥ����


bool read_sq(const size_t filenum, size_t&file_id, FILE*&index_fp, FILE*&sq_fp, char*sq_str); //��sq�ļ��ж�ȡһ��sq����
int ans_enzyme(const char*sq_str);													//ͳ��sq��ø��λ�����Ŀ
static bool isBelong(const char a);													//�ж�ĳ���ַ��Ƿ���ø�е�
static int create_prime_array();													//������������
static double Cal_Godel_value(const char*str, const int len);						//�����¶����룬�����Ķε���ʼλ�úͳ���
void printf_pepDic(const Peptide_Dictionary *Dic, const int dic_len);				//����Ķδʵ���ȫ���Ķ���Ϣ��Peptide_Dictionary��
void printf_pepDicFinal(const Peptide_Dictionary_FINAL*dic, const int dic_len);		//����Ķδʵ���ȫ���Ķ���Ϣ(Peptide_Dictionary_FINAL)



static void sort_word(Peptide_Dictionary *array, size_t n, int word, Peptide_Dictionary*buf);	//��double���ĸ��ֽ�����
void fsort_mass(Peptide_Dictionary *array, size_t n);											//��������
static void adjust(int ls[], Peptide_Dictionary b[], int k, int s);								//����������
static void create_loser_tree(int ls[], Peptide_Dictionary b[], int k);							//����������


#endif // ! _PEPTIDE_INDEX_H_