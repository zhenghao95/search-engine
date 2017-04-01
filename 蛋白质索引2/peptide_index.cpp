#include"peptide_index.h"
#include<cmath>
#include <cstdlib>
#include <cstring>
#include <cstdbool>

typedef unsigned long long uint64;
typedef long long int64;
typedef unsigned short uint16;

/*�ⲿȫ�ֱ���-from file_operation.c */
extern char output_dir_result[300];
extern char enzyme_aa[5];
extern char enzyme_terminal;
extern double min_mass, max_mass;
extern int min_pep_length, max_pep_length, max_miss_site;



/*ȫ�ֱ���*/
const int max_ezy_num = 800;		//���ø��λ����
#define prime_num 120				//��������Ŀ��120������
#define prime_range 600				//����600���ڵ�����
double log_prime[prime_num];			
const int max_pepDic_size = max_mem_size / sizeof(Peptide_Dictionary);
const int max_pepDicFINALsize = max_mem_size / sizeof(Peptide_Dictionary_FINAL);  //ÿ���Ķδʵ��ļ����Ķε�����
const int max_same_pep = 200;		//	����ͬ�������Ķ����������

double aa_mass[26] = { 71.03711,166.99836,103.00919,		//������������
115.02694,129.04259,147.06841,57.02146,
137.05891,113.08406,181.01401,128.09496,
113.08406,131.04048,114.04293,243.02965,
97.05276,128.05858,156.10111,87.03203 ,
101.04768 ,0.00000,99.06841,186.07931,
113.08406,163.06332,128.55059 };



/*
Function:		create_peptide_index
Description:    �����Ķ�����
Calls:			create_prime_array
				Digestion
				k_merge
Input:			��
Return:         �Ķδʵ��ļ���
*/
int create_peptide_index()
{
	create_prime_array();

	int tmpFileNum = Digestion();	
	printf("��ʱ�ļ���Ϊ��%d\n", tmpFileNum);

	int *ls = new int[tmpFileNum+3];
	Peptide_Dictionary*b = new Peptide_Dictionary[tmpFileNum + 3];
	
	int dic_num = k_merge(ls, b, tmpFileNum);		
	//printf("dic�ļ���Ϊ��%d\n", dic_num);

	delete[]ls;
	delete[]b;

	//return dic_num;
	return 0;
}



/*
Function:		Digestion
Description:    �����ݿ��е����е�����ø��
Calls:			read_sq
				create_basic_pepDic
				specific_enzyme
				fsort_mass
Input:			��
Return:         ������ʱ�ļ�pepTmp����Ŀ
*/
int Digestion()
{
	Protein_Meta pro_meta;
	size_t file_id = 0;
	FILE*index_fp, *sq_fp, *meta_fp;
	char*sq_str = new char[max_sq_len+1];

	meta_fp = fopen(new_file_name("pro", ".meta", output_dir_result).c_str(), "rb");
	check_fp(meta_fp);

	fread(&pro_meta, sizeof(Protein_Meta), 1, meta_fp);
	fclose(meta_fp);

	char*suffix = new char[8];
	sprintf(suffix, "%d", ++file_id);
	index_fp = fopen(new_file_name("proIndex", suffix, output_dir_result).c_str(), "rb");
	check_fp(index_fp);

	sq_fp = fopen(new_file_name("sq", suffix, output_dir_result).c_str(), "rb");
	check_fp(sq_fp);
	delete[]suffix;

	int proID = -1;
	Peptide_Dictionary *basic_dic = new Peptide_Dictionary[max_ezy_num];

	Peptide_Dictionary*Container = new Peptide_Dictionary[max_pepDic_size];
	int container_pos = -1;

	int tmpFileID = -1;			//��ʱ�ļ�ָ��
	char suff[6];
	FILE*fp_tmp;

	/*test*/
	int total_pep_num = 0;	//�Ķ�����

	while (read_sq(pro_meta.File_Num, file_id, index_fp, sq_fp, sq_str))
	{
		
		proID++;
		//printf("proID:%d	%s\n",proID, sq_str);
		int basic_len = create_basic_pepDic(enzyme_terminal, proID, sq_str, aa_mass, basic_dic);//�������ʵ�

		if (container_pos + (1 + basic_len)*(max_miss_site + 1) + 1 > max_pepDic_size)//����˴�ø��������ķ��գ�����д����ʱ�ļ�
		{
		//	printf("����ǰ�Ľ��Ϊ��\n");
		//	printf_pepDic(Container, container_pos + 1);
			//����
			fsort_mass(Container, container_pos + 1);
		//	printf("�����Ľ��Ϊ��\n");
			//printf_pepDic(Container, container_pos + 1);
			

			//д����ʱ�ļ�
			sprintf(suff, "%d", ++tmpFileID);
			fp_tmp = fopen(new_file_name("pepTmp", suff, output_dir_result).c_str(), "wb");//��һ����ʱ�ļ�
			fwrite(Container, sizeof(Peptide_Dictionary), container_pos + 1, fp_tmp);
			
			fclose(fp_tmp);
			//�������
			container_pos = -1;
		}

		bool mFlag = 0;//�ж��Ƿ���Ҫ����M��
		if (enzyme_terminal == 'N'&&sq_str[0] == 'M')mFlag = 1;
		int pepNum = specific_enzyme(sq_str, mFlag,  Container, container_pos, basic_dic, basic_len);	//����ø��
		/*test*/
		total_pep_num += pepNum;
		//printf("proID: %d cut pep num is ��%d \n", proID,pepNum);
	}

	/*test*/
	printf("after digestion, total pep num is:%d\n", total_pep_num);

	//����ʣ����Ķ�
	if (container_pos != -1)
	{
		//printf("����ǰ�Ľ��Ϊ��\n");
		//printf_pepDic(Container, container_pos + 1);
	
		fsort_mass(Container, container_pos + 1);
	//	printf("�����Ľ��Ϊ��\n");
		//printf_pepDic(Container, container_pos + 1);

		sprintf(suff, "%d", ++tmpFileID);
		fp_tmp = fopen(new_file_name("pepTmp", suff, output_dir_result).c_str(), "wb");
		check_fp(fp_tmp);

		fwrite(Container, sizeof(Peptide_Dictionary), container_pos + 1, fp_tmp);
		fclose(fp_tmp);

		container_pos = -1;
	}

	delete[]Container;
	delete[]basic_dic;
	delete[]sq_str;

	return tmpFileID + 1;          //������ʱ�ļ�����
}

/*
Function:		k_merge
Description:    ��ø���Ķι鲢ȥ����
Calls:			create_loser_tree
Input:			ls��Ԫ��λ��
b�����Ԫ��
k:�鲢��Ŀ
Return:         int�ͣ��Ķδʵ��ļ���
*/
int k_merge(int ls[], Peptide_Dictionary b[], int k)
{

	FILE* *fp_in = new FILE*[k];		//��ʱ�ļ�
	char suff[5];
	for (int i = 0; i < k; i++)
	{
		sprintf(suff, "%d", i);
		fp_in[i] = fopen(new_file_name("pepTmp", suff, output_dir_result).c_str(), "rb");
		check_fp(fp_in[i]);
	}

	int dic_ID = -1;
	sprintf(suff, "%d", ++dic_ID);
	FILE*dic_fp = fopen(new_file_name("pepDictionary", suff, output_dir_result).c_str(), "wb");//�Ķδʵ��ļ�
	check_fp(dic_fp);

	FILE* Ivt_fp = fopen(new_file_name("index", "Invert", output_dir_result).c_str(), "wb");	//���������ļ�
	check_fp(Ivt_fp);

	Peptide_Dictionary_FINAL  *dic = new Peptide_Dictionary_FINAL[max_pepDicFINALsize];
	int dic_pos = -1;

	for (int i = 0; i < k; ++i)
		fread(&b[i], sizeof(Peptide_Dictionary), 1, fp_in[i]);

	create_loser_tree(ls, b, k);		//����������
	int q = ls[0];
	Peptide_Dictionary f = b[q];		//�õ��Ѷ���Ԫ�أ�����f

	int*	pro_id = new int[max_same_pep];				    //��ŵ�����ID������
	int pro_id_pos = -1;
	pro_id[++pro_id_pos] = f.InvertFilePos;	//����f��proID

	fread(&b[q], sizeof(Peptide_Dictionary), 1, fp_in[q]);		//����һ������
	create_loser_tree(ls, b, k);		//����������

	/*test*/
	int pep_num_total = 0;

	while (b[ls[0]].Mass != (numeric_limits<double>::max)())
	{
		q = ls[0];
		Peptide_Dictionary l = b[q];

		if (l.Mass == f.Mass && l.Godel == f.Godel)			//�����ͬһ���Ķ�
		{
			if (l.InvertFilePos != f.InvertFilePos)
			{
				_ASSERTE(_CrtCheckMemory());
				pro_id[++pro_id_pos] = l.InvertFilePos;
				_ASSERTE(_CrtCheckMemory());
			}
		}
		else                          //���˵l,f����ͬһ���Ķ�
		{
			//�õ����ű��λ��
			fpos_t pos;
			fgetpos(Ivt_fp, &pos);

			//��pro_id_pos,pro_id����д�뵹�ű�
			int proNum = pro_id_pos + 1;
			fwrite(&proNum, sizeof(int), 1, Ivt_fp);
			fwrite(pro_id, sizeof(int), pro_id_pos + 1, Ivt_fp);
			pro_id_pos = -1;

			//���dicд���ˣ���dicд���ļ���
			if (dic_pos == max_pepDicFINALsize - 1)
			{
				/*test*/
				pep_num_total += dic_pos + 1;

				//�����������ȥ����������
				//printf("�鲢����ȥ����������Ϊ��\n");
				//printf_pepDicFinal(dic, max_pepDicFINALsize);

				fwrite(dic, sizeof(Peptide_Dictionary_FINAL), max_pepDicFINALsize, dic_fp);
				dic_pos = -1;
				fclose(dic_fp);

				sprintf(suff, "%d", ++dic_ID);
				dic_fp = fopen(new_file_name("pepDictionary", suff, output_dir_result).c_str(), "wb");//�Ķδʵ��ļ�
				check_fp(dic_fp);
			}
			//��dic�����һ������
			_ASSERTE(_CrtCheckMemory());
			dic[++dic_pos].InvertFilePos = pos;
			dic[dic_pos].Mass = f.Mass;
			dic[dic_pos].Len = f.Len;
			dic[dic_pos].MissCleaveSites = f.MissCleaveSites;
			dic[dic_pos].Pos = f.Pos;
			_ASSERTE(_CrtCheckMemory());

			//��f����Ϊl
			f = l;
			pro_id[++pro_id_pos] = f.InvertFilePos;
		}

		if (!fread(&b[q], sizeof(Peptide_Dictionary), 1, fp_in[q]))		//����һ������
			b[q].Mass = (numeric_limits<double>::max)();

		adjust(ls, b, k, q);				//����������
	}

	//�������һ������
	fpos_t pos;
	fgetpos(Ivt_fp, &pos);

	int proNum = pro_id_pos + 1;
	fwrite(&proNum, sizeof(int), 1, Ivt_fp);
	fwrite(pro_id, sizeof(int), pro_id_pos + 1, Ivt_fp);
	pro_id_pos = -1;

	if (dic_pos == max_pepDicFINALsize - 1)
	{
		/*test*/
		pep_num_total += dic_pos + 1;
		//�����������ȥ����������
		//printf("�鲢����ȥ����������Ϊ��\n");
		//printf_pepDicFinal(dic, max_pepDicFINALsize);

		fwrite(dic, sizeof(Peptide_Dictionary_FINAL), max_pepDicFINALsize, dic_fp);
		dic_pos = -1;
		fclose(dic_fp);

		sprintf(suff, "%d", ++dic_ID);
		dic_fp = fopen(new_file_name("pepDictionary", suff, output_dir_result).c_str(), "wb");//�Ķδʵ��ļ�
		check_fp(dic_fp);
	}
	_ASSERTE(_CrtCheckMemory());
	dic[++dic_pos].InvertFilePos = pos;
	dic[dic_pos].Mass = f.Mass;
	dic[dic_pos].Len = f.Len;
	dic[dic_pos].MissCleaveSites = f.MissCleaveSites;
	dic[dic_pos].Pos = f.Pos;
	_ASSERTE(_CrtCheckMemory());


	//��ʣ�������д���Ķ��ļ�,���ر��Ķ��ļ�
	if (dic_pos > -1)
	{
		/*test*/
		pep_num_total += dic_pos + 1;

		fwrite(dic, sizeof(Peptide_Dictionary_FINAL), dic_pos + 1, dic_fp);
		//printf("�鲢�������Ϊ��\n");
		//printf_pepDicFinal(dic, dic_pos + 1);
		dic_pos = -1;
	}

	/*test*/
	printf("ȥ�������Ķ�����%d\n", pep_num_total);

	fclose(dic_fp);


	//�ر���ʱ�ļ�
	for (int i = 0; i < k; i++)
		fclose(fp_in[i]);
	//ɾ����ʱ�ļ�
	for (int i = 0; i < k; i++)
	{
		sprintf(suff, "%d", i);
		if (remove(new_file_name("pepTmp", suff, output_dir_result).c_str()) != 0)
			perror("Error deleting file!");
	}
	//�رյ��ű��ļ�
	fclose(Ivt_fp);

	//�ͷ�����ռ�
	delete[]fp_in;
	delete[]pro_id;
	return dic_ID + 1;		//�����Ķδʵ��ļ���Ŀ
}
/*
Function:		specific_enzyme
Description:    ��һ��sq������ø��
Calls:			��
Input:			str��ָ��Ҫø�е�sq����
				M_flag:�Ƿ������Ҫ����N��M�����
				pep_dic�����ø�н��������
				pep_dic_pos��ָʾpep_dic��λ��
				basic_dic������ø�в������Ķλ��ʵ�
				basic_len��ָʾbasic_dic�ĳ���
Return:        int�ͣ��˴�ø���ܹ��������Ķ���Ŀ
*/

/*test*/
void printf_dic_pepStr(Peptide_Dictionary* pep_dic, int  pep_dic_pos,const char*str)
{
	for (int i = 0; i <= pep_dic_pos; i++)
	{
		for (int j = pep_dic[i].Pos; j < pep_dic[i].Pos + pep_dic[i].Len; j++)
			printf("%c", str[j]);
		printf("\n");
	}
}
int specific_enzyme(const char*str,const bool M_flag,Peptide_Dictionary* pep_dic, int & pep_dic_pos, const Peptide_Dictionary*basic_dic, const int basic_len)
{
	int start_pos = pep_dic_pos+1;
	
	int proID = basic_dic[0].InvertFilePos;

	for (int i = 0; i < basic_len; i++)
	{
		int tmp = i, tmpLen = 0, start = 0;
		double tmpMass = 0.0;

 		if (i <= max_miss_site)       //�Կ�ͷλ��Ϊ��㣬ö��
		{
			while (tmp >= 0)
			{
				tmpLen += basic_dic[tmp].Len;
				tmpMass += basic_dic[tmp].Mass;
				tmp--;
			}
			tmpMass += H2O;		//	�Ķ�����Ҫ����ˮ��������

			if (tmpLen>=min_pep_length&&tmpLen<=max_pep_length && tmpMass>=min_mass&&tmpMass<=max_mass)
			{
				_ASSERTE(_CrtCheckMemory());
				pep_dic[++pep_dic_pos].Pos = 0;
				pep_dic[pep_dic_pos].Len = tmpLen;
				pep_dic[pep_dic_pos].Mass = tmpMass;
				pep_dic[pep_dic_pos].MissCleaveSites = i;
				pep_dic[pep_dic_pos].InvertFilePos = proID;
				//�����¶�����
				pep_dic[pep_dic_pos].Godel = Cal_Godel_value(str, tmpLen);
				_ASSERTE(_CrtCheckMemory());
			}

			if (M_flag)             //�����N��M,�����г�M�����
			{
				tmpLen--;
				tmpMass -= aa_mass['M' - 'A'];		//ǰ���Ѿ�����ˮ���ӣ��˴�����Ҫ��

				if (tmpLen >= min_pep_length&&tmpLen <= max_pep_length && tmpMass >= min_mass&&tmpMass <= max_mass)
				{
					_ASSERTE(_CrtCheckMemory());
					pep_dic[++pep_dic_pos].Pos = 1;			//�г�M�����1
					pep_dic[pep_dic_pos].Len = tmpLen;
					pep_dic[pep_dic_pos].Mass = tmpMass;
					pep_dic[pep_dic_pos].MissCleaveSites = i;
					pep_dic[pep_dic_pos].InvertFilePos = proID;
					pep_dic[pep_dic_pos].Godel = Cal_Godel_value(str+1, tmpLen);	//�����¶�����
					_ASSERTE(_CrtCheckMemory());
				}
			}

		}
		
		for (int j = i + 1; j <= i + 1 + max_miss_site&&j < basic_len; j++)//��iΪ��㣬ö��
		{
			tmp = j; tmpLen = 0; tmpMass = 0.0; start = basic_dic[i + 1].Pos;
			while (tmp>i)
			{
				tmpLen += basic_dic[tmp].Len;
				tmpMass += basic_dic[tmp].Mass;
				tmp--;
			}
			tmpMass += H2O;			//�Ķ�����Ҫ����ˮ��������

			if (tmpLen >= min_pep_length&&tmpLen <= max_pep_length && tmpMass >= min_mass&&tmpMass <= max_mass)
			{
				_ASSERTE(_CrtCheckMemory());
				pep_dic[++pep_dic_pos].Pos = start;
				pep_dic[pep_dic_pos].Len = tmpLen;
				pep_dic[pep_dic_pos].Mass = tmpMass;
				pep_dic[pep_dic_pos].MissCleaveSites = j - i - 1;
				pep_dic[pep_dic_pos].InvertFilePos = proID;
				pep_dic[pep_dic_pos].Godel = Cal_Godel_value(str+start, tmpLen);	//�����¶�����
				_ASSERTE(_CrtCheckMemory());
			}
		}//end for j

	}//end for i
	/*test*/
	//printf("����ø�еĽ��Ϊ��\n");
	//printf_dic_pepStr(pep_dic, pep_dic_pos, str);
	return pep_dic_pos - start_pos+1;		//���ز������Ķ���Ŀ
}

/*
Function:		create_basic_pepDic
Description:    �����Ķλ��ʵ䣨����ø�е��ø�У�
Calls:			isBelong
Input:			ezy_terminal��ø�����ͣ�N��orC��
				proID����ǰsq���ݶ�Ӧ�ĵ�����ID
				sq_str����ø�е�sq����
				aa:�������
				basic_dic����Ż�ø�в������Ķ�
Return:         int�ͣ�basic_dic�ĳ���
*/

void  printf_a_pep(const char*sq_str,int start,int len)
{
	for (int i = start; i < start + len; i++)
		printf("%c",sq_str[i]);
	printf("\n");
}

int create_basic_pepDic(const char ezy_terminal, const int proID, const char*sq_str, const double aa[], Peptide_Dictionary* basic_dic)
{
	int sq_pos = 0;
	int i = -1;

	if (ezy_terminal == 'N')			//�����N��ø��
	{
		while (sq_str[sq_pos] != '\0')
		{
			basic_dic[++i].Pos = sq_pos;
			basic_dic[i].InvertFilePos = proID;			//���ű��λ����ʱдproID
			double ans_mass = aa[sq_str[sq_pos] - 'A'];

			while (!isBelong(sq_str[++sq_pos]) && sq_str[sq_pos] != '\0')
				ans_mass += aa[sq_str[sq_pos] - 'A'];

			basic_dic[i].Mass = ans_mass;
			basic_dic[i].Len = sq_pos - basic_dic[i].Pos;

			/*test*/
		//	printf_a_pep(sq_str, basic_dic[i].Pos, basic_dic[i].Len);
		}
	}
	else									//�����C��ø��
	{
		while (sq_str[sq_pos] != '\0')
		{
			basic_dic[++i].Pos = sq_pos;
			basic_dic[i].InvertFilePos = proID;
			double ans_mass = aa[sq_str[sq_pos] - 'A'];

			while (!isBelong(sq_str[sq_pos]) && sq_str[sq_pos] != '\0')
			{
				sq_pos++;
				ans_mass += aa[sq_str[sq_pos] - 'A'];
			}
			if (sq_str[sq_pos] != '\0')		//Ҫ�������һ��λ�õ����⴦��
			{
				basic_dic[i].Mass = ans_mass + aa[sq_str[sq_pos] - 'A'];
				sq_pos++;
				basic_dic[i].Len = sq_pos - basic_dic[i].Pos;

			}
			else
			{
				basic_dic[i].Mass = ans_mass;
				basic_dic[i].Len = sq_pos - basic_dic[i].Pos;
				break;
			}
			/*test*/
			//printf_a_pep(sq_str, basic_dic[i].Pos, basic_dic[i].Len);
		}

	}

	return i + 1;
}




/*
Function:		sort_word
Description:    ��double��ĳ���ֽ����򣬶�array����
Calls:			��
Input:			array�������������
				n�����������ݵĸ���
				word����double��word�ֽڵĴ�С����
				buf����ʱ������
Return:         void
*/
static void sort_word(Peptide_Dictionary *array, size_t n, int word,Peptide_Dictionary*buf)
{

	int *cnt = new int[65536];
	//����
	for (int i = 0; i < 65536; i++)
		cnt[i] = 0;
	//����
	for (size_t i = 0; i < n; i++)
	{
		//uint16 tmp = ((uint16 *)(array + i))[word];
		uint16 t[4];
		double*mass = (double*)t;
		*mass = (array + i)->Mass;
		uint16 tmp = t[word];
		cnt[tmp]++;
	}
	//������ʼλ��
	for (size_t i = 1; i < 65536; i++)
	{
		cnt[i] += cnt[i - 1];
	}

	memcpy(buf, array, sizeof(Peptide_Dictionary)*n);

	for (int i = (int)n - 1; i >= 0; i--)
	{
		//uint16 tmp = ((uint16 *)&(buf[i].Mass))[word];
		uint16 t[4];
		double*mass = (double*)t;
		*mass = (buf + i)->Mass;
		uint16 tmp = t[word];

		array[cnt[tmp] - 1] = buf[i];
		cnt[tmp]--;
	}

	delete[]cnt;
}

/*
Function:		fsort_mass
Description:    ��������array��������
Calls:			sort_word
Input:			array�������������
				n�����������ݵĸ���
Return:         void
*/
void fsort_mass(Peptide_Dictionary *array, size_t n)//�������������򣬲���Ҫ���Ǹ��������
{

	Peptide_Dictionary *buf = (Peptide_Dictionary *)malloc(sizeof(Peptide_Dictionary)*n);

	sort_word(array, n, 0, buf);
	sort_word(array, n, 1, buf);
	sort_word(array, n, 2, buf);
	sort_word(array, n, 3, buf);

	memcpy(buf, array, sizeof(Peptide_Dictionary)*n);

	memcpy(array, buf, sizeof(Peptide_Dictionary)*n);

	free(buf);

}


/*
Function:		adjust
Description:    �������ؽ�
Calls:			��
Input:			ls��Ԫ��λ�� 
				b�����Ԫ��
				k:�鲢��Ŀ
				s:
Return:         void
*/
static void adjust(int ls[], Peptide_Dictionary b[], int k, int s)
{
	int i, t;
	t = (s + k) >> 1;
	while (t>0)
	{
		if (b[s]>b[ls[t]])
		{
			i = s;
			s = ls[t];
			ls[t] = i;
		}
		t = t >> 1;
	}
	ls[0] = s;
}

/*
Function:		create_loser_tree
Description:    ����������
Calls:			adjust
Input:			ls��Ԫ��λ��
				b�����Ԫ��
				k:�鲢��Ŀ
Return:         void
*/
static void create_loser_tree(int ls[], Peptide_Dictionary b[], int k)
{
	b[k].Mass = (numeric_limits<double>::lowest)();
	int i;
	for (i = 0; i<k; i++)
		ls[i] = k;
	for (int i = k - 1; i >= 0; i--)
		adjust(ls, b, k, i);
}




/*
Function:		read_sq
Description:    ��sq�ļ��ж�ȡһ��sq����
Calls:			��
Input:			filenum:sq�����ļ�����Ŀ
				file_id:sq�ļ�id
				index_fp:ָ�򵰰���index�ļ���ָ��
				sq_fp:ָ��sq�����ļ���ָ��
				sq_str:��Ŷ�ȡ��sq����
Return:         bool�ͣ�����ɹ���ȡһ�����ݣ�����true;����false
*/
bool read_sq(const size_t filenum, size_t&file_id, FILE*&index_fp, FILE*&sq_fp, char*sq_str)
{
	Protein_Index_Table pro_index;
	if (file_id <= filenum)
	{
		fread(&pro_index, sizeof(Protein_Index_Table), 1, index_fp);
		if (feof(index_fp))
		{
			fclose(index_fp);
			fclose(sq_fp);
			++file_id;
			if (file_id > filenum)
				return false;
			else
			{
				char suffix[5];
				sprintf(suffix, "%d", file_id);
				index_fp = fopen(new_file_name("proIndex", suffix, output_dir_result).c_str(), "rb");
				check_fp(index_fp);

				sq_fp = fopen(new_file_name("sq", suffix, output_dir_result).c_str(), "rb");
				check_fp(index_fp);

				fread(&pro_index, sizeof(Protein_Index_Table), 1, index_fp);
			}
		}

		fread(sq_str, sizeof(char), pro_index.SQ_Len, sq_fp);
		_ASSERTE(_CrtCheckMemory());
		sq_str[pro_index.SQ_Len] = '\0';
		_ASSERTE(_CrtCheckMemory());

		return true;
	}
	return false;
}





/*
Function:		ans_enzyme
Description:    ͳ��һ��sq����ø�е����Ŀ
Calls:			��
Input:			sq_str:sq����
Return:         int�ͣ�sq��ø�е���
*/
int ans_enzyme(const char*sq_str)    
{
	int i = -1;          
	int ans = 0;
	while (sq_str[++i] != '\0')
		if (isBelong(sq_str[i])) ans++;

	return ans;


}

/*
Function:		isBelong
Description:    �ж�ĳ���ַ��Ƿ���ø��λ��
Calls:			��
Input:			a:���жϵ��ַ�
Return:         bool�ͣ�a��ø�е㣬����true;����false
*/
static bool isBelong(const char a)		//�ж�ĳ���ַ��Ƿ���ø�е�
{
	int i = -1;                       //ע�ⲻҪ�ı�enzyme��ֵ
	while (enzyme_aa[++i] != '\0')
		if (enzyme_aa[i] == a) return true;

	return false;
}


/*
Function:		printf_pepDic
Description:    ����Ķδʵ�������ݵ�����̨��Peptide_Dictionary��
Calls:			��
Input:			Dic:�����������
				dic_len��Dic�ĳ���
Return:         void
*/
void printf_pepDic(const Peptide_Dictionary *Dic, const int dic_len)
{
	for (int i = 0; i < dic_len; i++)
		printf("proID(InvertPos):%-15dpos: %-15dlen%-15dmass: %-15.5fmissSite:%-15dG: %-15.5f\n",
			Dic[i].InvertFilePos,
			Dic[i].Pos,
			Dic[i].Len,
			Dic[i].Mass,
			Dic[i].MissCleaveSites,
			Dic[i].Godel);
}

/*
Function:		printf_pepDicFinal�� Peptide_Dictionary_FINAL��
Description:    ����Ķδʵ�������ݵ�����̨
Calls:			��
Input:			Dic:�����������
				dic_len��Dic�ĳ���
Return:         void
*/
void printf_pepDicFinal(const Peptide_Dictionary_FINAL*Dic, const int dic_len)
{
	for (int i = 0; i < dic_len; i++)
		printf("proID:%-15dpos: %-15dlen%-15dmass: %-15.5fmissSite:%-15d\n",
			Dic[i].InvertFilePos,
			Dic[i].Pos,
			Dic[i].Len,
			Dic[i].Mass,
			Dic[i].MissCleaveSites);
}

/*
Function:		create_prime_array
Description:    ������������
Calls:			��
Input:			��
Return:         void
*/
static int create_prime_array()					//������������
{
	int count = 0;
	bool *temp = new bool[prime_range];
	for (int i = 0; i != prime_range; ++i)
		temp[i] = true;
	temp[2] = true;
	for (int i = 2; i != prime_range / 2; ++i)
	{
		if (temp[i])
		{
			int j = 2;
			while (i*j<prime_range)//�����ı�������������
			{
				temp[i*j] = false;
				++j;
			}
		}
	}

	count = 0;
	for (int i = 2; i != prime_range; ++i)
	{
		if (temp[i])
		{
			
			log_prime[count++] = log10((double)i);
			/*test*/
			//printf("i = %d  log_prime[%d] = %.5f\n ",i,count-1,log_prime[count-1]);
		}
	}

	delete[] temp;
	return count;        
}

/*
Function:		Cal_Godel_value
Description:    �����¶�����
Calls:			��
Input:			str��ָ������ַ�����ʼλ��
				len���ַ�������
Return:         double,str��Ӧ�ĸ�¶�ֵ
*/
static double Cal_Godel_value(const char*str, const int len)//�����¶�����
{
	double G = 0.0;
	for (int i = 0; i < len; i++)
		G += (str[i] - 'A')*log_prime[i];
	return G;
}



