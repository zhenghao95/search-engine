#include"peptide_index.h"
#include<cmath>
#include <cstdlib>
#include <cstring>
#include <cstdbool>

typedef unsigned long long uint64;
typedef long long int64;
typedef unsigned short uint16;

/*外部全局变量-from file_operation.c */
extern char output_dir_result[300];
extern char enzyme_aa[5];
extern char enzyme_terminal;
extern double min_mass, max_mass;
extern int min_pep_length, max_pep_length, max_miss_site;



/*全局变量*/
const int max_ezy_num = 800;		//最大酶切位点数
#define prime_num 120				//素数的数目在120个以内
#define prime_range 600				//产生600以内的素数
double log_prime[prime_num];			
const int max_pepDic_size = max_mem_size / sizeof(Peptide_Dictionary);
const int max_pepDicFINALsize = max_mem_size / sizeof(Peptide_Dictionary_FINAL);  //每个肽段词典文件中肽段的条数
const int max_same_pep = 200;		//	有相同质量的肽段最多两百条

double aa_mass[26] = { 71.03711,166.99836,103.00919,		//氨基酸质量表
115.02694,129.04259,147.06841,57.02146,
137.05891,113.08406,181.01401,128.09496,
113.08406,131.04048,114.04293,243.02965,
97.05276,128.05858,156.10111,87.03203 ,
101.04768 ,0.00000,99.06841,186.07931,
113.08406,163.06332,128.55059 };



/*
Function:		create_peptide_index
Description:    创建肽段索引
Calls:			create_prime_array
				Digestion
				k_merge
Input:			无
Return:         肽段词典文件数
*/
int create_peptide_index()
{
	create_prime_array();

	int tmpFileNum = Digestion();	
	printf("临时文件数为：%d\n", tmpFileNum);

	int *ls = new int[tmpFileNum+3];
	Peptide_Dictionary*b = new Peptide_Dictionary[tmpFileNum + 3];
	
	int dic_num = k_merge(ls, b, tmpFileNum);		
	//printf("dic文件数为：%d\n", dic_num);

	delete[]ls;
	delete[]b;

	//return dic_num;
	return 0;
}



/*
Function:		Digestion
Description:    对数据库中的所有蛋白质酶切
Calls:			read_sq
				create_basic_pepDic
				specific_enzyme
				fsort_mass
Input:			无
Return:         生成临时文件pepTmp的数目
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

	int tmpFileID = -1;			//临时文件指针
	char suff[6];
	FILE*fp_tmp;

	/*test*/
	int total_pep_num = 0;	//肽段总数

	while (read_sq(pro_meta.File_Num, file_id, index_fp, sq_fp, sq_str))
	{
		
		proID++;
		//printf("proID:%d	%s\n",proID, sq_str);
		int basic_len = create_basic_pepDic(enzyme_terminal, proID, sq_str, aa_mass, basic_dic);//产生基词典

		if (container_pos + (1 + basic_len)*(max_miss_site + 1) + 1 > max_pepDic_size)//如果此次酶切有溢出的风险，排序写到临时文件
		{
		//	printf("排序前的结果为：\n");
		//	printf_pepDic(Container, container_pos + 1);
			//排序
			fsort_mass(Container, container_pos + 1);
		//	printf("排序后的结果为：\n");
			//printf_pepDic(Container, container_pos + 1);
			

			//写到临时文件
			sprintf(suff, "%d", ++tmpFileID);
			fp_tmp = fopen(new_file_name("pepTmp", suff, output_dir_result).c_str(), "wb");//打开一个临时文件
			fwrite(Container, sizeof(Peptide_Dictionary), container_pos + 1, fp_tmp);
			
			fclose(fp_tmp);
			//清空容器
			container_pos = -1;
		}

		bool mFlag = 0;//判断是否需要考虑M端
		if (enzyme_terminal == 'N'&&sq_str[0] == 'M')mFlag = 1;
		int pepNum = specific_enzyme(sq_str, mFlag,  Container, container_pos, basic_dic, basic_len);	//特异酶切
		/*test*/
		total_pep_num += pepNum;
		//printf("proID: %d cut pep num is ：%d \n", proID,pepNum);
	}

	/*test*/
	printf("after digestion, total pep num is:%d\n", total_pep_num);

	//处理剩余的肽段
	if (container_pos != -1)
	{
		//printf("排序前的结果为：\n");
		//printf_pepDic(Container, container_pos + 1);
	
		fsort_mass(Container, container_pos + 1);
	//	printf("排序后的结果为：\n");
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

	return tmpFileID + 1;          //返回临时文件总数
}

/*
Function:		k_merge
Description:    对酶切肽段归并去冗余
Calls:			create_loser_tree
Input:			ls：元素位置
b：存放元素
k:归并数目
Return:         int型，肽段词典文件数
*/
int k_merge(int ls[], Peptide_Dictionary b[], int k)
{

	FILE* *fp_in = new FILE*[k];		//临时文件
	char suff[5];
	for (int i = 0; i < k; i++)
	{
		sprintf(suff, "%d", i);
		fp_in[i] = fopen(new_file_name("pepTmp", suff, output_dir_result).c_str(), "rb");
		check_fp(fp_in[i]);
	}

	int dic_ID = -1;
	sprintf(suff, "%d", ++dic_ID);
	FILE*dic_fp = fopen(new_file_name("pepDictionary", suff, output_dir_result).c_str(), "wb");//肽段词典文件
	check_fp(dic_fp);

	FILE* Ivt_fp = fopen(new_file_name("index", "Invert", output_dir_result).c_str(), "wb");	//倒排索引文件
	check_fp(Ivt_fp);

	Peptide_Dictionary_FINAL  *dic = new Peptide_Dictionary_FINAL[max_pepDicFINALsize];
	int dic_pos = -1;

	for (int i = 0; i < k; ++i)
		fread(&b[i], sizeof(Peptide_Dictionary), 1, fp_in[i]);

	create_loser_tree(ls, b, k);		//创建败者树
	int q = ls[0];
	Peptide_Dictionary f = b[q];		//得到堆顶的元素，存入f

	int*	pro_id = new int[max_same_pep];				    //存放蛋白质ID的数组
	int pro_id_pos = -1;
	pro_id[++pro_id_pos] = f.InvertFilePos;	//存入f的proID

	fread(&b[q], sizeof(Peptide_Dictionary), 1, fp_in[q]);		//读入一条数据
	create_loser_tree(ls, b, k);		//创建败者树

	/*test*/
	int pep_num_total = 0;

	while (b[ls[0]].Mass != (numeric_limits<double>::max)())
	{
		q = ls[0];
		Peptide_Dictionary l = b[q];

		if (l.Mass == f.Mass && l.Godel == f.Godel)			//如果是同一条肽段
		{
			if (l.InvertFilePos != f.InvertFilePos)
			{
				_ASSERTE(_CrtCheckMemory());
				pro_id[++pro_id_pos] = l.InvertFilePos;
				_ASSERTE(_CrtCheckMemory());
			}
		}
		else                          //如果说l,f不是同一条肽段
		{
			//得到倒排表的位置
			fpos_t pos;
			fgetpos(Ivt_fp, &pos);

			//将pro_id_pos,pro_id数据写入倒排表
			int proNum = pro_id_pos + 1;
			fwrite(&proNum, sizeof(int), 1, Ivt_fp);
			fwrite(pro_id, sizeof(int), pro_id_pos + 1, Ivt_fp);
			pro_id_pos = -1;

			//如果dic写满了，将dic写到文件中
			if (dic_pos == max_pepDicFINALsize - 1)
			{
				/*test*/
				pep_num_total += dic_pos + 1;

				//输出整体排序去冗余后的数据
				//printf("归并排序去冗余后的数据为：\n");
				//printf_pepDicFinal(dic, max_pepDicFINALsize);

				fwrite(dic, sizeof(Peptide_Dictionary_FINAL), max_pepDicFINALsize, dic_fp);
				dic_pos = -1;
				fclose(dic_fp);

				sprintf(suff, "%d", ++dic_ID);
				dic_fp = fopen(new_file_name("pepDictionary", suff, output_dir_result).c_str(), "wb");//肽段词典文件
				check_fp(dic_fp);
			}
			//往dic中添加一条数据
			_ASSERTE(_CrtCheckMemory());
			dic[++dic_pos].InvertFilePos = pos;
			dic[dic_pos].Mass = f.Mass;
			dic[dic_pos].Len = f.Len;
			dic[dic_pos].MissCleaveSites = f.MissCleaveSites;
			dic[dic_pos].Pos = f.Pos;
			_ASSERTE(_CrtCheckMemory());

			//将f更新为l
			f = l;
			pro_id[++pro_id_pos] = f.InvertFilePos;
		}

		if (!fread(&b[q], sizeof(Peptide_Dictionary), 1, fp_in[q]))		//读入一条数据
			b[q].Mass = (numeric_limits<double>::max)();

		adjust(ls, b, k, q);				//调整败者树
	}

	//处理最后一组数据
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
		//输出整体排序去冗余后的数据
		//printf("归并排序去冗余后的数据为：\n");
		//printf_pepDicFinal(dic, max_pepDicFINALsize);

		fwrite(dic, sizeof(Peptide_Dictionary_FINAL), max_pepDicFINALsize, dic_fp);
		dic_pos = -1;
		fclose(dic_fp);

		sprintf(suff, "%d", ++dic_ID);
		dic_fp = fopen(new_file_name("pepDictionary", suff, output_dir_result).c_str(), "wb");//肽段词典文件
		check_fp(dic_fp);
	}
	_ASSERTE(_CrtCheckMemory());
	dic[++dic_pos].InvertFilePos = pos;
	dic[dic_pos].Mass = f.Mass;
	dic[dic_pos].Len = f.Len;
	dic[dic_pos].MissCleaveSites = f.MissCleaveSites;
	dic[dic_pos].Pos = f.Pos;
	_ASSERTE(_CrtCheckMemory());


	//将剩余的数据写到肽段文件,并关闭肽段文件
	if (dic_pos > -1)
	{
		/*test*/
		pep_num_total += dic_pos + 1;

		fwrite(dic, sizeof(Peptide_Dictionary_FINAL), dic_pos + 1, dic_fp);
		//printf("归并后的数据为：\n");
		//printf_pepDicFinal(dic, dic_pos + 1);
		dic_pos = -1;
	}

	/*test*/
	printf("去冗余后的肽段数：%d\n", pep_num_total);

	fclose(dic_fp);


	//关闭临时文件
	for (int i = 0; i < k; i++)
		fclose(fp_in[i]);
	//删除临时文件
	for (int i = 0; i < k; i++)
	{
		sprintf(suff, "%d", i);
		if (remove(new_file_name("pepTmp", suff, output_dir_result).c_str()) != 0)
			perror("Error deleting file!");
	}
	//关闭倒排表文件
	fclose(Ivt_fp);

	//释放申请空间
	delete[]fp_in;
	delete[]pro_id;
	return dic_ID + 1;		//返回肽段词典文件数目
}
/*
Function:		specific_enzyme
Description:    对一条sq做特异酶切
Calls:			无
Input:			str：指向要酶切的sq数据
				M_flag:是否存在需要考虑N端M的情况
				pep_dic：存放酶切结果的容器
				pep_dic_pos：指示pep_dic的位置
				basic_dic：基础酶切产生的肽段基词典
				basic_len：指示basic_dic的长度
Return:        int型，此次酶切总共产生的肽段数目
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

 		if (i <= max_miss_site)       //以开头位置为起点，枚举
		{
			while (tmp >= 0)
			{
				tmpLen += basic_dic[tmp].Len;
				tmpMass += basic_dic[tmp].Mass;
				tmp--;
			}
			tmpMass += H2O;		//	肽段质量要加上水分子质量

			if (tmpLen>=min_pep_length&&tmpLen<=max_pep_length && tmpMass>=min_mass&&tmpMass<=max_mass)
			{
				_ASSERTE(_CrtCheckMemory());
				pep_dic[++pep_dic_pos].Pos = 0;
				pep_dic[pep_dic_pos].Len = tmpLen;
				pep_dic[pep_dic_pos].Mass = tmpMass;
				pep_dic[pep_dic_pos].MissCleaveSites = i;
				pep_dic[pep_dic_pos].InvertFilePos = proID;
				//计算哥德尔质量
				pep_dic[pep_dic_pos].Godel = Cal_Godel_value(str, tmpLen);
				_ASSERTE(_CrtCheckMemory());
			}

			if (M_flag)             //如果是N端M,考虑切除M的情况
			{
				tmpLen--;
				tmpMass -= aa_mass['M' - 'A'];		//前面已经加了水分子，此处不需要加

				if (tmpLen >= min_pep_length&&tmpLen <= max_pep_length && tmpMass >= min_mass&&tmpMass <= max_mass)
				{
					_ASSERTE(_CrtCheckMemory());
					pep_dic[++pep_dic_pos].Pos = 1;			//切除M起点是1
					pep_dic[pep_dic_pos].Len = tmpLen;
					pep_dic[pep_dic_pos].Mass = tmpMass;
					pep_dic[pep_dic_pos].MissCleaveSites = i;
					pep_dic[pep_dic_pos].InvertFilePos = proID;
					pep_dic[pep_dic_pos].Godel = Cal_Godel_value(str+1, tmpLen);	//计算哥德尔质量
					_ASSERTE(_CrtCheckMemory());
				}
			}

		}
		
		for (int j = i + 1; j <= i + 1 + max_miss_site&&j < basic_len; j++)//以i为起点，枚举
		{
			tmp = j; tmpLen = 0; tmpMass = 0.0; start = basic_dic[i + 1].Pos;
			while (tmp>i)
			{
				tmpLen += basic_dic[tmp].Len;
				tmpMass += basic_dic[tmp].Mass;
				tmp--;
			}
			tmpMass += H2O;			//肽段质量要加上水分子质量

			if (tmpLen >= min_pep_length&&tmpLen <= max_pep_length && tmpMass >= min_mass&&tmpMass <= max_mass)
			{
				_ASSERTE(_CrtCheckMemory());
				pep_dic[++pep_dic_pos].Pos = start;
				pep_dic[pep_dic_pos].Len = tmpLen;
				pep_dic[pep_dic_pos].Mass = tmpMass;
				pep_dic[pep_dic_pos].MissCleaveSites = j - i - 1;
				pep_dic[pep_dic_pos].InvertFilePos = proID;
				pep_dic[pep_dic_pos].Godel = Cal_Godel_value(str+start, tmpLen);	//计算哥德尔质量
				_ASSERTE(_CrtCheckMemory());
			}
		}//end for j

	}//end for i
	/*test*/
	//printf("特意酶切的结果为：\n");
	//printf_dic_pepStr(pep_dic, pep_dic_pos, str);
	return pep_dic_pos - start_pos+1;		//返回产生的肽段数目
}

/*
Function:		create_basic_pepDic
Description:    创建肽段基词典（遇到酶切点就酶切）
Calls:			isBelong
Input:			ezy_terminal：酶切类型，N端orC端
				proID：当前sq数据对应的蛋白质ID
				sq_str：被酶切的sq数据
				aa:氨基酸表
				basic_dic：存放基酶切产生的肽段
Return:         int型，basic_dic的长度
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

	if (ezy_terminal == 'N')			//如果是N端酶切
	{
		while (sq_str[sq_pos] != '\0')
		{
			basic_dic[++i].Pos = sq_pos;
			basic_dic[i].InvertFilePos = proID;			//倒排表的位置暂时写proID
			double ans_mass = aa[sq_str[sq_pos] - 'A'];

			while (!isBelong(sq_str[++sq_pos]) && sq_str[sq_pos] != '\0')
				ans_mass += aa[sq_str[sq_pos] - 'A'];

			basic_dic[i].Mass = ans_mass;
			basic_dic[i].Len = sq_pos - basic_dic[i].Pos;

			/*test*/
		//	printf_a_pep(sq_str, basic_dic[i].Pos, basic_dic[i].Len);
		}
	}
	else									//如果是C端酶切
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
			if (sq_str[sq_pos] != '\0')		//要考虑最后一个位置的特殊处理
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
Description:    按double的某个字节排序，对array排序
Calls:			无
Input:			array：待排序的数组
				n：待排序数据的个数
				word：按double第word字节的大小排序
				buf：临时缓冲区
Return:         void
*/
static void sort_word(Peptide_Dictionary *array, size_t n, int word,Peptide_Dictionary*buf)
{

	int *cnt = new int[65536];
	//归零
	for (int i = 0; i < 65536; i++)
		cnt[i] = 0;
	//计数
	for (size_t i = 0; i < n; i++)
	{
		//uint16 tmp = ((uint16 *)(array + i))[word];
		uint16 t[4];
		double*mass = (double*)t;
		*mass = (array + i)->Mass;
		uint16 tmp = t[word];
		cnt[tmp]++;
	}
	//计算起始位置
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
Description:    按质量对array基数排序
Calls:			sort_word
Input:			array：待排序的数组
				n：待排序数据的个数
Return:         void
*/
void fsort_mass(Peptide_Dictionary *array, size_t n)//对质量进行排序，不需要考虑负数的情况
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
Description:    败者树重建
Calls:			无
Input:			ls：元素位置 
				b：存放元素
				k:归并数目
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
Description:    创建败者树
Calls:			adjust
Input:			ls：元素位置
				b：存放元素
				k:归并数目
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
Description:    从sq文件中读取一条sq数据
Calls:			无
Input:			filenum:sq数据文件的数目
				file_id:sq文件id
				index_fp:指向蛋白质index文件的指针
				sq_fp:指向sq数据文件的指针
				sq_str:存放读取的sq数据
Return:         bool型，如果成功读取一条数据，返回true;否则，false
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
Description:    统计一条sq数据酶切点的数目
Calls:			无
Input:			sq_str:sq数据
Return:         int型，sq总酶切点数
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
Description:    判断某个字符是否是酶切位点
Calls:			无
Input:			a:待判断的字符
Return:         bool型，a是酶切点，返回true;否则，false
*/
static bool isBelong(const char a)		//判断某个字符是否是酶切点
{
	int i = -1;                       //注意不要改变enzyme的值
	while (enzyme_aa[++i] != '\0')
		if (enzyme_aa[i] == a) return true;

	return false;
}


/*
Function:		printf_pepDic
Description:    输出肽段词典里的数据倒控制台（Peptide_Dictionary）
Calls:			无
Input:			Dic:待输出的数组
				dic_len：Dic的长度
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
Function:		printf_pepDicFinal（ Peptide_Dictionary_FINAL）
Description:    输出肽段词典里的数据倒控制台
Calls:			无
Input:			Dic:待输出的数组
				dic_len：Dic的长度
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
Description:    创建素数数组
Calls:			无
Input:			无
Return:         void
*/
static int create_prime_array()					//创建素数数组
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
			while (i*j<prime_range)//素数的倍数都不是素数
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
Description:    计算哥德尔编码
Calls:			无
Input:			str：指向计算字符串起始位置
				len：字符串长度
Return:         double,str对应的哥德尔值
*/
static double Cal_Godel_value(const char*str, const int len)//计算哥德尔编码
{
	double G = 0.0;
	for (int i = 0; i < len; i++)
		G += (str[i] - 'A')*log_prime[i];
	return G;
}



