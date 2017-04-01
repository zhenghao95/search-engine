#include"protein_index.h"
#include<Windows.h>
#include<tchar.h>

/*外部全局变量*/
extern char input_path_fasta[300], output_dir_result[300];

/*全局变量*/
const int  max_ac_len = 500;
const int  max_de_len = 1000;


/*
Function:		read_proInfo
Description:    读取一条蛋白质数据
Calls:          无
Input:          PDATA：指向待读取字符的起始位置
				AC_str:存放AC数据
				DE_str:存放DE数据
				SQ_str:存放SQ数据
				AC_len:存AC数据的长度
				DE_len:存放DE数据的长度
				SQ_len:存放SQ数据的长度
Return:         bool型;成功读取一条蛋白质数据，返回true;否则，返回false
*/
bool read_proInfo(char*&PDATA, char*AC_str, char*DE_str, char*SQ_str, int&AC_len, int &DE_len, int&SQ_len)
{
	char*p = PDATA;
	if (*p == '\0')return false;

	p++;
	_ASSERTE(_CrtCheckMemory());
	while (*p != ' '&&*p!='\t'&&*p!='\n'&&*p!='\v'&&*p!='\f'&&*p!='\r')AC_str[AC_len++] = *p++;
	//AC_str[AC_len++] = '.';
	_ASSERTE(_CrtCheckMemory());

	*p++;
	_ASSERTE(_CrtCheckMemory());
	while (*p != '\n'&&*p!='\r')DE_str[DE_len++] = *p++;
	//DE_str[DE_len++] = '.';
	_ASSERTE(_CrtCheckMemory());

	_ASSERTE(_CrtCheckMemory());
	while (*p != '>'&&*p != '\0')
	{
		if (*p >= 'A'&&*p <= 'Z')
			SQ_str[SQ_len++] = *p;
		p++;

	}
	//SQ_str[SQ_len++] = '.';
	_ASSERTE(_CrtCheckMemory());
	PDATA = p;

	//DE_str[DE_len++] = '\0';
	//AC_str[AC_len++] = '\0';
	//SQ_str[SQ_len++] = '\0';
	return true;
}

/*
Function:		create_pro_index
Description:    创建蛋白质索引
Calls:          read_proInfo
Input:          无
Return:         int型;总共处理的蛋白质数目
*/
int create_pro_index()
{
	HANDLE hMapFile;      // handle for the file's memory-mapped region
	HANDLE hFile;         // the file handle
	BOOL bFlag;           // a result holder
	DWORD dwFileSize;     // temporary storage for file sizes
	LPVOID lpMapAddress;  // pointer to the base address of the
						  // memory-mapped region
	char * pData;         // pointer to the data

						  //convert char to wchar
	DWORD dwPathLen = MultiByteToWideChar(CP_ACP, 0, input_path_fasta, -1, NULL, 0);
	wchar_t*lpcTheFile = new wchar_t[dwPathLen];
	MultiByteToWideChar(CP_ACP, 0, input_path_fasta, -1, lpcTheFile, dwPathLen);

	hFile = CreateFile(lpcTheFile, GENERIC_READ, 0, NULL, OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	if (hFile == INVALID_HANDLE_VALUE)
	{
		_tprintf(TEXT("hFile is NULL\n"));
		_tprintf(TEXT("Target file is %s\n"),
			lpcTheFile);
		return 4;
	}

	dwFileSize = GetFileSize(hFile, NULL);
	_tprintf(TEXT("hFile size: %10d\n"), dwFileSize);


	hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
	if (hMapFile == NULL)
	{
		_tprintf(TEXT("hMapFile is NULL: last error: %d\n"), GetLastError());
		return (2);
	}

	lpMapAddress = MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, dwFileSize);
	if (lpMapAddress == NULL)
	{
		_tprintf(TEXT("lpMapAddress is NULL: last error: %d\n"), GetLastError());
		return 3;
	}
	pData = (char *)lpMapAddress;             //获取指向文件首字符的指针


	Protein_Meta pro_meta;			//meta datda
	pro_meta.File_Num = 0;
	pro_meta.Pro_Num = 0;
	pro_meta.offset = sizeof(Protein_Meta);
	bool first_write_to_meta = 1;

	P pro_id;
	pro_id.first = -1;
	pro_id.second = -1;

	int max_idxTable_size = max_mem_size / (sizeof(Protein_Index_Table)) + 5;
	Protein_Index_Table *pro_idx_table = new Protein_Index_Table[max_idxTable_size];		//index datda
	int pro_idx_table_len = 0;

	char*ac_str = new char[max_ac_len];							//ac,de,sq data
	char*de_str = new char[max_de_len];
	char*sq_str = new char[max_sq_len];
	int ac_len = 0, de_len = 0, sq_len = 0;


	FILE*meta_fp, *index_fp, *ac_fp, *de_fp, *sq_fp;			//file pointer
	meta_fp = fopen(new_file_name("pro", ".meta",  output_dir_result).c_str(), "wb");
	check_fp(meta_fp);

	char suffix[6];
	sprintf(suffix, "%d", ++pro_meta.File_Num);
	index_fp = fopen(new_file_name("proIndex", suffix,  output_dir_result).c_str(), "wb");
	check_fp(index_fp);

	ac_fp = fopen(new_file_name("ac", suffix,  output_dir_result).c_str(), "wb+");
	check_fp(ac_fp);
	de_fp = fopen(new_file_name("de", suffix,  output_dir_result).c_str(), "wb+");
	check_fp(de_fp);
	sq_fp = fopen(new_file_name("sq", suffix,  output_dir_result).c_str(), "wb+");
	check_fp(sq_fp);

	while (read_proInfo(pData, ac_str, de_str, sq_str, ac_len, de_len, sq_len))
	{
		pro_meta.Pro_Num++;
		if (pro_idx_table_len >= max_idxTable_size)
		{
			fwrite(pro_idx_table, sizeof(Protein_Index_Table), pro_idx_table_len, index_fp);
			pro_idx_table_len = 0;
			fclose(index_fp);
			fclose(ac_fp);
			fclose(de_fp);
			fclose(sq_fp);

			sprintf(suffix, "%d", ++pro_meta.File_Num);

			index_fp = fopen(new_file_name("proIndex", suffix, output_dir_result).c_str(), "wb");
			check_fp(index_fp);

			ac_fp = fopen(new_file_name("ac", suffix, output_dir_result).c_str(), "wb+");
			check_fp(ac_fp);
			de_fp = fopen(new_file_name("de", suffix, output_dir_result).c_str(), "wb+");
			check_fp(de_fp);
			sq_fp = fopen(new_file_name("sq", suffix, output_dir_result).c_str(), "wb+");
			check_fp(sq_fp);

			if (first_write_to_meta)
			{
				first_write_to_meta = false;
				fseek(meta_fp, sizeof(Protein_Meta), SEEK_SET);
			}
			pro_id.second = pro_meta.Pro_Num - 1;
			pro_id.first++;

			fwrite(&pro_id, sizeof(pro_id), 1, meta_fp);
			pro_id.first = pro_id.second;
		}//end if(idex_size>mem_linit_size)
		//添加一条数据到索引文件中
		pro_idx_table[pro_idx_table_len].AC_Pos = ftell(ac_fp);
		pro_idx_table[pro_idx_table_len].DE_Pos = ftell(de_fp);
		pro_idx_table[pro_idx_table_len].SQ_Pos = ftell(sq_fp);

		pro_idx_table[pro_idx_table_len].ac_Len = ac_len;
		pro_idx_table[pro_idx_table_len].DE_Len = de_len;
		pro_idx_table[pro_idx_table_len].SQ_Len = sq_len;
		pro_idx_table_len++;

		//将ac,de，sq写到对应的文件中
		fwrite(ac_str, sizeof(char), ac_len, ac_fp);
		fwrite(de_str, sizeof(char), de_len, de_fp);
		fwrite(sq_str, sizeof(char), sq_len, sq_fp);
		//printf("蛋白质ID： %d   ac：%s\n", pro_meta.Pro_Num - 1, ac_str);
		ac_len = 0;
		de_len = 0;
		sq_len = 0;
		
	}//end while(read a protein data)


	if (pro_idx_table_len!=0)
	{
		fwrite(pro_idx_table, sizeof(Protein_Index_Table), pro_idx_table_len, index_fp);
		pro_idx_table_len = 0;

		if (first_write_to_meta)
		{
			first_write_to_meta = false;
			fseek(meta_fp, sizeof(Protein_Meta), SEEK_SET);
		}

		pro_id.second = pro_meta.Pro_Num - 1;
		pro_id.first++;
		fwrite(&pro_id, sizeof(P), 1, meta_fp);
	}

	fclose(index_fp);
	fclose(ac_fp);
	fclose(de_fp);
	fclose(sq_fp);


	fseek(meta_fp, 0, SEEK_SET);
	fwrite(&pro_meta, sizeof(Protein_Meta), 1, meta_fp);
	fclose(meta_fp);


	// Close the file mapping object and the open file
	bFlag = UnmapViewOfFile(lpMapAddress);
	bFlag = CloseHandle(hMapFile);

	if (!bFlag)
	{
		_tprintf(TEXT("\nError %ld occurred closing the mapping object!"),
			GetLastError());
	}

	bFlag = CloseHandle(hFile);

	if (!bFlag)
	{
		_tprintf(TEXT("\nError %ld occurred closing the file!"),
			GetLastError());
	}

	delete[]lpcTheFile;
	delete[]pro_idx_table;
	delete[]ac_str;
	delete[]de_str;
	delete[]sq_str;

	return pro_meta.Pro_Num;
}

/*
Function:		search_pro_index
Description:    蛋白质索引查询
Calls:          
Input:         
Return:         bool型;成功读取一条蛋白质数据，返回true;否则，返回false
*/
void search_pro_index()
{

}

void test_create_pro_idex()
{
	FILE*meta_fp, *idx_fp, *ac_fp, *de_fp, *sq_fp;

	Protein_Meta pro_meta;
	P pro_id;
	size_t file_id = 0;

	Protein_Index_Table pro_idx_table;

	char *ac_str = new char[max_ac_len];
	char*de_str = new char[max_de_len];
	char*sq_str = new char[max_sq_len];

	meta_fp = fopen(new_file_name("pro", ".meta",  output_dir_result).c_str(), "rb");
	check_fp(meta_fp);

	fread(&pro_meta, sizeof(Protein_Meta), 1, meta_fp);
	printf("total pronum:%d\n", pro_meta.Pro_Num);
	printf("total filenum:%d\n", pro_meta.File_Num);

	/*while (file_id <= pro_meta.File_Num)
	{
	fread(&pro_id, sizeof(P), 1, meta_fp);
	char suffix[10];
	sprintf(suffix, "%d", file_id);
	idx_fp = fopen(new_file_name("proIndex", suffix,  output_dir_result).c_str(), "rb");
	ac_fp = fopen(new_file_name("ac", suffix,  output_dir_result).c_str(), "rb");
	de_fp = fopen(new_file_name("de", suffix,  output_dir_result).c_str(), "rb");
	sq_fp = fopen(new_file_name("sq", suffix,  output_dir_result).c_str(), "rb");

	for (int i = pro_id.first; i <= pro_id.second; i++)
	{
	fread(&pro_idx_table, sizeof(Protein_Index_Table), 1, idx_fp);
	fread(ac_str, sizeof(char), pro_idx_table.ac_Len, ac_fp);
	fread(de_str, sizeof(char), pro_idx_table.DE_Len, de_fp);
	fread(sq_str, sizeof(char), pro_idx_table.SQ_Len, sq_fp);

	ac_str[pro_idx_table.ac_Len] = '\0';
	de_str[pro_idx_table.DE_Len] = '\0';
	sq_str[pro_idx_table.SQ_Len] = '\0';

	printf("pro.%d\n", i);
	printf("ac:%s\n", ac_str);
	printf("de:%s\n", de_str);
	printf("sq:%s\n\n", sq_str);

	}
	fclose(idx_fp);
	fclose(ac_fp);
	fclose(de_fp);
	fclose(sq_fp);
	file_id++;
	}*/

	fclose(meta_fp);
	delete[]ac_str;
	delete[]de_str;
	delete[]sq_str;

}
