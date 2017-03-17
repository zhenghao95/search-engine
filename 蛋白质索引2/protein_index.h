#include"data_struct.h"
#include"file_operation.h"
#include<Windows.h>
#include<tchar.h>

const int  max_mem_size = 1024 * 1024 * 256;
const int max_idxTable_size = max_mem_size / (sizeof(size_t) * 3) + 5;
const int max_ac_len = 50;
const int max_de_len = 500;
const int max_sq_len = 8000;

bool read_proInfo(char*&PDATA,char*AC_str,char*DE_str,char*SQ_str,int&AC_len,int &DE_len,int&SQ_len)
{
	char*p = PDATA;
	if (*p == '\0')return false;
	
	*p++;
	while (*p != ' ')AC_str[AC_len++] = *p++;
	//AC_str[AC_len++] = '.';
	
	*p++;
	while (*p != '\r')DE_str[DE_len++] = *p++;
	//DE_str[DE_len++] = '.';

	while (*p != '>'&&*p != '\0')
	{
		if (*p >='A'&&*p<='Z')
			SQ_str[SQ_len++] = *p;
		p++;
	}
	//SQ_str[SQ_len++] = '.';
	PDATA = p;

	//DE_str[DE_len++] = '\0';
	//AC_str[AC_len++] = '\0';
	//SQ_str[SQ_len++] = '\0';
	return true;
}

int create_pro_index(const char*input_path,const char*output_dir)
{
	HANDLE hMapFile;      // handle for the file's memory-mapped region
	HANDLE hFile;         // the file handle
	BOOL bFlag;           // a result holder
	DWORD dwFileSize;     // temporary storage for file sizes
	LPVOID lpMapAddress;  // pointer to the base address of the
						  // memory-mapped region
	char * pData;         // pointer to the data

	//convert char to wchar
	DWORD dwPathLen = MultiByteToWideChar(CP_ACP,0,input_path, -1, NULL, 0);
	wchar_t*lpcTheFile = new wchar_t[dwPathLen];
	MultiByteToWideChar(CP_ACP, 0, input_path, -1, lpcTheFile, dwPathLen);

	hFile = CreateFile(lpcTheFile,GENERIC_READ,0,NULL,OPEN_ALWAYS,FILE_ATTRIBUTE_NORMAL,NULL);
	if (hFile == INVALID_HANDLE_VALUE)
	{
		_tprintf(TEXT("hFile is NULL\n"));
		_tprintf(TEXT("Target file is %s\n"),
			lpcTheFile);
		return 4;
	}

	dwFileSize = GetFileSize(hFile, NULL);
	_tprintf(TEXT("hFile size: %10d\n"), dwFileSize);


	hMapFile = CreateFileMapping(hFile,NULL,PAGE_READONLY,0,0,NULL);
	if (hMapFile == NULL)
	{
		_tprintf(TEXT("hMapFile is NULL: last error: %d\n"), GetLastError());
		return (2);
	}
	
	lpMapAddress = MapViewOfFile(hMapFile,FILE_MAP_READ,0,0,dwFileSize);
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
	bool first_write = 1;
	bool file_end = 0;
	
	P pro_id;		
	pro_id.first = -1;
	pro_id.second = -1;

	Protein_Index_Table *pro_idx_table = new Protein_Index_Table[max_idxTable_size];		//index datda
	int pro_idx_table_len = 0;

	char*ac_str = new char[max_ac_len];							//ac,de,sq data
	char*de_str = new char[max_de_len];
	char*sq_str = new char[max_sq_len];
	int ac_len = 0, de_len = 0, sq_len = 0;
	//int ac_file_size = 0, de_file_size = 0, sq_file_size = 0;

	FILE*meta_fp, *index_fp, *ac_fp, *de_fp, *sq_fp;			//file pointer
	meta_fp = fopen(new_file_name("pro", ".meta", output_dir).c_str(), "wb");
	check_fp(meta_fp);

	char suffix[6];
	sprintf(suffix, "%d", pro_meta.File_Num);
	index_fp = fopen(new_file_name("proIndex", suffix, output_dir).c_str(), "wb");
	check_fp(index_fp);

	ac_fp = fopen(new_file_name("ac", suffix, output_dir).c_str(), "wb+");
	check_fp(ac_fp);
	de_fp = fopen(new_file_name("de", suffix, output_dir).c_str(), "wb+");
	check_fp(de_fp);
	sq_fp = fopen(new_file_name("sq", suffix, output_dir).c_str(), "wb+");
	check_fp(sq_fp);

	int mem_limit_size = max_mem_size / sizeof(Protein_Index_Table);   //能容index结构体的个数
	while (read_proInfo(pData, ac_str, de_str, sq_str, ac_len, de_len, sq_len))
	{
		pro_meta.Pro_Num++;
		if (pro_idx_table_len > mem_limit_size)
		{
			fclose(index_fp);
			fclose(ac_fp);
			fclose(de_fp);
			fclose(sq_fp);

			pro_id.second = pro_meta.Pro_Num - 1;
			if (first_write)
			{
				first_write = false;
				fseek(meta_fp, sizeof(Protein_Meta), SEEK_SET);
			}
			pro_id.first++;
			fwrite(&pro_id, sizeof(pro_id), 1, meta_fp);
			pro_id.first = pro_id.second;

			if (*pData == '\0')
			{
				file_end = 1;
				break;
			}

			pro_meta.File_Num++;
			sprintf(suffix, "%d",pro_meta.File_Num);

			index_fp = fopen(new_file_name("proIndex", suffix, output_dir).c_str(), "wb");
			check_fp(index_fp);

			ac_fp = fopen(new_file_name("ac", suffix, output_dir).c_str(), "wb+");
			check_fp(ac_fp);
			de_fp = fopen(new_file_name("de", suffix, output_dir).c_str(), "wb+");
			check_fp(de_fp);
			sq_fp = fopen(new_file_name("sq", suffix, output_dir).c_str(), "wb+");
			check_fp(sq_fp);
		}//end if(idex_size>mem_linit_size)
		else
		{
			pro_idx_table[pro_idx_table_len].AC_Pos = ftell(ac_fp);
			pro_idx_table[pro_idx_table_len].DE_Pos = ftell(de_fp);
			pro_idx_table[pro_idx_table_len].SQ_Pos = ftell(sq_fp);

			fwrite(pro_idx_table + pro_idx_table_len, sizeof(Protein_Index_Table), 1, index_fp);
			pro_idx_table_len++;

			fwrite(ac_str, sizeof(char), ac_len, ac_fp);
			fwrite(de_str, sizeof(char), de_len, de_fp);
			fwrite(sq_str, sizeof(char), sq_len, sq_fp);

			ac_len = 0;
			de_len = 0;
			sq_len = 0;
		}
	}//end while(read a protein data)

	
	if (!file_end)
	{
		fclose(index_fp);
		fclose(ac_fp);
		fclose(de_fp);
		fclose(sq_fp);

		if (first_write)
		{
			first_write = false;
			fseek(meta_fp, sizeof(Protein_Meta), SEEK_SET);
		}

		pro_id.second = pro_meta.Pro_Num - 1;
		pro_id.first++;
		fwrite(&pro_id, sizeof(P), 1, meta_fp);
	}


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

	return 0;
}

void search_pro_index()
{

}