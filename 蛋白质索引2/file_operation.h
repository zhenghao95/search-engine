#ifndef _FILE_OPERATION_H
#define   _FILE_OPERATION_H

#include<string>
#include<iostream>
using namespace std;



void check_fp(FILE*fp);																//���fp�쳣
string new_file_name(const char*pre, const char*suff, const char*out_dir_result);   //�½��ļ���
void param_config(char*para, FILE*fp);												//����һ��index�ļ��е�param����
void Init_Config(FILE*fp_index);													//����param�ļ�


#endif // !_FILE_OPERATION_H