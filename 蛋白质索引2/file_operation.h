#ifndef _FILE_OPERATION_H
#define   _FILE_OPERATION_H

#include<string>
#include<iostream>
using namespace std;



void check_fp(FILE*fp);																//检查fp异常
string new_file_name(const char*pre, const char*suff, const char*out_dir_result);   //新建文件名
void param_config(char*para, FILE*fp);												//配置一条index文件中的param数据
void Init_Config(FILE*fp_index);													//配置param文件


#endif // !_FILE_OPERATION_H