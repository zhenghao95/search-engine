#ifndef  _PRO_INDEX_H_
#define  _PRO_INDEX_H_

#include"file_operation.h"
#include"data_struct.h"


bool read_proInfo(char*&PDATA, char*AC_str, char*DE_str, char*SQ_str, int&AC_len, int &DE_len, int&SQ_len); //读取一条蛋白质
int create_pro_index();		//创建蛋白质索引
void search_pro_index();	//蛋白质索引查询
void test_create_pro_idex();		//测试函数


#endif // ! _PRO_INDEX_H_