#ifndef  _PRO_INDEX_H_
#define  _PRO_INDEX_H_

#include"file_operation.h"
#include"data_struct.h"


bool read_proInfo(char*&PDATA, char*AC_str, char*DE_str, char*SQ_str, int&AC_len, int &DE_len, int&SQ_len); //��ȡһ��������
int create_pro_index();		//��������������
void search_pro_index();	//������������ѯ
void test_create_pro_idex();		//���Ժ���


#endif // ! _PRO_INDEX_H_