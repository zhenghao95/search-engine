#include<iostream>
#include<cstdio>
#include"protein_index.h"
#include"peptide_index.h"
using namespace std;



int main()
{
	FILE*fp_index = fopen("index.param", "r");
	check_fp(fp_index);
	
	Init_Config(fp_index);

	fclose(fp_index);

	int proNum = create_pro_index();
	printf("蛋白质的总数为：%d\n", proNum);

	//test_create_pro_idex(output_dir_result);

	create_peptide_index();
	//int pepNum = Digestion();
	//printf("肽段数目为：%d\n", pepNum);
	//for (int i = 0; i < 26; i++)
	//{
	//	printf("%.5f\n", a[i]);
	//}


	return 0; 
}

