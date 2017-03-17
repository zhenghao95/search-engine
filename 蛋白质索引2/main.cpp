#include<iostream>
#include<cstdio>
#include"protein_index.h"
//#include"file_operation.h"
using namespace std;


//index�ļ�����
int enzyme_type,
	max_miss_site,
	min_pep_length, max_pep_length,
	min_mass, max_mass;

char input_path_fasta[300], output_dir_result[300],
	 enzyme_terminal,
	 enzyme_aa[5];

void Init_Config(FILE*fp_index);


int main()
{


	FILE*fp_index = fopen("index.param", "r");
	check_fp(fp_index);
	
	Init_Config(fp_index);

	fclose(fp_index);

	create_pro_index(input_path_fasta, output_dir_result);
	

	return 0;
}

void Init_Config(FILE*fp_index)
{
	cout << "configuring..." << endl;
	char enzyme_type_str[3], enzyme_terminal_str[3], max_miss_site_str[5],
		min_pep_length_str[5], max_pep_length_str[5],
		min_mass_str[10], max_mass_str[10];

	//input output
	configure(input_path_fasta, fp_index);
	configure(output_dir_result, fp_index);

	//ø������
	configure(enzyme_type_str, fp_index);
	enzyme_type = atoi(enzyme_type_str);

	//ø�ж�
	configure(enzyme_terminal_str, fp_index);
	enzyme_terminal = *enzyme_terminal_str;

	//ø��λ��
	configure(enzyme_aa, fp_index);

	//�����©λ��
	configure(max_miss_site_str, fp_index);
	max_miss_site = atoi(max_miss_site_str);

	//min max�Ķγ���
	configure(min_pep_length_str, fp_index);
	configure(max_pep_length_str, fp_index);

	min_pep_length = atoi(min_pep_length_str);
	max_pep_length = atoi(max_pep_length_str);

	//min max �Ķ�����
	configure(min_mass_str, fp_index);
	configure(max_mass_str, fp_index);
	min_mass = atoi(min_mass_str);
	max_mass = atoi(max_mass_str);

}