#include"file_operation.h"

/*ȫ�ֱ���*/
int enzyme_type,
max_miss_site,
min_pep_length, max_pep_length;
double min_mass, max_mass;
char input_path_fasta[300], output_dir_result[300],
enzyme_terminal,
enzyme_aa[5];

/*
Function:		check_fp
Description:    ����ļ����쳣
Calls:			��
Input:			��
Return:         void
*/
void check_fp(FILE*fp)
{
	if (fp == NULL)
		perror("Error opening file!");

}

/*
Function:		new_file_name
Description:    �½�һ���ļ���
Calls:			��
Input:			pre���ļ���ǰ׺
				suff���ļ�����׺
				out_dir_result�����·��
Return:         string,�½����ļ���
*/
string new_file_name(const char*pre, const char*suff, const char*out_dir_result)
{
	char filename[300];
	strcpy(filename, out_dir_result);
	strcat(filename, pre);
	strcat(filename, suff);
	return filename;
}


/*
Function:		param_config
Description:    ����һ������
Calls:			��
Input:			para����Ų�����Ϣ
				fp:���ò������ڵ��ļ�ָ��
Return:         void
*/
void param_config(char*para, FILE*fp)
{
	char p[300];
	fgets(p, 300, fp);

	char*tmp = p;
	while (*tmp != '=')tmp++;
	tmp++;

	char*str = para;
	while (*tmp != '\n'&&*tmp != '\r'&&*tmp != ' ')*str++ = *tmp++;
	*str = '\0';
}

/*
Function:		Init_Config
Description:    �����ļ�����
Calls:			param_config
Input:			fp_index:ָ�������ļ���ָ��
Return:         void
*/
void Init_Config(FILE*fp_index)
{
	cout << "configuring..." << endl;
	char enzyme_type_str[3], enzyme_terminal_str[3], max_miss_site_str[5],
		min_pep_length_str[5], max_pep_length_str[5],
		min_mass_str[10], max_mass_str[10];

	param_config(input_path_fasta, fp_index);
	param_config(output_dir_result, fp_index);

	param_config(enzyme_type_str, fp_index);
	enzyme_type = atoi(enzyme_type_str);

	//ø�ж�
	param_config(enzyme_terminal_str, fp_index);
	enzyme_terminal = *enzyme_terminal_str;

	//ø��λ��
	param_config(enzyme_aa, fp_index);

	//�����©λ��
	param_config(max_miss_site_str, fp_index);
	max_miss_site = atoi(max_miss_site_str);

	//min max�Ķγ���
	param_config(min_pep_length_str, fp_index);
	param_config(max_pep_length_str, fp_index);

	min_pep_length = atoi(min_pep_length_str);
	max_pep_length = atoi(max_pep_length_str);

	//min max �Ķ�����
	param_config(min_mass_str, fp_index);
	param_config(max_mass_str, fp_index);
	min_mass = strtod(min_mass_str,NULL);
	max_mass = strtod(max_mass_str,NULL);


	printf("input path:%s\n", input_path_fasta);
	printf("output path:%s\n", output_dir_result);
	printf("enzyme type is;%d (0���⣬1��ʾ�����죬2��ʾ������)\n", enzyme_type);
	printf("enzyme terminal is : %c\n", enzyme_terminal);
	printf("enzyme amino acids: %s\n", enzyme_aa);
	printf("max miss site is: %d\n", max_miss_site);
	printf("length range of peptide is: %d  -   %d\n", min_pep_length, max_pep_length);
	printf("mass range of peptide is: %f    -   %f\n", min_mass, max_mass);

}



