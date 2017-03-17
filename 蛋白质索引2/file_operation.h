#include<string>

using namespace std;

void check_fp(FILE*fp)
{
	if (fp == NULL)
		perror("Error opening file!");

}

string new_file_name(const char*pre,const char*suff,const char*out_dir_result)
{
	char filename[300];
	//strcat(filename, out_dir_result);
	strcpy(filename, out_dir_result);
	strcat(filename, pre);
	strcat(filename, suff);
	return filename;
}

void configure(char*para, FILE*fp)
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


//配置氨基酸的质量
void an_ini_config()
{

}