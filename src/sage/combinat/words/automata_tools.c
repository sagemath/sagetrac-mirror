#include <stdio.h>
#include <stdlib.h>
//#include <system.h>

const char *pass (const char *txt)
{
	while(*txt != '{')
		txt++;
	return txt+1;
}

void TikZ (const char *data, const char *graph_name)
{
	int verb = 0;
	const char *name = "/Users/mercat/Desktop/a.dot";
	char tamp[1024];
	FILE *f = fopen(name, "w");
	if (!f)
	{
		printf("Impossible d'ouvrir le fichier a.dot !\n");
		return;
	}
	//printf("len = %d\n", len);
	//printf("%s\n", data);
	if (verb)
		printf("pass...\n");
	data = pass(data);
	if (verb)
		printf("start...\n");
	fprintf(f, "digraph %s\n{\n"\
	"	node[fontsize=20]"\
	"	edge[fontsize=20]"\
	"	rankdir = LR;\n"\
	"	size = \"10,10\";\n"\
	"	center = 1;\n"\
	"	nodesep = \"0.2\"\n", graph_name);
//	"	ranksep = \"0.4 equally\";\n", graph_name);
//	"	rotate = -90\n"\
//	"	orientation=landscape\n"\
//	"orientation = Landscape\n");
	//data = sscanf(data, "digraph {\n");
	if (verb)
		printf("write...\n");
	//fwrite(data, len, 1, f);
	fprintf(f, "%s", data);
	fclose(f);
	if (verb)
		printf("draw...\n");
	sprintf(tamp, "dot %s -Gname -Tsvg > output.svg", name);
	system(tamp);
}
