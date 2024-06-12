// aclust.cpp : Defines the entry point for the console application.
//

#include <string.h>
#include <stdio.h>


#define MAXCLUST 5000


int nClust, nObj, nAllClust;

bool *bClust_obj[MAXCLUST];

int lson[MAXCLUST], rson[MAXCLUST];
bool UsedClust[MAXCLUST], centr[MAXCLUST];
float cl_levels[MAXCLUST];

char SiteName[MAXCLUST][16];
int SiteNum[MAXCLUST];

float *dist[MAXCLUST / 2];

int find_param(int argc, char* argv[], char *key, bool bOne){
	for(int i = 1; i < argc; i++)
		if(!strcmp(key, argv[i])){
			if(bOne) return i;
			if(i < argc - 1) return i + 1;
			return 0;
			}
	return 0;
}

void Init_ClustObj(){
	int i, j;
	bool *pb = new bool[nAllClust * nObj];
	for(i = 0; i < nAllClust; i++){
		bClust_obj[i] = pb;
		for(j = 0; j < nObj; j++, pb++) *pb = false;
		}
	for(i = 0; i < nObj; i++) bClust_obj[i][i] = true;

}

void add_son(int nClust, int nSon){
	if(nSon >= nObj){
		for(int i = 0; i < nObj; i++)
			if(bClust_obj[nSon][i]) bClust_obj[nClust][i] = true;
	}
	else bClust_obj[nClust][nSon] = true; 
}

void Fill_ClustObj(){
	for(int i = nObj; i < nAllClust; i++){
		add_son(i, lson[i]);
		add_son(i, rson[i]);
		}
}

int read_clust(FILE *in){
	int i, j;

	fscanf(in, "%d", &nClust);
	nObj = nClust + 1;
	nAllClust = nObj + nClust;
	//Init_ClustObj();
	for(i = 0; i < nClust; i++){
		j = i + nObj;
		fscanf(in, "%d %d %f", &lson[j], &rson[j], &cl_levels[j]);
		lson[j]--;
		rson[j]--;
		}
	Init_ClustObj();
	Fill_ClustObj();
	return nClust;
}

void get_clust_setN(int n){
	int i, cnt;

	for(i = 1; i < nAllClust; i++) SiteNum[i] = -1;
	for(i = 0; i < nAllClust; i++) UsedClust[i] = false;

	cnt = 0;
	for(i = nAllClust - 1; i >= 0; i--){
		UsedClust[lson[i]] = true;
		UsedClust[rson[i]] = true;

		if(UsedClust[i]) UsedClust[i] = false;
		else ++cnt;
		if(++cnt >= n) break;
		}
}

void get_clust_setLevel(double level){
	int i, cnt;

	for(i = 1; i < nAllClust; i++) SiteNum[i] = -1;
	for(i = 0; i < nAllClust; i++) UsedClust[i] = false;

	cnt = 0;
	for(i = nAllClust - 1; i >= 0; i--)
		if(cl_levels[i] >= level){
			UsedClust[lson[i]] = true;
			UsedClust[rson[i]] = true;

			if(UsedClust[i]) UsedClust[i] = false;
		}
}

int read_sites(char *FName){
	int n, i;
	FILE *in = fopen(FName, "r");
	if(in){
		char pbuff[256], buff[256], tmp[8];
		fgets(buff, 255, in);
		fgets(buff, 255, in);
		
		fscanf(in, "%d", &n);
		for(i = 0; i < n; i++){
			do{
				strcpy(pbuff, buff);
				fgets(buff, 255, in);
			}while(strstr(buff, "site") != buff);

			if(pbuff[0] >= '0' && pbuff[0] <= '9'){
				sscanf(pbuff, "%s %s", tmp, SiteName[i]);
				SiteNum[i] = 0;
				}
			else{
				strcpy(SiteName[i], SiteName[i - 1]);
				SiteNum[i] = SiteNum[i - 1] + 1;
				}

			}
		fclose(in);

		for(i = 1; i < n; i++)
			if(SiteNum[i] > 0) SiteNum[i - 1]++;
			else if(SiteNum[i - 1] > 0) SiteNum[i - 1]++;
		}
	return n;
}

void init_dist(){
	float *p = new float[(nObj * (nObj + 1)) / 2];
	for(int i = 0; i < nObj; i++){
		dist[i] = p;
		p += i + 1;
	}
}

int read_dist(char *FName){
	FILE *in = fopen(FName, "r");
	if(in){
		int n, i, j;
		fscanf(in, "%d", &n);
		if(n != nObj){
			fclose(in);
			return -1;
			}
		init_dist();
		for(i = 1; i < n; i++){
			for(j = 0; j < i; j++)
				fscanf(in, "%f", &dist[i][j]);
			dist[i][i] = 0;
		}
		fclose(in);
		return n;
	}
	return 0;
}

int calc_centr(bool *clust, double *cdist, double *avdist){
	int clust_nums[MAXCLUST];
	int i, j, n, i1, nc;
	double d, min_dist, sum;

	for(n = i = 0; i < nObj; i++) 
		if(clust[i]) clust_nums[n++] = i;

	nc = -1;
	min_dist = 1e20;
	sum = 0;
	for(i = 0; i < n; i++){
		d = 0;
		i1 = clust_nums[i];
		for(j = 0; j < i; j++) d += dist[i1][clust_nums[j]];
		for(j++; j < n; j++) d += dist[clust_nums[j]][i1];
		sum += d;
		if(d < min_dist){
			min_dist = d;
			nc = i1;
			}
		}
	if(nc >= 0){
		*cdist = min_dist / (n - 1);
		*avdist = sum / (n * (n - 1));
		}
	return nc;
}

void output_clust(FILE *outf, bool bCentr){
	int n, nc, i, j;
	double dst, mdist;

	for(n = i = 0; i < nAllClust; i++) if(UsedClust[i]) n++;
	fprintf(outf, "%d clusters\n", n);
	for(i = 0; i < nAllClust; i++)
		if(UsedClust[i]){
			for(n = j = 0; j < nObj; j++) 
				if(bClust_obj[i][j]) n++;

			if(n >= 2) fprintf(outf, "cluster %d - %d level %-9.4f", i + 1, n, cl_levels[i]);
			else fprintf(outf, "cluster %d - 1", i + 1);

			nc = -1;
			if(n > 2){
				if(bCentr){
					nc = calc_centr(bClust_obj[i], &dst, &mdist);
					fprintf(outf, "mean %-9.4f", mdist);
					}
				}
			fprintf(outf, "\n");

			for(n = j = 0; j < nObj; j++) 
				if(bClust_obj[i][j]){ 
					if(SiteNum[j] < 0) fprintf(outf, "%d", j + 1);
					else if(SiteNum[j] > 0) fprintf(outf, "%d: %s-%d", j + 1, SiteName[j], SiteNum[j]);
					else fprintf(outf, "%d: %s", j + 1, SiteName[j]);

					if(j == nc) fprintf(outf, "\t\t dist 0.0000 centroid %-9.4lf\n", dst);
					else if(nc >= 0) fprintf(outf, "\t\t dist %-9.4lf\n", (j > nc) ? dist[j][nc] : dist[nc][j]);
					else fprintf(outf, "\n"); 

					rson[j] = i;
					lson[j] = nc; // centroid
				}
			}
}


void output_clust2(FILE *outf, bool bCentr){
	int i, j;

	fprintf(outf, "\n\n");
	for(i = 0; i < nObj; i++){
		j = i;

		if(SiteNum[j] > 0) fprintf(outf, "%d: %s-%d", j + 1, SiteName[j], SiteNum[j]);
		else fprintf(outf, "%d: %s", j + 1, SiteName[j]);

		j = lson[i]; //cntroid

		fprintf(outf, "\t\t");
		if(SiteNum[j] > 0) fprintf(outf, "%s-%d\t%d", SiteName[j], SiteNum[j], j + 1);
		else fprintf(outf, "%s\t%d", SiteName[j], j + 1);	
		fprintf(outf, "\n");
		}
}

char* Params[] = { "Cluster", "Sites", "Matrix", "Output", "Level", "Number"};
char ClusterFile[300], SitesFile[300], MatrixFile[300], OutFile[300];
float ClustLevel;
int NumClust;


void ReadConfigFile(char* FName) {
	char buff[400], * pos;
	FILE* inf = fopen(FName, "r");
	if (inf) {
		while (fgets(buff, 400, inf) != NULL)
		{
			int len = strlen(buff);
			if (len > 2) {
				if (buff[len - 1] == '\n') buff[len - 1] = '\0';
				for (int i = 0; i < (sizeof(Params) / sizeof(char*)); i++)
					if ((pos = strstr(buff, Params[i])) != NULL) {
						pos = strchr(pos, '=');
						if (pos != NULL) {
							pos++;
							switch (i) {
							case 0: strcpy(ClusterFile, pos); break;
							case 1: strcpy(SitesFile, pos); break;
							case 2: strcpy(MatrixFile, pos); break;
							case 3: strcpy(OutFile, pos); break;
							case 4: sscanf(pos, "%f", &ClustLevel);  break;
							case 5: sscanf(pos, "%d", &NumClust);  break;

							default:;
							}
							break;
						}
					}
			}
		}
		fclose(inf);
	}
	if (ClustLevel < 1e-6 && NumClust == 0) NumClust = 5;
	if (ClusterFile[0] == '\0') strcpy(ClusterFile, "cluster.txt");
	if (SitesFile[0] == '\0') strcpy(SitesFile, "sites.txt");
	if (MatrixFile[0] == '\0') strcpy(MatrixFile, "dist_matr.txt");
}

int main(int argc, char* argv[])
{
	int i;

	if (argc < 2) {
		printf("There must be an agrument with config file name!\n");
		return 1;
	}

	ReadConfigFile(argv[1]);

	FILE *in = fopen(ClusterFile, "r");
	if(in){
		read_clust(in);
		fclose(in);

		i = find_param(argc, argv, "-l", false);
		if(ClustLevel > 1e-6)
			get_clust_setLevel(ClustLevel);
		else
			get_clust_setN(NumClust);

		read_sites(SitesFile);

		read_dist(MatrixFile);

		FILE *outf;
		if (OutFile[0] != '\0') outf = fopen(OutFile, "w");
		else outf = stdout;

		output_clust(outf, (dist[0] != NULL));
		output_clust2(outf, (dist[0] != NULL));

		if(outf != stdout) fclose(outf);
		}
	else printf("No input file!\n");

	//scanf("%d", &num_clust);
	_fcloseall();

	return 0;
}
