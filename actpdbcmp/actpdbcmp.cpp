// actpdbcmp.cpp : Defines the entry point for the console application.
//

#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>

#ifdef _WIN32
	#include<io.h>
	#include<windows.h>
#else
	#include<sys/types.h>
	#include<dirent.h>
#endif

#include "Compare.h"
#include "pdb.h"
#include "claster.h"
#define float double
#define MAX_SITES 2000
#define MAX_FILES 1000
#define MAX_THREADS 100

void InitCriticalSection();
void StartComparingTwoSites(int i, int j);
void CompareTwoSites(int i, int j);
void WaitForRemainingThreads();

CpdbSite pdbSites[MAX_SITES];
double clevel[MAX_SITES];
int iclson[MAX_SITES], icrson[MAX_SITES], iptr[MAX_SITES];
int iopt[2] = {1, 0};
double *dist_matr, res_impact[MAXRES], max_res_impact[MAXRES], min_res_impact[MAXRES];
int MinCase[MAXRES][2], MaxCase[MAXRES][2];
int nSites, nFiles, res_cnt[MAXRES];
char *FNames[MAX_FILES];

bool bCA = false;
int PI_steps = 3;
bool bOutPDB = false;
char TemplateFile[260];
char pdbext[400];

int nMaxThreads = 4;

char *Params[] = { "Template", "CA", "PI_steps", "OutPDB", "Threads", "OnlyBackbone", "pdbext"};


void print_site_info(int argc, char* argv[]){
	int i, j, n, n2;
	FILE *outf = fopen("sites.txt", "w");
	for(i = 0; i < argc; i++) fprintf(outf, " %s", argv[i]);
	fprintf(outf, "\n%d\n%d\n", nFiles, nSites);
	for(n = i = 0; i < nSites; i++){
		if(!i || strcmp(pdbSites[i].FName, pdbSites[i - 1].FName)){
			fprintf(outf, "%d %s\n", ++n, pdbSites[i].FName);
			n2 = 0;
			}
		fprintf(outf, "site %d : atoms %d %s\n", i + 1, pdbSites[i].m_nAtoms, pdbSites[i].SiteName);
		for(j = 0; j < pdbSites[i].m_nResidues; j++)
			fprintf(outf, "%s %d %s\n", pdbSites[i].Residues[j].Name,
			  pdbSites[i].Residues[j].num,
			  pdbSites[i].Residues[j].ChainName);
		}
	fclose(outf);
}

void ReadConfigFile(char *FName) {
	FILE* inf = fopen(FName, "r");
	char buff[300], *pos;


	while (fgets(buff, 300, inf) != NULL)
	{
		int len = strlen(buff);
		if (len > 2) {
			if (buff[len - 1] == '\n') buff[len - 1] = '\0';
			for (int i = 0; i < (sizeof(Params)/sizeof(char*)); i++)
				if ((pos = strstr(buff, Params[i])) != NULL) {
					pos = strchr(pos, '=');
					if (pos != NULL) {
						pos++;
						switch (i) {
							case 0: strcpy(TemplateFile, pos); break;
							case 1: if(strchr(pos, '1') != NULL) bCA = true; break;
							case 2: sscanf(pos, "%d", &PI_steps);  break;
							case 3: if (strchr(pos, '1') != NULL) bOutPDB = true; break;
							case 4: sscanf(pos, "%d", &nMaxThreads);  break;
							case 5: if (strchr(pos, '1') != NULL) CpdbSite::m_bOnlyBackbone = true; break;
#ifdef _WIN32
							case 6: strcpy(pdbext, "*."); strcpy(pdbext + 2, pos); break;
#else
							case 6:  strcpy(pdbext, "."); strcpy(pdbext + 1, pos); break;
#endif
							default: ;
						}
						break;
					}
				}
		}
	}
	fclose(inf);
	if (nMaxThreads < 1) nMaxThreads = 1;
	if (PI_steps < 1) PI_steps = 1;
#ifdef _WIN32
	if(pdbext[0] == '\0') strcpy(pdbext, "*.pdb");
#else
	if (pdbext[0] == '\0') strcpy(pdbext, ".pdb");
#endif
}

void print_dist_matr(){
	FILE *outf = fopen("dist_matr.txt", "w");
	int i, j, n;
	double sum_dist[MAX_SITES];
	double mean = 0;

	fprintf(outf, "%d", nSites);
	for(i = 0; i < nSites; i++) sum_dist[i] = 0;
	for(n = i = 0; i < nSites; i++, n++){
		for(j = 0; j < i; j++){
			fprintf(outf, "%-9.5f ", dist_matr[n]);
			sum_dist[i] += dist_matr[n];
			sum_dist[j] += dist_matr[n];
			mean += dist_matr[n++];
			}
		fprintf(outf, "\n");
		}
	i = (nSites * (nSites - 1)) / 2;
	fprintf(outf, "count = %d\nmean = %-9.4lf\n", i, mean / i);
	for(i = j = 0; i < nSites; i++)
		if(sum_dist[i] < sum_dist[j]) j = i;
	fprintf(outf, "centroid = %d %-9.4lf\n", j + 1, sum_dist[j] / (nSites - 1));
	fclose(outf);
}


void print_clust_info(){
	FILE *outf = fopen("cluster.txt", "w");
	int i, n = nSites - 1;
	fprintf(outf, "%d\n", n);
	for(i = 0; i < n; i++)
		fprintf(outf, "%d %d %-9.5f\n", iclson[i], icrson[i], clevel[i]);
	fclose(outf);
}

int cmpImp( const void *arg1, const void *arg2){
	return (res_impact[*(int*)arg1] > res_impact[*(int*)arg2]) ? -1 : 
			(res_impact[*(int*)arg1] < res_impact[*(int*)arg2]) ? 1 : 0;
}

void print_res_impact(CpdbSite &TplSite, float num){
	FILE *outf = fopen("res_impact.txt", "w");
	int nums[MAXRES];
	int n, i, n1;
	double mean = 0;
	int cnt = 0;

	for(i = 0; i < TplSite.m_nResidues; i++){
		nums[i] = i;
		if(res_cnt[i]){
			mean += res_impact[i];
			cnt += res_cnt[i];
			res_impact[i] = sqrt(res_impact[i] / res_cnt[i]);
			}
		}
	qsort(nums, TplSite.m_nResidues, sizeof(int), cmpImp);
	fprintf(outf, "%d\n", TplSite.m_nResidues);
	n1 = -1;
	for(i = 0; i < TplSite.m_nResidues; i++){
		n = nums[i];
		if(res_cnt[n]){
			if(n1 < 0) n1 = n;
			fprintf(outf, "%s%s\t%-8.2lf\t%-8.1f\t%d\t%-8.1f\t%-8.2f\t%s\t%s\t%-8.2f\t%s\t%s\n", TplSite.Residues[n].Name,
				TplSite.Residues[n].NumStr, res_impact[n], (res_impact[n] * 100) / res_impact[n1], res_cnt[n], res_cnt[n] / num,
				max_res_impact[n], pdbSites[MaxCase[n][0]].SiteName, pdbSites[MaxCase[n][1]].SiteName,
				min_res_impact[n], pdbSites[MinCase[n][0]].SiteName, pdbSites[MinCase[n][1]].SiteName);
			}
		}
	fprintf(outf, "count = %d\nmean = %-9.4lf\n", cnt, sqrt(mean / cnt));
	fclose(outf);
}

int compare_fname( const void *arg1, const void *arg2 )
{
   /* Compare all of both strings: */
   return strcmp( * ( char** ) arg1, * ( char** ) arg2 );
}

#ifdef _WIN32
CRITICAL_SECTION CriticalSection;

HANDLE  hThreadArray[MAX_THREADS];
DWORD   dwThreadIdArray[MAX_THREADS];

DWORD WINAPI ThreadFunction(LPVOID lpParam) {
	int n = (int)lpParam;
	CompareTwoSites(HIWORD(n), LOWORD(n));
	return 0;
}


int WaitFunction() {
	int i;
	for (i = 0; i < nMaxThreads; i++)
		if (hThreadArray[i] == 0) break;
	if (i < nMaxThreads) return i;

	DWORD n = WaitForMultipleObjects(nMaxThreads, hThreadArray, FALSE, INFINITE);
	if (n >= WAIT_OBJECT_0 && n < WAIT_OBJECT_0 + nMaxThreads) {
		n -= WAIT_OBJECT_0;
		CloseHandle(hThreadArray[n]);
		hThreadArray[n] = 0;
		return n;
	}
	return -1;
}

void InitCriticalSection() {
	InitializeCriticalSection(&CriticalSection);
}

/*
void WaitForRemainingThreads() {
	int i, n;
	for (n = i = 0; i < nMaxThreads; i++)
		if (hThreadArray[i])
			hThreadArray[i - n] = hThreadArray[i];
		else ++n;

	n = nMaxThreads - n;

	WaitForMultipleObjects(n, hThreadArray, TRUE, INFINITE);
	for (i = 0; i < n; i++) CloseHandle(hThreadArray[i]);
	DeleteCriticalSection(&CriticalSection);
}
*/

void WaitForRemainingThreads() {
	WaitForMultipleObjects(nMaxThreads, hThreadArray, TRUE, INFINITE);
	for (int i = 0; i < nMaxThreads; i++)
		if (hThreadArray[i] != 0)
			CloseHandle(hThreadArray[i]);
	DeleteCriticalSection(&CriticalSection);
}

void StartComparingTwoSites(int i, int j) {

	int nThread = WaitFunction();

	hThreadArray[nThread] = CreateThread(
		NULL,                   // default security attributes
		0,                      // use defalt stack size = 1M
		ThreadFunction,			// thread function name
		(LPVOID)((i << 16) + j),   // argument to thread function 
		0,                      // use default creation flags 
		&dwThreadIdArray[nThread]);   // returns the thread identifier 
}

#else
#endif

void CompareTwoSites(int i, int j) {
	CCompare comp1, comp2;
	double d1, d2;
	int n1, n2, n = ((i * (i + 1)) / 2) + j;

	comp1.bCA = bCA;
	comp1.PI_steps = PI_steps;

	comp2.bCA = bCA;
	comp2.PI_steps = PI_steps;

	n1 = comp1.GetAtomPairs(pdbSites[i], pdbSites[j]);
	d1 = comp1.SobOptim();
	n2 = comp2.GetAtomPairs(pdbSites[j], pdbSites[i]);
	d2 = comp2.SobOptim();

	if (bOutPDB) {
		comp1.SaveToFile(i, j, pdbSites[i], pdbSites[j]);
		comp2.SaveToFile(i, j, pdbSites[j], pdbSites[i]);
	}

#ifdef _WIN32
	EnterCriticalSection(&CriticalSection);
#else
#endif

	//printf("%-8.3f %-8.3f\n", d1, d2);
	if (d2 < d1) {
		dist_matr[n] = d2;
		n1 = n2;
		comp2.AddToGlobalImpact(res_impact, max_res_impact, min_res_impact, res_cnt, MinCase, MaxCase, i, j);
		//if(bOutPDB) comp2.SaveToFile(i, j, pdbSites[j], pdbSites[i]);
	}
	else {
		dist_matr[n] = d1;
		comp1.AddToGlobalImpact(res_impact, max_res_impact, min_res_impact, res_cnt, MinCase, MaxCase, i, j);
		//if(bOutPDB) comp1.SaveToFile(i, j, pdbSites[i], pdbSites[j]);
	}
#ifdef _WIN32
	LeaveCriticalSection(&CriticalSection);
#else
#endif
}

int main(int argc, char* argv[])
{
	int i, j, n, st, st_all;
	bool bNums = false;

	if (argc < 2) {
		printf("There must be an agrument with config file name!\n");
		return 1;
	}

	ReadConfigFile(argv[1]);
	if (TemplateFile[0] == '\0') {
		printf("There must be template file name (Template=) in the config file!\n");
		return 2;
	}

#ifdef _WIN32
	struct _finddata_t pdb_file;
    long hFile;

	if ((hFile = _findfirst(pdbext, &pdb_file)) == -1L) {
		fprintf(stderr, "No %s files in current directory!\n", pdbext);
		return 3;
	}
	else{
		do{
			if(_stricmp(TemplateFile, pdb_file.name)){
				FNames[nFiles] = new char[strlen(pdb_file.name) + 1];
				strcpy(FNames[nFiles++], pdb_file.name);
				}
			}while( _findnext( hFile, &pdb_file ) == 0 );
		_findclose( hFile );
	}
#else
    DIR *dp;
    struct dirent *dirp;

    if((dp  = opendir(".")) == NULL) {
        return -1;
    }

    while ((dirp = readdir(dp)) != NULL) {
		if(strstr(dirp->d_name, pdbext) == dirp->d_name + strlen(dirp->d_name) - 4){
			if (stricmp(TemplateFile, dirp->d_name)) {
				FNames[nFiles] = new char[strlen(dirp->d_name) + 1];
				strcpy(FNames[nFiles++], dirp->d_name);
			}
		}
    }

    closedir(dp);
#endif
	CPdbReader *pdbRead = new CPdbReader();
	pdbRead->ReadTemplate(TemplateFile);
	qsort(FNames, nFiles, sizeof(char *), compare_fname);
	FILE *outfb = fopen("bad_files.txt", "w");
	for(i = 0; i < nFiles; i++){
		if(!pdbRead->ReadPDBFile(FNames[i])){
			if(pdbRead->m_nTmpSites <= 0){
					printf("File %d: %s does not have sites!\n", i + 1, FNames[i]);
					fprintf(outfb, "%d: %s\n", i + 1, FNames[i]);
					j = -pdbRead->m_nTmpSites;
					if(j > 0){
						CResidue prs = pdbRead->Template.Residues[j++];
						printf("Problem residue %d: %s%s\n", j, prs.Name, prs.NumStr);
						fprintf(outfb, "Problem residue %d: %s%s\n", j, prs.Name, prs.NumStr);
						}
					}
			for(j = 0; j < pdbRead->m_nTmpSites; j++)
					pdbSites[nSites++] = pdbRead->TmpSites[j];
			}
	}
	fclose(outfb);
		
	for(i = 0; i < nSites; i++) pdbSites[i].SortAtoms();
	for (i = 0; i < nSites; i++) pdbSites[i].GetSiteName();
	print_site_info(argc, argv);

	st_all = (nSites * (nSites + 1)) / 2;
	dist_matr = new double[st_all];
	st_all -= nSites;
	n = 0;
	printf("Comparing %d sites\n", nSites);

	for(i = 0; i < pdbRead->Template.m_nResidues; i++) min_res_impact[i] = 1e6;
	//for(i = 0; i < nSites; i++) pdbSites[i].GetSiteName();

	InitCriticalSection();
	n = st = 0;
	for(i = 0; i < nSites; i++){
		for(j = 0; j < i; j++, n++){	
			StartComparingTwoSites(i, j);
			//CompareTwoSites(i, j); //DEBUG
			//printf("%d - %d = %-9.4f pairs %d\n", i + 1, j + 1, dist_matr[n++], n1);
			if(++st % 10 == 0) fprintf(stderr, "%d of %d\r", st, st_all);
		}
		dist_matr[n++] = 0; //diagonal 
	}
	WaitForRemainingThreads();
	fprintf(stderr, "%d of %d\n", st, st_all);
	print_dist_matr();

	clust(nSites, iopt, dist_matr, clevel, iclson, icrson, iptr);
	delete dist_matr;

	print_clust_info();
	print_res_impact(pdbRead->Template, st_all);
	delete pdbRead;

	return 0;
}

