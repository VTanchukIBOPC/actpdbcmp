// pdb.cpp: implementation of the Cpdb class.
//
//////////////////////////////////////////////////////////////////////


#include "pdb.h"
#include <string.h>
#include <stdlib.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
bool CpdbSite::m_bOnlyBackbone = false;

CpdbSite::CpdbSite()
{
	m_nAtoms = m_nResidues = 0;
	Atoms = NULL;
	Residues = NULL;
}

CpdbSite::CpdbSite(int nAt, int nRes)
{
	m_nAtoms = nAt;
	m_nResidues = nRes;
	AllocResidues();
	AllocAtoms();
}

void CpdbSite::operator=(CpdbSite &st){
//CpdbSite::CpdbSite(const CpdbSite &st){
	int i;
	m_nAtoms = st.m_nAtoms;
	m_nResidues = st.m_nResidues;
	m_nSiteNum = st.m_nSiteNum;
	FName = st.FName;
	Conform = st.Conform;
	AllocResidues();
	AllocAtoms();

	
	for(i = 0; i < m_nResidues; i++) Residues[i] = st.Residues[i];
	for(i = 0; i < m_nAtoms; i++){
		Atoms[i] = st.Atoms[i];
		Atoms[i].m_Res = &Residues[Atoms[i].m_Res->m_nNumInSite];
		}
	//return *this;
}

CpdbSite::~CpdbSite()
{
	if(Atoms) delete Atoms;
	delete Residues;
}


bool CPdbReader::ParseStr(char conform)
{
	CurAtom.ParseAtomStr(buff, conform);
	return (CurAtom.num > 0);
}


void CAtom::ParseAtomStr(char* str, char conform)
{
	num = 0;
	if (strstr(str, "ATOM") == str && (str[16] == conform || str[16] == ' ')) {
		Conform = str[16];
		if (str[16] != ' ') str[16] = ' ';
		sscanf(str + 6, "%7d%4s%4s %c %s %lf %lf %lf", &num, Name, m_Res->Name, m_Res->ChainName,
			m_Res->NumStr, &xyz[0], &xyz[1], &xyz[2]);
		if (Name[0] == 'H') num = 0;
		m_Res->ChainName[1] = '\0';
	}
}

int CPdbReader::ReadPDBFile(char *FName)
{
	static char Conformations[] = { ' ', 'A', 'B', 'C', 'D' };

	m_nTmpAtoms = m_nTmpSites = 0;
	CurAtom.m_Res = &CurResidue;
	CurAtom.num = 0;

	FILE *in = fopen(FName, "r");
	int n_num = 0, nRes, n, i;
	char PrevStr[8];
	bool bGetAtom;

	if (in) {
		CpdbSite* pSites = TmpSites;
		for (int cc = 0; cc < 5; cc++) {
			fseek(in, 0, SEEK_SET);
			n = ReadResidues(in, Conformations[cc]);
			CpdbSite::SortResidues(TmpResidues, n);

			int ns = Template.FindItself(TmpResidues, n, pSites);

			n = 0;
			PrevStr[0] = '\0';
			if (ns > 0) {
				bool hasLetters = false;
				fseek(in, 0, SEEK_SET);
				bGetAtom = false;
				int nAtoms = m_nTmpAtoms;
				while (fgets(buff, 250, in) != NULL) {
					buff[77] = 0;
					if (ParseStr(Conformations[cc])) {
						if (strcmp(CurResidue.NumStr, PrevStr)) {
							strcpy(PrevStr, CurResidue.NumStr);
							bGetAtom = false;
							for (i = 0; i < ns; i++)
								if (pSites[i].HasResidue(CurResidue) >= 0) {
									bGetAtom = true;
									break;
								}
							n++;
						}

						if (bGetAtom)
							for (i = 0; i < ns; i++) {
								nRes = pSites[i].HasResidue(CurResidue);
								if (nRes >= 0) {
									if (CurAtom.Conform != ' ') {
										pSites[i].Conform = CurAtom.Conform;
										hasLetters = true;
									}
									TmpAtoms[m_nTmpAtoms] = CurAtom;
									TmpAtoms[m_nTmpAtoms].m_nSiteNum = i + m_nTmpSites;
									TmpAtoms[m_nTmpAtoms].m_nAtNum = pSites[i].m_nAtoms++;
									TmpAtoms[m_nTmpAtoms++].m_Res = &pSites[i].Residues[nRes];
								}
							}
					}
				}
				if (cc > 0 && !hasLetters) {
					m_nTmpAtoms = nAtoms;
					break;
				}
				m_nTmpSites += ns;
				pSites += ns;
			}
		}
		fclose(in);
	}
	else return 1;

	for (i = 0; i < m_nTmpSites; i++) {
		TmpSites[i].AllocAtoms();
		TmpSites[i].FName = FName;
	}

	for (i = 0; i < m_nTmpAtoms; i++)
		TmpSites[TmpAtoms[i].m_nSiteNum].Atoms[TmpAtoms[i].m_nAtNum] = TmpAtoms[i];

	return 0;
}

int CpdbSite::FindItself(CResidue Res[], int nRes, CpdbSite Sites[])
{
	int maxStop, n, nSh, i, j;
	int nSRes[MAXRES]; 
	Sites[0].m_nResidues = maxStop = 0;
	for(n = i = 0; i < nRes; i++) 
		if(!strcmp(Residues[0].Name, Res[i].Name)){
			nSRes[0] = i;
			nSh = Res[i].num - Residues[0].num;
			
			for(j = 1; j < m_nResidues; j++){
				CResidue *pRes = FindResidue(Residues[j], Res, nRes, nSh, Res[i].ChainName);				
				if(pRes) nSRes[j] = pRes - Res;
				else break;
				}

			if(j >= m_nResidues){
				Sites[n].m_nResidues = m_nResidues;
				Sites[n].AllocResidues(); 
				for(j = 0; j < m_nResidues; j++){
					Sites[n].Residues[j] = Res[nSRes[j]];
					Sites[n].Residues[j].m_nNumInSite = j;
					}
				Sites[n].m_nResidues = m_nResidues;
				Sites[n].m_nAtoms = 0;
				Sites[n].m_nSiteNum = n;
				Sites[++n].m_nResidues = 0;
				}
			else if(j > maxStop) maxStop = j;
			}
		return (n > 0) ? n : -maxStop;
}



void CpdbSite::AllocAtoms()
{
	if(m_nAtoms > 0)	Atoms = new CAtom[m_nAtoms];
	else Atoms = NULL;	
}

int CpdbSite::ReadFromFile(char *FName)
{
	FILE *in = fopen(FName, "r");
	if(in){
		fscanf(in, "%d", &m_nResidues);
		Residues = new CResidue[m_nResidues];
		Atoms = NULL;
		for(int i = 0; i < m_nResidues; i++){
			fscanf(in, "%s %d", Residues[i].Name, &Residues[i].num);
			strcpy(Residues[i].ChainName, "A");
		}

		fclose(in);
		return 0;
	}
	return 1;
}

bool CAtom::operator>(CAtom &At2){
	if(m_Res->m_nNumInSite > At2.m_Res->m_nNumInSite) return true;
	if(m_Res->m_nNumInSite == At2.m_Res->m_nNumInSite){
		if(strcmp(Name, At2.Name) > 0) return true;
		}
	return false;
}


bool CAtom::operator==(CAtom &At2){
	return (m_Res->m_nNumInSite == At2.m_Res->m_nNumInSite && !strcmp(Name, At2.Name));
}

int cmpAtoms( const void *arg1, const void *arg2 )
{
    if( (*(CAtom*)arg1) == (*(CAtom*)arg2)) return 0;
    if( (*(CAtom*)arg1) > (*(CAtom*)arg2)) return 1;
	return -1;
}


void CpdbSite::SortAtoms()
{
	if(m_bOnlyBackbone) LeaveOnlyBackbone();
	qsort(Atoms, m_nAtoms, sizeof(CAtom), cmpAtoms);
}

int CPdbReader::ReadTemplate(char *FName)
{
	FILE *in = fopen(FName, "r");
	if(in){
		Template.m_nResidues = ReadResidues(in);
		Template.AllocResidues();
		fclose(in);

		for(int i = 0; i < Template.m_nResidues; i++) Template.Residues[i] = TmpResidues[i];
		return 0;
	}
	return 1;
}


int CPdbReader::ReadResidues(FILE *in, char conform)
{
	int n;
	char PrevStr[8];
	CurAtom.m_Res = &CurResidue;
	strcpy(CurResidue.ChainName, "A");
	CurAtom.num = 0;
	n = 0;
	PrevStr[0] = '\0';

	if (in) {
		while (fgets(buff, 250, in) != NULL)
			if (ParseStr(conform)) {
				if (strcmp(CurResidue.NumStr, PrevStr)) {
					strcpy(PrevStr, CurResidue.NumStr);
					CurResidue.GetNum();
					CurResidue.m_nNumInSite = n;
					TmpResidues[n++] = CurResidue;
				}
			}
	}
	return n;
}

void CpdbSite::AllocResidues()
{
	if(m_nResidues > 0)	Residues = new CResidue[m_nResidues];
	else Residues = NULL;
}

void CpdbSite::LeaveOnlyBackbone()
{
	int i, n;
	for(n = i = 0; i < m_nAtoms; i++)
		if(!Atoms[i].IsBackboneAtom()) n++;
		else if(n) Atoms[i - n] = Atoms[i];
	m_nAtoms -= n;
}

bool CAtom::IsBackboneAtom()
{
	static char *BackboneAtoms[4] = {"CA", "C", "N", "O"};
	for(int i = 0; i < 4; i++)
		if(!strcmp(Name, BackboneAtoms[i])) return true;
	return false;
}

int cmpRes( const void *arg1, const void *arg2){
	int n = strcmp(((CResidue*)arg1)->ChainName, ((CResidue*)arg2)->ChainName);
	if(!n){
		n = strcmp(((CResidue*)arg1)->Name, ((CResidue*)arg2)->Name);
		if(!n) return strcmp(((CResidue*)arg1)->NumStr, ((CResidue*)arg2)->NumStr);
	}
	return n;
}

void CpdbSite::SortResidues(CResidue ResArr[], int nRes)
{
	qsort(ResArr, nRes, sizeof(CResidue), cmpRes);	
}

int CpdbSite::HasResidue(CResidue &Res)
{
	for(int i = 0; i < m_nResidues; i++)
		if(!cmpRes(&Res, &Residues[i])) return i;
	return -1;
}


CResidue* CpdbSite::FindResidue(CResidue &Res, CResidue ResArr[], int nRes, int nShift, char *ChName)
{
	if(ChName){
		CResidue sRes = Res;
		sRes.ModifyName(nShift);
		strcpy(sRes.ChainName, ChName);
		return (CResidue*)bsearch(&sRes, ResArr, nRes, sizeof(CResidue), cmpRes);
	}
	return (CResidue*)bsearch(&Res, ResArr, nRes, sizeof(CResidue), cmpRes);
}

void CpdbSite::GetSiteName(){
	int n;
	char *p = strrchr(FName, '.');
	if(p) n = p - FName; 
	else n = strlen(FName);
	strncpy(SiteName, FName, n);
	if(m_nSiteNum > 0) sprintf(SiteName + n, "-%d", m_nSiteNum + 1);
	if (Conform != ' ') {
			int pos = strlen(SiteName);
			SiteName[pos] = Conform;
			SiteName[++pos] = '\0';
	}
}

void CResidue::ModifyName(int nSh)
{
	char buff[8], *p;
	int n = GetNum(&p);
	sprintf(buff, "%d%s", n + nSh, p);
	strcpy(NumStr, buff);
}

int CResidue::GetNum(char **pStr)
{
	int n;
	sscanf(NumStr, "%d%n", &num, &n);
	if(pStr) *pStr = NumStr + n;
	return num;
}
