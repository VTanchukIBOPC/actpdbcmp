// Compare.cpp: implementation of the CCompare class.
//
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "Compare.h"
#include "PowellOptim.h"
#include "Sobseq.h"
#include <string.h>
#include <stdio.h>
#ifdef _WIN32
	#include <direct.h>
#else
	#include <sys/types.h>
	#include <sys/stat.h>
#endif

#define PI 3.141592653589793


void calc_move(void *p, int n, float par[], float *y, float dyda[], int na){
*y = ((CCompare*)p)->CalcMove(n, par, &dyda[1]);
}


CCompare::CCompare()
{
	int i;
	double *pd = fcov;
	for(i = 0; i < 7; i++, pd += 7) covar[i] = pd;
	pd = falf;
	for(i = 0; i < 7; i++, pd += 7) alpha[i] = pd;
	PI_steps = 3;	
	bCA = false;
}

CCompare::~CCompare()
{

}

double CAtomPair::CalcPairMove(CTransform &tra, double dyda[], bool bSave)
{
	double txyz[3], der[3][3], res;
	int i, j;

	for(i = 0; i < 3; i++) txyz[i] = xyz2[i];

	tra.Rotate(txyz, dyda ? der : NULL);
	tra.Move(txyz);
	
	res = 0;
	if(dyda){
		double tmp[3];
		for(i = 0; i < 3; i++){
			tmp[i] = (xyz1[i] - txyz[i]);
			dyda[i + 3] = -2.0 * tmp[i];
			res += tmp[i] * tmp[i];
			}
			
		for(i = 0; i < 3; i++){
			dyda[i] = 0;
			for(j = 0; j < 3; j++)
				dyda[i] += dyda[j + 3] * der[j][i];
			}
		}
	else{
		for(i = 0; i < 3; i++)
			res += (xyz1[i] - txyz[i]) * (xyz1[i] - txyz[i]);
		}
	if(bSave){
		for(i = 0; i < 3; i++) xyz2[i] = txyz[i];
		}
	diff = res;
	return res;

}


void CTransform::SimpleRotate(double xyz[])
{
	double x1, y1, z1;

	x1 = xyz[0] * cz - xyz[1] * sz;
	y1 = xyz[0] * sz + xyz[1] * cz;

	xyz[1] = y1 * cx - xyz[2] * sx;
	z1 = y1 * sx + xyz[2] * cx;

	xyz[0] = z1 * sy + x1 * cy;
	xyz[2] = z1 * cy - x1 * sy;
}


void CTransform::Rotate(double xyz[], double der[][3])
{
	double x0, y0, z0, x1, y1, z1;

	x0 = xyz[0];
	y0 = xyz[1];
	z0 = xyz[2];

	x1 = x0 * cz - y0 * sz;
	y1 = x0 * sz + y0 * cz;

	xyz[1] = y1 * cx - z0 * sx;
	z1 = y1 * sx + z0 * cx;

	xyz[0] = z1 * sy + x1 * cy;
	xyz[2] = z1 * cy - x1 * sy;

	if(der){
		double dx1dAz, dy1dAz, dz1dAx;
		dx1dAz = -x0 * sz - y0 * cz;
		dy1dAz = x0 * cz - y0 * sz;
		dz1dAx = y1 * cx - z0 * sx;

		der[0][0] = dz1dAx * sy;
		der[0][1] = z1 * cy - x1 * sy;
		der[0][2] = dx1dAz * cy;

		der[1][0] = -y1 * sx - z0 * cx;
		der[1][1] = 0;
		der[1][2] = dy1dAz * cx;

		der[2][0] = dz1dAx * cy;
		der[2][1] = -z1 * sy - x1 * cy;
		der[2][2] = -dx1dAz * sy;
		}
}

void CTransform::Move(double xyz[])
{
	for(int i = 0; i < 3; i++) xyz[i] += dxyz[i];
}

void CTransform::Init(double par[])
{
	sz = sin(par[2]);
	cz = cos(par[2]);
	sx = sin(par[0]);
	cx = cos(par[0]);
	sy = sin(par[1]);
	cy = cos(par[1]);
	for(int i = 0; i < 3; i++) dxyz[i] = par[i + 3];
}

/*
double CCompare::CalcMove(int n, double par[], double dyda[])
{
	int i, j;
	if(n == 1){
		tra.Init(&par[1]);
		double par2[6];

		for(i = 0; i < 6; i++){
			for(j = 0; j < 6; j++)	par2[j] = par[j + 1];
			par2[i] += 0.01;
			dtra[i].Init(par2);
			}
		}
	double d2[6], res = Pairs[n - 1].CalcMove(tra, d2, false);
	for(i = 0; i < 6; i++)
		//dyda[i + 1] = d2[i];
		dyda[i + 1] = (Pairs[n - 1].CalcMove(dtra[i], NULL, false) - res) / 0.01;
	//dyda[6] = d2[3];
	//dyda[4] = d2[5];
	return res;
}
*/

double CCompare::CalcMove(int n, double par[], double dyda[], bool bSave)
{
	if(n == 1) tra.Init(&par[1]);
	return Pairs[n - 1].CalcPairMove(tra, dyda, bSave);
}

int CCompare::GetAtomPairs(CpdbSite &s1, CpdbSite &s2)
{
	int i, j, k;
	m_nResidues = m_nAtomPairs = i = j = 0;
	for(; i < s1.m_nAtoms; i++){
		while(j < s2.m_nAtoms && s1.Atoms[i] > s2.Atoms[j]) j++;
		if(j < s2.m_nAtoms && s1.Atoms[i] == s2.Atoms[j]){
			for(k = 0; k < 3; k++){
				Pairs[m_nAtomPairs].xyz1[k] = s1.Atoms[i].xyz[k];
				Pairs[m_nAtomPairs].xyz2[k] = s2.Atoms[j].xyz[k];
				}
			strcpy(Pairs[m_nAtomPairs].Name, s1.Atoms[i].Name);
			Pairs[m_nAtomPairs++].res_num = s1.Atoms[i].m_Res->m_nNumInSite;
			if(s1.Atoms[i].m_Res->m_nNumInSite >= m_nResidues) 
				m_nResidues = s1.Atoms[i].m_Res->m_nNumInSite + 1;
			j++;
			}
		}
	PrepareSymmetric(s1);
	return m_nAtomPairs;
}

double CCompare::Optim(){
	double res, best, step = PI / PI_steps;
	int i, j, k, l;
	PrepareToOptim();

	best = 1e20;
	for(i = 0;  i < PI_steps; i++){
		
		for(j = 0; j < PI_steps; j++){
			for(k = 0; k < PI_steps; k++){
				param[1] = i * step;
				param[2] = j * step;
				param[3] = k * step;
				res = RealOptim2();
				if(res < best){
					for(l = 1; l < 7; l++) best_param[l] = param[l];	
					best = res;
					}
			}
		}
	}
	
	return CalcRMSD(best_param, NULL, true);
}

double CCompare::SobOptimReal(int nSteps){
	double res, best;
	int i, j, n1 = (-1), n2 = 3;
	double p[4];
	CSobseq sob;

	sob.sobseq(&n1, p);
	best = 1e20;
	for(i = 0;  i < nSteps; i++){
		sob.sobseq(&n2, p);
		for(j = 1; j <= 3; j++) param[j] = p[j] * 2 * PI;
		res = RealOptim2();
		if(res < best){
			for(j = 1; j <= 6; j++) best_param[j] = param[j];	
			best = res;
			}
	}
	
	return CalcRMSD(best_param, NULL);
}

double CCompare::SobOptim(){
	int i, n1 = (-1), n2 = 3;

	PrepareToOptim();

	if(bCA){
		for(i = 0; i < m_nAtomPairs; i++) 
			bUsePair[i] = (Pairs[i].Name[0] == 'C' && Pairs[i].Name[1] == 'A');
		SobOptimReal(10 * PI_steps);

		for(i = 1; i <= 6; i++) param[i] = best_param[i];
		for(i = 0; i < m_nAtomPairs; i++) bUsePair[i] = true;

		RealOptim2(); 
		for(i = 1; i <= 6; i++) best_param[i] = param[i];
		}
	else{
		for(i = 0; i < m_nAtomPairs; i++) bUsePair[i] = true;
		SobOptimReal(5 * PI_steps);
		}
	return CalcRMSD(best_param, NULL, true);
}


double rmsd_func(void* pobj, double par[]){
	return ((CCompare*)pobj)->CalcRMSD(par, NULL, NULL);
}

/*
void rmsd_dfunc(double par[], double der[]){
	cur_comp->CalcRMSD(par, der, NULL);
}
*/

double CCompare::RealOptim2(){
	double fret;
	int iter, i, j;

	for(i = 1; i <= 6; i++){
		for(j = 1; j <= 6; j++) alpha[i][j] = 0;
		alpha[i][i] = 1.0;
		}
	CPowellOptim powell(this, rmsd_func);
	powell.powell(param, alpha, 6, 1.0e-6, &iter, &fret);
	return fret;
}

void CCompare::PrepareToOptim()
{
	int i;

	CenterBoth();

	for(i = 1; i <= 6; i++){
		ia[i] = 1;
		param[i] = 0;
		}
	//for(i = 4; i <= 6; i++) ia[i] = 0;

	for(i = 1; i <= m_nAtomPairs; i++){
		x[i] = i;
		y[i] = 0;
		sig[i] = 1.0;
		bUsePair[i - 1] = true;
		}

}

void CCompare::CenterBoth()
{
	int i, j;

	double cxyz1[3], cxyz2[3];
	for(j = 0; j < 3; j++){
			cxyz1[j] = 0;
			cxyz2[j] = 0;
			}
	for(i = 0; i < m_nAtomPairs; i++)
		for(j = 0; j < 3; j++){
			cxyz1[j] += Pairs[i].xyz1[j];
			cxyz2[j] += Pairs[i].xyz2[j];
			}

	for(j = 0; j < 3; j++)
			cxyz2[j] = (cxyz1[j] - cxyz2[j]) / m_nAtomPairs;
	
	for(i = 0; i < m_nAtomPairs; i++)
		for(j = 0; j < 3; j++)
			Pairs[i].xyz2[j] += cxyz2[j];

	for(i = 0; i < m_nAltPairs; i++)
		for(j = 0; j < 3; j++){
			AltPairs[i][0].xyz2[j] += cxyz2[j];
			AltPairs[i][1].xyz2[j] += cxyz2[j];
			}
}

double CCompare::CalcRMSD(double par[], double der[], bool bFinal)
{
	int nPairs, i, j, k;
	double dyda[7], res = 0;
	if(!bUsePair[0]) tra.Init(&par[1]);
	if(der)
		for(i = 1; i <= 6; i++) der[i] = 0;
		if(bFinal){
			for(i = 0; i < m_nResidues; i++){
				ResidueImpact[i] = 0;
				ResCnt[i] = 0;
			}
		}

	for(nPairs = i = 0; i < m_nAtomPairs; i++)
		if(bUsePair[i]){	
			if(der){
				res += CalcMove(i + 1, par, dyda, bFinal);
				for(j = 1; j <= 6; j++) der[j] += dyda[j - 1];
			}
			else res += CalcMove(i + 1, par, NULL, bFinal);
			++nPairs;
		}

	if(!der){
		int n1, n2, n2c, n1c = -1;
		double PairImpact = 0;
		bool bConjugPairFirst;

		for(i = 0; i < m_nAltPairs; i++){
			bConjugPairFirst = false;
			n1 = AltPairs[i][0].nAltPair;
			n2 = AltPairs[i][1].nAltPair;

			if(n1c < 0){
				if(n1 < 0){
					n1 = -(n1 + 1); //conjugated pair
					n1c = n1;
					n2c = n2;
					bConjugPairFirst = true;
					}
			}


			if(!bUsePair[n1]) Pairs[n1].CalcPairMove(tra, NULL, bFinal);
			if(!bUsePair[n2]) Pairs[n2].CalcPairMove(tra, NULL, bFinal);
			PairImpact += AltPairs[i][0].CalcPairMove(tra, NULL, bFinal) + AltPairs[i][1].CalcPairMove(tra, NULL, bFinal) -
						 Pairs[n1].diff - Pairs[n2].diff;
			
			if(!bConjugPairFirst){
				if(PairImpact < 0){
					res += PairImpact;
					if(bFinal){
						k = i;
						for(;;){
							for(j = 0; j < 3; j++){
								Pairs[n1].xyz2[j] = AltPairs[k][0].xyz2[j];
								Pairs[n2].xyz2[j] = AltPairs[k][1].xyz2[j];
								}
							Pairs[n1].diff = AltPairs[k][0].diff;
							Pairs[n2].diff = AltPairs[k][1].diff;
							if(n1c >= 0){
								n1 = n1c;
								n2 = n2c;
								n2c = n1c = -1;
								k--;
								}
							else break;
						}
					}
				}
			n2c = n1c = -1;
			PairImpact = 0;
			}
		}
	}

	if(bFinal){
		for(j = i = 0; i < m_nAtomPairs; i++)
			if(bUsePair[i]){
				ResidueImpact[Pairs[i].res_num] += Pairs[i].diff;
				ResCnt[Pairs[i].res_num]++;
			}
		}

	return sqrt(res / nPairs);
}

void CCompare::AddToGlobalImpact(double GlobResImp[], double MaxImpact[], double MinImpact[], int GlobResCnt[], int MinCase[][2], int MaxCase[][2], int n1, int n2)
{
	for(int i = 0; i < m_nResidues; i++){
		GlobResImp[i] += ResidueImpact[i]; 
		double r = sqrt(ResidueImpact[i] / ResCnt[i]); 
		if(r > MaxImpact[i]){
			MaxImpact[i] = r;
			MaxCase[i][0] = n1;
			MaxCase[i][1] = n2;
			}
		if(r < MinImpact[i]){
			MinImpact[i] = r;
			MinCase[i][0] = n1;
			MinCase[i][1] = n2;
			}
		GlobResCnt[i] += ResCnt[i];
		}
}

void CCompare::SaveToFile(int n1, int n2, CpdbSite &s1, CpdbSite &s2){
	FILE *outf;
	char FName[256];
	int i;
	double rmsd = 0;

#ifdef _WIN32
	_mkdir(".\\pdbout");
	sprintf(FName, ".\\pdbout\\%d_%d_(%s_%s).pdb", ++n1, ++n2, s1.SiteName, s2.SiteName);
#else
	i = mkdir("./pdbout", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	sprintf(FName, "./pdbout/%d_%d_(%s_%s).pdb", ++n1, ++n2, s1.SiteName, s2.SiteName);
#endif
	outf = fopen(FName, "w");
	for(i = 0; i < m_nAtomPairs; i++){
		fprintf(outf, "ATOM  %5d  %-4s%3s A%4d    %8.3f%8.3f%8.3f %8.4f %d\n", i + 1, Pairs[i].Name, s1.Residues[Pairs[i].res_num].Name, s1.Residues[Pairs[i].res_num].num,
			Pairs[i].xyz1[0], Pairs[i].xyz1[1], Pairs[i].xyz1[2], sqrt(Pairs[i].diff), 0);
		rmsd += Pairs[i].diff;
		}

	rmsd = sqrt(rmsd / m_nAtomPairs);

	for(i = 0; i < m_nAtomPairs; i++){
		fprintf(outf, "ATOM  %5d  %-4s%3s B%4d    %8.3f%8.3f%8.3f %8.4f %d\n", i + 1, Pairs[i].Name, s1.Residues[Pairs[i].res_num].Name, s1.Residues[Pairs[i].res_num].num,
			Pairs[i].xyz2[0], Pairs[i].xyz2[1], Pairs[i].xyz2[2], rmsd, 0);
		}
	fclose(outf);
}

void CCompare::PrepareSymmetric(CpdbSite &s1){
	int n1, cur_res, cur_cnt, i, j, k;
	char *str[4];
	bool bConjug;

	m_nAltPairs = 0;
	for(i = 0; i < m_nAtomPairs; i++) Pairs[i].nAltPair = -1;
	cur_res = n1 = -1;
	for(i = 0; i < m_nAtomPairs; i++){
		 bConjug = false;
		if(Pairs[i].res_num != cur_res){
			cur_res = Pairs[i].res_num;
			n1 = -1;
			str[2] = NULL;
			cur_cnt = 0;

			if(strstr(s1.Residues[Pairs[i].res_num].Name, "ARG")){
				str[0] = "NH1";
				str[1] = "NH2";
				}
			else if(strstr(s1.Residues[Pairs[i].res_num].Name, "ASP")){
				str[0] = "OD1";
				str[1] = "OD2";
				}
			else if(strstr(s1.Residues[Pairs[i].res_num].Name, "GLU")){
				str[0] = "OE1";
				str[1] = "OE2";
				}
			else if(strstr(s1.Residues[Pairs[i].res_num].Name, "LEU")){
				str[0] = "CD1";
				str[1] = "CD2";
				}
			else if(strstr(s1.Residues[Pairs[i].res_num].Name, "VAL")){
				str[0] = "CG1";
				str[1] = "CG2";
				}
			else if(strstr(s1.Residues[Pairs[i].res_num].Name, "TYR")){
				str[0] = "CD1";
				str[1] = "CD2";

				str[2] = "CE1";
				str[3] = "CE2";
				bConjug = true;
				}
			else if(strstr(s1.Residues[Pairs[i].res_num].Name, "PHE")){
				str[0] = "CD1";
				str[1] = "CD2";

				str[2] = "CE1";
				str[3] = "CE2";
				bConjug = true;
				}
			else str[0] = NULL;
		}

		for(int l = 0; l < 4; l += 2)
		 if(str[l]){
			for(j = 0; j < 2; j++)
				if(strstr(Pairs[i].Name, str[l + j])){
					if(n1 >= 0){
						Pairs[n1].nAltPair = n1;
						AltPairs[m_nAltPairs][0] = Pairs[n1];
						for(k = 0; k < 3; k++) AltPairs[m_nAltPairs][0].xyz2[k] = Pairs[i].xyz2[k];

						Pairs[i].nAltPair = i;
						AltPairs[m_nAltPairs][1] = Pairs[i];
						for(k = 0; k < 3; k++) AltPairs[m_nAltPairs][1].xyz2[k] = Pairs[n1].xyz2[k];

						if(cur_cnt++) 
							AltPairs[m_nAltPairs - 1][0].nAltPair = -(AltPairs[m_nAltPairs - 1][0].nAltPair + 1);
						m_nAltPairs++;
						n1 = -1;
						str[l] = NULL;
						}
					else n1 = i;

					break;
					}
				}
		}
}


