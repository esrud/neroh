/*
 * Program: neroh
 *
 * Description: Ne Estimation by Runs Of Homozygosity
 *
 * Authors: Enrique Santiago and Carlos Köpke
 *
 * License: TBD
 */

#include <algorithm>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <omp.h>
#include <random>
#include <sstream>
#include <cstdlib>
#include <string>
#include <climits>

#include "lib/progress.hpp"

static constexpr int MAXLOCI  = 1000000;
static constexpr int MAXIND   = 1000;
static constexpr int MAXIND2  = 2000;
static constexpr int MAXCROMO = 1000;
static constexpr int MAXBINS  = 1000;


// ---------------------------------------------------------------------------
// Global data arrays
// ---------------------------------------------------------------------------
char indi[2][MAXIND][MAXLOCI]={'\0'}; // Diploid genotypes: 0=one allele, 1=the other, 3=unassigned
bool sexo[MAXIND];
char base[MAXLOCI]={'\0'}; // Reference bases for each locus
int rangocromo[MAXCROMO],rangocromo2[MAXCROMO],numcromo2,p[MAXLOCI];
int quitarango[MAXCROMO];
double frec[MAXLOCI]={0}, homo[MAXLOCI]={0};
int nfrec[MAXLOCI]={0};
double Het[MAXIND];
bool valid[MAXLOCI]={false};
double posicM[MAXLOCI];
long int posibp[MAXLOCI];
int *pj, *pi, *prefloc,*prefind,*pk,*pfinloc,*pfinind;
int eneind,eneindX,dobleinds,eneloc,enemachos,enehembras;
int n_sample,n_SNPs,n_threads,ff;
double acun,acupq,acuPp2;
double n,fs,fp,cenmor,cenmor2,mucM=0,errorcM;
double acuHetreal=0,Het_med=0,Het_total, cortel,Het_esp=0;
double tini,tpas;
double basebinsize,cMbase, cMtope;
bool flagleecMMbenmap={true}; // By default, assumes cM data comes from the .map file
double crec;
double dH5_total;

// ---------------------------------------------------------------------------
// Application parameters
// ---------------------------------------------------------------------------
struct AppParams {
    int diphap;
    int numThreads;
    int numSample;
    int numSNPs;
    double mu;
    double error;
    double hc;
    double lc;
    double cMMb;
    double binsize; // Bin width (-b option)
    double mindist; // Minimum distance between consecutive markers (-d option)
    double maxcM; // Maximum number of SNPs per cM (-c option)
    double MAF; // MAF cutoff threshold (-M option)
    bool HoHet; // Use Ho/Het method instead of heterozygosity
    bool onlyanalysis;
    bool unphase; // Unphase entries that are phased
    bool printToStdOut;
    bool flagX;
    bool flaghc;
    ProgressStatus progress;
};

struct AppParams params =
{
    .diphap = 2, // 0: unphased, 1: haploids, 2: phased, 4: two diploids
    .numThreads = 0,
    .numSample = 0,
    .numSNPs = 0,
    .mu=0,
    .error=0,
    .hc = 30,
    .lc = 1.0,
    .cMMb = 0,
    .binsize=0.125, // Bin width (-b option)
    .mindist = 0, // Minimum distance between consecutive markers (-d option)
    .maxcM=100000, // Maximum number of SNPs per cM (-c option)
    .MAF=0.05, // MAF cutoff threshold (-M option)
    .HoHet=false, // By default uses the heterozygosity method
    .onlyanalysis = false,
    .unphase = false,
    .printToStdOut = false,
    .flagX=false,
    .flaghc=false,
    .progress = ProgressStatus()
};

// ---------------------------------------------------------------------------
// Population information
// ---------------------------------------------------------------------------
struct PopulationInfo
{
    int numIndividuals;
    int numMachos;
    int numHembras;
    int numLoci;
    int fichLoci;
    int numcromo;
    double Mtot;
    double Mbtot;
};

// NEROH variables
// ---------------------------------------------------------------------------
// NEROH estimation constants and arrays
// ---------------------------------------------------------------------------
static constexpr int nlinmax     = 1000;
static constexpr int ngenmax     = 2000;
static constexpr int nbloquesmax = 100;

int neroh(const std::string& fichentra);

double SC;
double nbin[nlinmax];
double rohval[nlinmax];
double rohobs[nlinmax];
double xnbin[nlinmax]={0};
double ynbin[nlinmax]={0};
double xpeso[nlinmax]={0};
double xrohval[nlinmax]={0};
double xrohobs[nlinmax]={0};
double yrohobs[nlinmax]={0};
double rohprd[nlinmax];
double Neblock[nlinmax];
double Nerohs[nlinmax];
double edadrohs[nlinmax];
int segdesdehasta[nlinmax];
double Ne[ngenmax];
double sumNe[ngenmax];
double acusumNe[ngenmax]={0};
double acurohprd[nlinmax]={0};
double Necons, Nemed;
int nlin,maxnlin,gmax,gmax2,gmax3, nsegmentos,minedadrohs;
int muestrasalida=10;
int hp=2;
int conta,conta2,conta0=0;
double sample_size=0, sample_size_h=0;
bool flagporbloques=false, flagsolape=true;
double maxdouble=std::numeric_limits< double >::max();
unsigned long long int semilla=0;
int resolucion, topeposigen,ndes,gen, maxgen,des, ind1,ind2, ind3, posigen, posiblock, ancho1, ancho2;
int ndesini, indexmaxSC,indexminSC, maxsegmentos, posi1,posi2,posi3, contagen, topeposiinicio, topeposiinicio1,topeposiinicio4;
int tercio1, tercio2, tercio12, nhijos,topeposiinicio2,topeposiinicio3;
double efectomut,efectomutlateral, frecmut, frecrec,frecnorec, frecinversion, frecmutlateral, efecto, maxSC, minSC;
double efectomutsuave, SCmed, SCmed2, nsegmed,nsegmed2, SCbest,nsegbest, SCmedanterior;
double tope,topesalto,invtopesalto,topesalto2,invtopesalto2;
bool flag;
struct serie { // ++
    double SCval; // ++
    int nseg; // ++
    int segbl[nbloquesmax]; // ++
    double Nebl[nbloquesmax]; // ++
}; // ++
serie bichoP[1000];// ++ Parents
serie bichoH[1000];// ++ Offspring
serie bicho;// ++

double CalculaSC(serie& bicho, double* rohprd);

double Nc=1000,Gc=999999; // NOT USED
double Nb=100,Gb=75; // NOT USED
double Na=1000,Ga=50; // NOT USED
double heterozygosity; // Heterozygosity among segregating sites
double mu; // Mutation rate per Morgan
double err; // Error rate per Morgan
// double mucM; // Mutation rate per cM
// double errorcM; // tasa de error por cM
double d; // In Morgans
double dcM; // In cM
double H;
double genomesize;
double binwidth; // In Morgans
double binwidthcM; // In cM
double wl,suma;
int ncrom;
// int diphap;
int NREPETICIONES=14, REPE;
double dH_vals[5];
double FactorHoHetH;

std::random_device seed;  // used to obtain a seed for the random number engine
std::mt19937_64 genera(seed()); // mersenne_twister_engine 64bit (very good but very big)
std::uniform_real_distribution<> uniforme01(0.0, 1.0); // uses the result of the engine to generate uniform dist



void readFile(std::string fichPed, std::string fichMap, char (&population)[2][MAXIND][MAXLOCI], PopulationInfo (&popInfo))
{
    /*
     * Takes as input the file name and a pointer to the population
     * matrix and returns the number of individuals in the file
     */
    char base1[1], base2[1];
    int contaLociBase = 0;
    int conta = 0, posi = 0, posi2 = 0, longi = 0, ncopiMAFinf, ncopiMAFsup;
    int conta2, contapos, contaposant, aa, bb, i;
    int eneblocks, n1cM, j, jj, k, imin, imin1, imin2,quita;
    double sizechr, sizeblock, diff, minimo;
    int contalines=0, quitarestos,quitados,parmaxcMexacto,quitaeste;
    double posicMant=0;

    std::string line;
    std::string cromocod;
    std::string cromocodback = "laksjhbqne";

    std::ifstream entrada;

     // READING .map DATA:
    entrada.open(fichMap, std::ios::in); // Read loop for .map file
    if (!entrada.good())
    {
        std::cerr << "Could not open \"" << fichMap << "\". Does the file exist?" << std::endl;
        exit(EXIT_FAILURE);
    }

    while (std::getline(entrada,line)){
        longi=int(line.length());
        if (longi<8){
            std::cerr << "Line too short in map file" << std::endl;
            exit(EXIT_FAILURE);
        }

        posi=int(line.find_first_of(" \t",0));
        if (posi <= 0) {
            std::cerr << "Empty line in map file" << std::endl;
            exit(EXIT_FAILURE);
        }
        cromocod = line.substr(0, posi);
        if (cromocod != cromocodback){
            cromocodback = cromocod;
            rangocromo[popInfo.numcromo] = contalines;
            ++ popInfo.numcromo;
        }

        ++posi; // SNP name field (not read)
        posi2=posi;
        posi=int(line.find_first_of(" \t",posi2));
        if (posi <= 0) {
            std::cerr << "Error in map file (1)" << std::endl;
            exit(EXIT_FAILURE);
        }

        ++posi; // cM position
        posi2=posi;
        posi=int(line.find_first_of(" \t",posi2));
        if (posi <= 0) {
            std::cerr << "Error in map file (2)" << std::endl;
            exit(EXIT_FAILURE);
        }
        posicM[contalines]=std::stod(line.substr(posi2, posi-posi2));
        if (params.flagX){ // For X chromosome, compress distances by a factor of 2/3
          posicM[contalines]*=0.6666666;
        }
        ++posi;

        if (posi > longi) {    // bp position
            std::cerr << "Error in map file (3)" << std::endl;
            exit(EXIT_FAILURE);
        }
        posibp[contalines]=std::stoi(line.substr(posi, longi-posi));

        if (!flagleecMMbenmap){
            posicM[contalines] = params.cMMb * posibp[contalines] / 1000000;
            if (params.flagX){ // For X chromosome, compress distances by a factor of 2/3
              posicM[contalines]*=0.6666666;
            }
        }
        ++contalines;
    }

    rangocromo[popInfo.numcromo] = contalines;
    popInfo.Mbtot=0;
    for (conta=0; conta<popInfo.numcromo; ++conta){
        popInfo.Mbtot += posibp[rangocromo[conta+1]-1]-posibp[rangocromo[conta]];
    }
    popInfo.Mbtot /= 1000000; // Convert to megabases (from here on we use Mb)
    entrada.close();

  for (i=0;i<popInfo.numLoci;++i){
    frec[i]=0;// Counter for non-reference allele copies
    nfrec[i]=0;
    homo[i]=0;
  }
  // READING .ped DATA:
    entrada.open(fichPed, std::ios::in); // Read loop for .ped file
    if (!entrada.good())
    {
        std::cerr << "Could not open \"" << fichPed << "\". Does the file exist?" << std::endl;
        exit(EXIT_FAILURE);
    }
    while (std::getline(entrada,line)){
        longi=int(line.length());
        if (longi<12){
            std::cerr << "Line too short in ped file" << std::endl;
            exit(EXIT_FAILURE);
        }
        conta=0;
        posi=0;
        while ((posi < longi) && (conta < 6)){
            posi2=posi;
            posi=int(line.find_first_of(" \t",posi2));
            if (posi < 0) {
                std::cerr << "Error in file .ped" << std::endl;
                exit(EXIT_FAILURE);
            }
            ++posi;
            ++conta;
            if (conta==4){// Sexo
                base1[0]=line.at(posi); // Reusing base1 variable; this is actually the sex field
                if (base1[0]=='2'){
                    sexo[popInfo.numIndividuals]=true; // true = female
                }
            }
        }

        if (conta==6){
            contapos = 0;
            while(posi < longi) { // asigna genot.
                base1[0]=line.at(posi);
                posi2=posi;
                posi=int(line.find_first_of(" \t",posi2));
                if (posi < 0) {break;}
                ++posi;
                base2[0]=line.at(posi);
                k=0;
                if (params.unphase){
                  if (uniforme01(genera) <0.5){
                    std::swap(base1[0],base2[0]);
                  }
                }
                if (base1[0]!='0'){ // primer cromosoma
                    if (base[contapos] == '\0'){
                        base[contapos] = base1[0];
                    }
                    if (base1[0]==base[contapos]){
                        population[0][popInfo.numIndividuals][contapos] = 0; // Reference allele
                    }
                    else{
                        population[0][popInfo.numIndividuals][contapos] = 1; // Alternative allele
                        ++k;
                        ++frec[contapos];
                    }
                    ++nfrec[contapos];
                }
                else{
                    population[0][popInfo.numIndividuals][contapos] = 3; // Unassigned (missing)
                }
                if (base2[0]!='0'){ // The other chromosome
                    if (base[contapos] == '\0'){
                        base[contapos] = base2[0];
                    }
                    if (base2[0]==base[contapos]){
                        population[1][popInfo.numIndividuals][contapos] = 0; // Reference allele
                    }
                    else{
                        population[1][popInfo.numIndividuals][contapos] = 1; // Alternative allele
                        if ((!params.flagX) || (sexo[popInfo.numIndividuals])){
                            ++k;
                            ++frec[contapos];
                        }
                    }
                    if ((!params.flagX) || (sexo[popInfo.numIndividuals])){
                        ++nfrec[contapos];
                    }
                }
                else{
                    population[1][popInfo.numIndividuals][contapos] = 3; // Unassigned (missing)
                }
                homo[contapos]+=(k==2); // Accumulate homozygotes
                posi2=posi;
                posi=int(line.find_first_of(" \t",posi2));
                if (posi < 0) {posi=longi;}
                ++posi;
                ++contapos;
                if (contapos >= MAXLOCI) {
                    std::cerr<<"Reached max number of loci (" << MAXLOCI << ")" << std::endl;
                    break;
                }
            }

            if (popInfo.numIndividuals == 0){
                contaLociBase = contapos;
            }

            if (contapos != contaLociBase){
                std::cerr << "Some genomes in the sample have different sizes" << std::endl;
                exit(EXIT_FAILURE);
            }

            if (params.flagX){
                if (sexo[popInfo.numIndividuals]){
                  ++popInfo.numHembras;
                }
                else{
                  ++popInfo.numMachos;
                }
            }

            popInfo.numIndividuals++;
            if (popInfo.numIndividuals > MAXIND){
                std::cerr << "Reached limit of sample size (" << MAXIND <<")" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    entrada.close();
    if (params.flagX && (params.diphap==2)){
        if (((popInfo.numHembras *2) + popInfo.numMachos) > MAXIND){
            std::cerr << "The number of genomes is too large " << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    popInfo.numLoci=contapos;
    if (popInfo.numLoci != contalines){
        std::cerr << "Different number of loci in ped and map files" << std::endl;
        exit(EXIT_FAILURE);
    }

    // ELIMINACION DE NO-VALIDOS:
    // Remove non-segregating, incomplete, and below-MAF loci:
    if (params.flagX){
        dobleinds=popInfo.numHembras*2+popInfo.numMachos;
    }
    else{
        dobleinds=popInfo.numIndividuals*2;
    }
    ncopiMAFinf=int(dobleinds*params.MAF);
    if (ncopiMAFinf<1){ncopiMAFinf=1;}
    ncopiMAFsup=dobleinds-ncopiMAFinf;
    popInfo.fichLoci = contapos;
    conta=0;
    for (i=0;i<contapos;++i){
        valid[i]=false;
        if (nfrec[i]==dobleinds){ // Only keep markers with complete information
            if ((frec[i]>=ncopiMAFinf) && (frec[i]<=ncopiMAFsup)){
                p[conta]=i;
                valid[i]=true;
                ++conta;
            }
        }
    }

    contalines=conta;
    // Readjust chromosome pointers to now point through the p[] indirection array:
    for (j=0; j<popInfo.numcromo;++j){ // Chromosome loop
        quitarango[j]=0;
        for (k=rangocromo[j]; k<rangocromo[j+1]; ++k){
            if (!valid[k]){++quitarango[j];}
        }
    }
    for (j=0; j<popInfo.numcromo;++j){ // Chromosome loop
        for (k=j+1; k<=popInfo.numcromo; k++){
            rangocromo[k]-=quitarango[j];
        }
    }

    for (j=0; j<popInfo.numcromo; ++j){
        j=j;
        if ((rangocromo[j+1] - rangocromo[j])<100){
           std::cerr << "Too few genotyping data for chromosome " <<j+1<< std::endl;
           exit(EXIT_FAILURE);
        }
    }

    conta=0; // Thin markers if density exceeds params.maxcM per cM
    for (j=0; j<popInfo.numcromo;++j){ // Chromosome loop
        sizechr= posicM[p[rangocromo[j+1]-1]]-posicM[p[rangocromo[j]]];
        eneblocks = (int)(sizechr/1.0)+1; // Integer number of blocks of ~1 cM
        sizeblock=sizechr/eneblocks; // Exact block size in cM (close to 1 cM)
        conta=rangocromo[j];
        conta2=conta;
        for (jj=0;jj<eneblocks;++jj){
            n1cM=0;
            conta=conta2;
            while ((posicM[p[conta2]]<(posicM[p[conta]]+sizeblock)) && (conta2<rangocromo[j+1])){ // Count markers in this block
                ++n1cM;
                ++conta2;
            }
            if (n1cM>params.maxcM){ // If marker count exceeds max, remove the excess
                quitarestos =n1cM-params.maxcM;
                for (k=0;k<(quitarestos);++k){ // These are the excess ones
                    minimo=99999;
                    for (i=conta+1;i<conta2-1;++i){ // Cual quito?
                        diff=floor((posicM[p[i+1]]-posicM[p[i-1]])*1000000)/1000000;
                        if (minimo>diff){
                            minimo=diff;
                            imin=i;
                        }
                    }
                    for (i=imin;i<conta2-1;++i){
                        p[i]=p[i+1];  // Remove the pointer to this site
                    }
                    --conta2;
                }
                contalines-=quitarestos;
                for (i=conta2;i<contalines;++i){// Readjust pointers in the following ranges
                    p[i]=p[i+quitarestos];  // Remove the pointer to this site
                }
                for (i=j+1;i<=popInfo.numcromo;++i){
                    rangocromo[i]-=quitarestos;
                }
            }
        }
    }

    i=0;
    if (params.mindist>0){ // Thin markers for inter-adjacent distance less than params.mindist
        for (j=0; j<popInfo.numcromo;++j){ // Chromosome loop
            posicMant=posicM[p[rangocromo[j]]];
            jj=rangocromo[j]+1;
            quitados=0;
            while (jj<rangocromo[j+1]){
                if (((posicM[p[jj]]-posicMant)) < params.mindist){
                    for (i=jj;i<rangocromo[j+1]-1;++i){
                        p[i]=p[i+1];  // Remove the pointer to this site
                    }
                    ++quitados;
                    --rangocromo[j+1];
                }
                else{
                    posicMant=posicM[p[jj]];
                    ++jj;
                }
            }
            contalines-=quitados;// Readjust pointers in the following chromosomes
            for (i=rangocromo[j+1];i<contalines;++i){
                p[i]=p[i+quitados];  // Remove the pointer to this site
            }
            for (i=j+2;i<=popInfo.numcromo;++i){
                rangocromo[i]-=quitados;
            }
        }
    }

    rangocromo[popInfo.numcromo] = contalines;
    popInfo.numLoci = contalines;
    entrada.close();

    popInfo.Mtot=0;
     for (conta=0; conta<popInfo.numcromo; ++conta){
        popInfo.Mtot += posicM[p[rangocromo[conta+1]-1]]-posicM[p[rangocromo[conta]]];
    }
    popInfo.Mtot /= 100; // In Morgans
}

void printHelp(char * appName)
{
    fprintf(stderr,
        "NEROH (v1.0 - April 2022)\n"
        "Authors: Enrique Santiago - Carlos Köpke\n"
        "\n"
        "Usage: %s [OPTIONS] <file_name>\n"
        "       where file_name has no extension .ped nor .map\n"
        "Examples: %s -o 2 -M 0.05 file\n"
        "          %s -l 0.5 -u 20 -c 50 file\n"
        "\n"
        "  OPTIONS:\n"
        "\n"
        "    -h     Print out this help\n"
        "    -t     Number of threads (default: %d)\n"
        "    -o     Diploids unphased:0, haploids:1, diploid phased:2 (default: 2)\n"
        "    -X     Linked to chromosome X\n"
        "    -m     Mutation rate per Morgan\n"
        "    -e     Genotyping error rate per site (default: 0)\n"
        "    -l     Lower bound distance to be considered in cM (default: 1.0)\n"
        "    -u     Upper bound distace to be considered in cM (default: 30)\n"
        "    -r     If specified, constant rec rate in cM/Mb across the genome\n"
        "    -b     Bin size in cM (for estimation of Ne using rohs - default 0.125)\n"
        "    -d     Minimum distance between adyacent markers (default 0.0)\n"
        "    -c     Cutoff for the maximum number of markers per cM (default all)\n"
        "    -M     MAF threshold (default: 0.05). Forced to 0 if mut. rate is given. \n"
//        "    -H     Usa frecuencia observada de sitios HoHet en vez de Heterocigosidad\n"
        "    -a     Performs only heterozigosity and marker density analysis\n"
//        "    -f     Unphases phased inputs.\n"
        "    -p     Print analysis to stdout. If not specified a file will be created\n"
        "\n",
        appName,
        appName,
        appName,
        omp_get_max_threads()
    );
}

int main(int argc, char *argv[]) {
  int i, ii, j, jj, j2, j3, conta, conta2, crom_i, crom_ii, ind_i,
      ind_ii, eneind2;
  int containdX, containdX2, contalocX, nsegments, eneblocks;
  int n5cM, nceros, desde_het, hasta_het;
  double sizechr, sizeblock, d5cM, H5cMA, H5cMG, H5cMAG, dH5cM;
  double f1fo, f1fe, sitedist, dh5_total;
  bool superadolimite;
  int binmax;
  int Matbinmax[MAXIND2]={0};
  int nbinroh[MAXBINS] = {0};
  double cMbinroh[MAXBINS] = {0};
  double matdH[MAXBINS];
  int tamavals,eneunits;
  double totHoHo, totHoHet, HoHet;
  bool flagmu;

  // Variables for OpenMP pragmas:
  int xj, xj3, xindi_i, xindi_ii, xcrom_i, xcrom_ii, xii, xconta, xconta2;
  bool xiniciocrom, xroh;
  double xcenmor;
  bool Hoant,Honow;
  long int contaHoHo[MAXIND2]={0}, contaHoHet[MAXIND2]={0};

  std::random_device rd;
  std::mt19937 g(rd());
  g.seed(rd());
  std::uniform_real_distribution<> uniforme01(0.0, 1.0);

  flagmu=false;
  for (;;) {
    switch (getopt(argc, argv, "ho:Xm:e:l:u:r:b:d:c:M:Ht:afp")) {
    case '?':
    case 'h':
    default:
      printHelp(argv[0]);
      return -1;
    case 'o':
      params.diphap = std::atoi(optarg);
      if ((params.diphap < 0) || (params.diphap > 2)) {
        std::cerr << "Invalid code haploid/diploid" << std::endl;
        return -1;
      }
      continue;
    case 'X':
      params.flagX = true;
      continue;
    case 'm':
      params.mu = std::atof(optarg);
      flagmu=true;
      if ((params.mu < 0)) {
        std::cerr << "Invalid mutation rate per Morgan" << std::endl;
        return -1;
      }
      continue;
    case 'e':
      params.error = std::atof(optarg);
      if ((params.error < 0) || (params.error > .1)) {
        std::cerr << "Invalid error rate" << std::endl;
        return -1;
      }
      continue;
    case 'l':
      params.lc = std::atof(optarg);
      if ((params.lc < 0) || (params.lc > 100)) {
        std::cerr << "Invalid minimum distance in centiMorgans" << std::endl;
        return -1;
      }
      continue;
    case 'u':
      params.hc = std::atof(optarg);
      params.flaghc=true;
      if ((params.hc < 0) || (params.hc > 100)) {
        std::cerr << "Invalid maximum distance in centiMorgans" << std::endl;
        return -1;
      }
      continue;
    case 'r':
      params.cMMb = std::atof(optarg);
      if (params.cMMb <= 0) {
        std::cerr << "Invalid ratio cM/Mb" << std::endl;
        return -1;
      }
      continue;
    case 'b':
      params.binsize = std::atof(optarg);
      if (params.binsize <= 0) {
        std::cerr << "Invalid bin size" << std::endl;
        return -1;
      }
      continue;
    case 'd':
      params.mindist = std::atof(optarg);
      if (params.mindist < 0) {
        std::cerr << "Invalid distance" << std::endl;
        return -1;
      }
      continue;
    case 'c':
      params.maxcM = std::atoi(optarg);
      if (params.maxcM < 2) {
        std::cerr << "Invalid number of markers/cM" << std::endl;
        return -1;
      }
      continue;
    case 'M':
      params.MAF = std::atof(optarg);
      if (params.MAF < 0) {
        std::cerr << "Invalid MAF cut-off" << std::endl;
        return -1;
      }
      continue;
    case 'H':
      params.HoHet = true;
      continue;
    case 't':
      params.numThreads = std::atoi(optarg);
      continue;
    case 'a':
      params.onlyanalysis = true;
      continue;
    case 'f':
      params.unphase = true;
      continue;
    case 'p':
      params.printToStdOut = true;
      continue;
    case -1:
      break;
    }
    break;
  }
  if (params.cMMb > 0) {
    flagleecMMbenmap = false;
  }
  if (flagmu) { // Si se da tasa de mutacion, entonces MAF=0
    params.MAF = 0;
  }

  std::string fich = "";
  std::string fichped = "";
  std::string fichmap = "";

  if (optind < argc) {
    fich = argv[optind];
  }

  if (fich == "") {
    std::cerr << "Missing data file name" << std::endl;
    return -1;
  }

  fichped = fich + ".ped";
  fichmap = fich + ".map";
  std::string fichProgress = fich + "_NEROH_progress.tmp";
  params.progress.InitTotalTasks(2, fichProgress.c_str());
  params.progress.SetCurrentTask(0, "Analysis of input data");
  params.progress.InitCurrentTask(7);
  params.progress.SetTaskProgress(0);
    // params.progress.PrintProgress();

  n_sample = params.numSample;
  n_SNPs = params.numSNPs;
  n_threads = params.numThreads;

  std::ifstream entrada;

  if (n_threads > 0) {
    omp_set_num_threads(n_threads);
  }

  std::cout << "  A progress report is stored at " << fichProgress << ".\n";
  std::cout << "  Check it by issuing 'cat " << fichProgress << "'.\n";
  std::cout << "  Reading files " << fich <<"."<< std::endl;

  tini = omp_get_wtime();
  // Read simulation/input data:
  // ---------------------------------------------------------------------------
// Population information
// ---------------------------------------------------------------------------
struct PopulationInfo popInfo = {
      .numIndividuals = 0,
      .numMachos=0,
      .numHembras=0,
      .numLoci = 0,
      .fichLoci = 0,
      .numcromo = 0,
      .Mtot = 0,
      .Mbtot = 0,
  };

  params.progress.SetTaskProgress(1);

  readFile(fichped, fichmap, indi, popInfo); // <<< Read input files

  double tProcessFile = (omp_get_wtime() - tini);
  std::cout << "  Reading .map and .ped files took " << std::fixed << std::setprecision(2) << tProcessFile << " sec." << std::endl;
  std::cout << "  Processing files ... " << std::endl;

  params.progress.SetTaskProgress(3);

  eneloc = popInfo.numLoci;

  eneind = popInfo.numIndividuals;
  enemachos = popInfo.numMachos;
  enehembras = popInfo.numHembras;

  // Compute non-reference allele frequency and homozygosity for all
  // variable loci in the sample
  params.progress.SetTaskProgress(4);
    // params.progress.PrintProgress();
  for (j = 0; j < popInfo.numLoci; ++j) {
    frec[p[j]] /= dobleinds;
    if (params.flagX){
      homo[p[j]] /= enehembras;
    }
    else{
      homo[p[j]] /= eneind;
    }
  }

  // Compute heterozygosity in 5 cM segments

  std::string fichsal3 = fich + "_NeROH_SITES.txt";
  std::ofstream outputFile3;
  std::stringstream salida3;

  outputFile3.open(fichsal3);
  salida3 << "Chromosome"
          << "\t"
          << "Block_No."
          << "\t"
          << "First_Het_Marker_No."
          << "\t"
          << "Last_Het_Marker_No."
          << "\t"
          << "First_Het_Marker_Loc(cM)"
          << "\t"
          << "Last_Het_Marker_Loc(cM)"
          << "\t"
          << "Number_of_Het_Markers"
          << "\t"
          << "Average_Distance(cM)"
          << "\t"
          << "Average_Heterozygosity"
          << "\t"
          << "Average_Dist/Het"
          << "\n";

  acuPp2 = 0;
  acupq = 0;
  nsegments = 0;
  Het_total = 0;
  sitedist = 0;
  dH5_total = 0;
  params.progress.SetTaskProgress(5);
    // params.progress.PrintProgress();
  acuHetreal = 0;
  nceros=0;
  numcromo2=0;
  bool estoyencromo;
  for (j = 0; j < popInfo.numcromo; ++j) { // Chromosome loop
    sizechr = posicM[p[rangocromo[j + 1] - 1]] - posicM[p[rangocromo[j]]];
    // Only chromosomes larger than 40cM are considered
    if (sizechr<40){
      continue;
    }
    eneblocks =(int)(sizechr / 5.0) + 1; // Integer number of blocks of ~5 cM
    sizeblock = sizechr / eneblocks; // Exact block size in cM (close to 5 cM)
    conta = rangocromo[j];
    conta2 = conta;
    estoyencromo=false;
    for (jj = 0; jj < eneblocks; ++jj) {
      n5cM = 0;
      d5cM = 0;
      H5cMA = 0;
      dH5cM = 0;
      H5cMG = 0;
      H5cMAG = 0;
      desde_het = -1;
      hasta_het = -1;
      // while ((conta2 < rangocromo[j + 1]) &&(posicM[p[conta2]]<((jj+1)*sizeblock+posicM[p[rangocromo[j]]]+ 0.00001))){
      while ((conta2 < rangocromo[j + 1]) && ((posicM[p[conta2]] - posicM[p[conta]]) <= sizeblock)) {
          if (desde_het < 0) {
            desde_het = conta2;
          }
          hasta_het = conta2;
          f1fo = (frec[p[conta2]]-homo[p[conta2]]); // Observed
          f1fe = (frec[p[conta2]] * (1.0 - frec[p[conta2]])); // Expected
          acuPp2 += (homo[p[conta2]] - frec[p[conta2]] * frec[p[conta2]]);
          acupq += f1fe;
          f1fe += f1fe;
          f1fo += f1fo;

          if (f1fo>0){ //
            // H5cMAG += f1fo; // Expected Het for arithmetic mean
            H5cMAG += log10(f1fo); // Expected Het for geometric mean
            ++n5cM;
          }

          ++conta2;
      }
      if ((n5cM > 0)) {
        if (!estoyencromo){
          estoyencromo=true;
          rangocromo2[numcromo2]=conta;
          numcromo2++;
        }
        H5cMAG /= n5cM;

        H5cMAG = pow(10, H5cMAG); // Only for geometric mean; disable if using arithmetic mean

        if (n5cM>1){
          d5cM = (posicM[p[hasta_het]] - posicM[p[desde_het]]) / (n5cM - 1);
        }
        else{
          d5cM=sizeblock;
        }
        dH5cM = d5cM / H5cMAG;
        for (i = 0; i < nsegments; ++i) { // Insert dH5cM in sorted order
          if (dH5cM < matdH[i]) {
            break;
          }
        }
        for (ii = nsegments; ii > i; --ii) {
          matdH[ii] = matdH[ii - 1];
        }
        matdH[i] = dH5cM;

        salida3 << std::fixed << std::setprecision(0);
        salida3 << j+1 << "\t" << jj+1 << "\t"  ;
        salida3 << p[conta] << "\t" << p[conta2 - 1]<< "\t" ;
        salida3 << std::fixed << std::setprecision(6);
        salida3<< posicM[p[conta]] << "\t" << posicM[p[conta2 - 1]] << "\t";
        salida3 << std::fixed << std::setprecision(0);
        salida3<< n5cM << "\t";
        salida3 << std::fixed << std::setprecision(6);
        salida3<< d5cM << "\t" << H5cMAG<< "\t" << dH5cM << "\n";
        // matdH[i] = dH5cM;
        Het_total += H5cMAG;
        sitedist += d5cM;
        dH5_total += dH5cM;
        nsegments++;
        if ((conta2 == rangocromo[j + 1])){
          rangocromo2[numcromo2]=conta2;
          estoyencromo=false;
        }
      }
      else{
        if (estoyencromo){
          rangocromo2[numcromo2]=conta2;
          estoyencromo=false;
        }

        salida3 << std::fixed << std::setprecision(0);
        salida3 << j+1 << "\t" << jj+1 << "\t"  ;
        salida3 << "(5cM_without_het_markers)"<<"\n";
      }
      conta = conta2;
    }
  }
  if (estoyencromo){
    rangocromo2[numcromo2]=conta;
  }
  outputFile3 << salida3.str();
  outputFile3.close();

  int numcromo1 = popInfo.numcromo;
  // Update new chromosome boundaries or breakpoints
  popInfo.numcromo = numcromo2;
  for (i=0;i<=numcromo2;++i){
    rangocromo[i]=rangocromo2[i];
  }

  if (params.onlyanalysis){
    std::cout << "  Done. " << std::endl;
    std::remove(fichProgress.c_str());
    return 0;
  }

  tamavals = int(nsegments / 5);
  for (i = 0; i < 5; ++i) {
    dH_vals[i] = 0;
    int basa1 = i * tamavals;
    int basa2 = basa1 + tamavals;
    for (ii = basa1; ii < basa2; ++ii) {
      dH_vals[i] += matdH[ii];
    }
    dH_vals[i] /= tamavals;
  }
  Het_total /= nsegments;
  sitedist /= nsegments;
  dH5_total /= nsegments;
  fs = acuPp2 / acupq;
  Het_esp = 2 * acupq / eneloc;

  cortel = 0.5 / Het_total * sitedist; // Correction to add to fragments
                                       // telom.
  mucM = params.mu / 100;       // Mutation rate per cM
  errorcM = params.error / sitedist; // tasa de error de genotipado por cM

  cMbase = params.lc - params.binsize / 2;
  if (cMbase < 0) {
    cMbase = 0;
  } // Lower class boundary

  cMtope = (floor((params.hc - cMbase) / params.binsize)) * params.binsize +
           params.binsize + cMbase;

  params.progress.SetTaskProgress(6);
    // params.progress.PrintProgress();


  if ((!params.flagX) && (params.diphap == 0)) { // DIPLOIDS UNPHASED
    binmax = 0; // DIPLOIDS UNPHASED
    eneunits = eneind;
    #pragma omp parallel for private(Hoant, Honow, xj, xj3, xconta, xconta2, xroh, xiniciocrom, xcenmor)
    for (i = 0; i < eneind; ++i) {
      Hoant=false;
      for (xj = 0; xj < popInfo.numcromo; ++xj) {
        xconta = rangocromo[xj];
        xconta2 = xconta;
        xroh = false;
        xiniciocrom = true;
        while (xconta2 < rangocromo[xj + 1]) {
          Honow=(indi[0][i][p[xconta2]] == indi[1][i][p[xconta2]]);
          if (Hoant){
            if (Honow){++contaHoHo[i];}
            else{++contaHoHet[i];}
          }
          Hoant=Honow;
          if (Honow) {
            if (!xroh) {
              xroh = true;
              xconta = xconta2;
              if (xconta > rangocromo[xj]) {
                --xconta;
              } // New criterion: distance between two het sites
            }
          }
          else {
            if (xroh) {
              xcenmor = fabs(posicM[p[xconta2]] - posicM[p[xconta]]);
              if (xiniciocrom) {
                xcenmor += cortel;
              } // Telomere correction at start
              if ((xcenmor >= cMbase) && (xcenmor <= cMtope)) {
                xj3 = floor(float((xcenmor - cMbase) / params.binsize));
                if (xj3 < MAXBINS) {
                  xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                               params.binsize +
                           params.binsize / 2.0 +
                           cMbase; // Use central bin mark
                  if (xj3>Matbinmax[i]){
                      Matbinmax[i] = xj3;
                  }
              // All these updates are atomic:
                    #pragma omp atomic
                    ++nbinroh[xj3];
                    #pragma omp atomic
                    cMbinroh[xj3] += xcenmor;
                }
              }
              xroh = false;
            }
            xiniciocrom = false;
          }
          ++xconta2;
        }
        if (xroh) {
          xcenmor = fabs(posicM[p[xconta2 - 1]] - posicM[p[xconta - 1]]);
          xcenmor += cortel + sitedist; // correcc. Telom. Final
          if (xcenmor >= cMbase && xcenmor <= cMtope) {
            xj3 = floor(float((xcenmor - cMbase) / params.binsize));
            if (xj3 < MAXBINS) {
              xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                           params.binsize +
                       params.binsize / 2.0 +
                       cMbase; // Use central bin mark
              if (xj3>Matbinmax[i]){
                  Matbinmax[i] = xj3;
              }
          // All these updates are atomic:
                #pragma omp atomic
                ++nbinroh[xj3];
                #pragma omp atomic
                cMbinroh[xj3] += xcenmor;
            }
          }
        }
      }
    }
    totHoHo=0;
    totHoHet=0;
    for (i = 0; i < eneind; ++i) {
        totHoHo += contaHoHo[i];
        totHoHet += contaHoHet[i];
        if (Matbinmax[i]>binmax){
          binmax=Matbinmax[i];
        }
    }
    HoHet = totHoHet/(totHoHet+totHoHo);
  }


  else if ((!params.flagX) && (params.diphap == 2)) { // DIPLOIDES FASE conocida
    binmax = 0; // DIPLOIDES FASE conocida
    eneind2 = eneind * 2;
    eneunits = (eneind2*(eneind2-1))/2;
    for (i = 0; i < eneind2 - 1; ++i) {
      #pragma omp parallel for private(Hoant, Honow, xj, xj3, xindi_i, xindi_ii, xcrom_i, xcrom_ii, xii, xconta, xconta2, xroh, xiniciocrom, xcenmor)
      for (xii = i + 1; xii < eneind2; ++xii) {
      xindi_i = i / 2;
      xcrom_i = i - xindi_i * 2;
        xindi_ii = xii / 2;
        xcrom_ii = xii - xindi_ii * 2;
        Hoant=false;
        for (xj = 0; xj < popInfo.numcromo; ++xj) {
          xconta = rangocromo[xj];
          xconta2 = xconta;
          xroh = false;
          xiniciocrom = true;
          while (xconta2 < rangocromo[xj + 1]) {
            Honow=(indi[xcrom_i][xindi_i][p[xconta2]] == indi[xcrom_ii][xindi_ii][p[xconta2]]);
            if (Hoant){
              if (Honow){++contaHoHo[xii];}
              else{++contaHoHet[xii];}
            }
            Hoant=Honow;
            if (Honow) {
              if (!xroh) {
                xroh = true;
                if (xconta > rangocromo[xj]) {
                  --xconta;
                } // New criterion: distance between two het sites
                xconta = xconta2;
              }
            } else {
              if (xroh) {
                xcenmor = fabs(posicM[p[xconta2]] - posicM[p[xconta]]);
                if (xiniciocrom) {
                  xcenmor += cortel;
                } // Telomere correction at start
                if (xcenmor >= cMbase && xcenmor <= cMtope) {
                  xj3 = floor(float((xcenmor - cMbase) / params.binsize));
                  if (xj3 < MAXBINS) {
                    xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                                 params.binsize +
                             params.binsize / 2.0 +
                             cMbase; // Use central bin mark
                  if (xj3>Matbinmax[xii]){
                      Matbinmax[xii] = xj3;
                  }
              // All these updates are atomic:
                    #pragma omp atomic
                    ++nbinroh[xj3];
                    #pragma omp atomic
                    cMbinroh[xj3] += xcenmor;
                  }
                }
                xroh = false;
              }
              xiniciocrom = false;
            }
            ++xconta2;
          }
          if (xroh) {
            xcenmor = fabs(posicM[p[xconta2 - 1]] - posicM[p[xconta - 1]]);
            xcenmor += cortel + sitedist; // correcc. Telom. Final
            if (xcenmor >= cMbase && xcenmor <= cMtope) {
              xj3 = floor(float((xcenmor - cMbase) / params.binsize));
              if (xj3 < MAXBINS) {
                xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                             params.binsize +
                         params.binsize / 2.0 +
                         cMbase; // Use central bin mark
                if (xj3>Matbinmax[xii]){
                    Matbinmax[xii] = xj3;
                }
         // All these updates are atomic:
                #pragma omp atomic
                ++nbinroh[xj3];
                #pragma omp atomic
                cMbinroh[xj3] += xcenmor;
              }
            }
          }
        }
      }
    }
    totHoHo=0;
    totHoHet=0;
    for (i = 0; i < eneind2 - 1; ++i) {
        totHoHo += contaHoHo[i];
        totHoHet += contaHoHet[i];
        if (Matbinmax[i]>binmax){
          binmax=Matbinmax[i];
        }
    }
    HoHet = totHoHet/(totHoHet+totHoHo);
  }


  else if ((!params.flagX) && (params.diphap == 1)) { // HAPLOIDES
    binmax = 0;  // HAPLOIDES
    eneunits = eneind*(eneind-1)/2;
    for (i = 0; i < eneind - 1; ++i) {
    #pragma omp parallel for private(Hoant, Honow, xj, xii, xj3, xconta, xconta2, xroh, xiniciocrom, xcenmor)
      for (xii = i + 1; xii < eneind; ++xii) {
        Hoant=false;
        for (xj = 0; xj < popInfo.numcromo; ++xj) {
          xconta = rangocromo[xj];
          xconta2 = xconta;
          xroh = false;
          xiniciocrom = true;
          while (xconta2 < rangocromo[xj + 1]) {
            Honow=(indi[0][i][p[xconta2]] == indi[0][xii][p[xconta2]]);
            if (Hoant){
              if (Honow){++contaHoHo[xii];}
              else{++contaHoHet[xii];}
            }
            Hoant=Honow;
            if (Honow) {
              if (!xroh) {
                xroh = true;
                xconta = xconta2;
                if (xconta > rangocromo[xj]) {
                  --xconta;
                } // New criterion: distance between two het sites
              }
            } else {
              if (xroh) {
                xcenmor = fabs(posicM[p[xconta2]] - posicM[p[xconta]]);
                if (xiniciocrom) {
                  xcenmor += cortel;
                } // Telomere correction at start
                if (xcenmor >= cMbase && xcenmor <= cMtope) {
                  xj3 = floor(float((xcenmor - cMbase) / params.binsize));
                  if (xj3 < MAXBINS) {
                    xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                                 params.binsize +
                             params.binsize / 2.0 +
                             cMbase; // Use central bin mark
                    if (xj3>Matbinmax[xii]){
                        Matbinmax[xii] = xj3;
                    }
                // All these updates are atomic:
                      #pragma omp atomic
                      ++nbinroh[xj3];
                      #pragma omp atomic
                      cMbinroh[xj3] += xcenmor;
                  }
                }
                xroh = false;
              }
              xiniciocrom = false;
            }
            ++xconta2;
          }
          if (xroh) {
            xcenmor = fabs(posicM[p[xconta2 - 1]] - posicM[p[xconta - 1]]);
            xcenmor += cortel + sitedist; // correcc. Telom. Final
            if (xcenmor >= cMbase && xcenmor <= cMtope) {
              xj3 = floor(float((xcenmor - cMbase) / params.binsize));
              if (xj3 < MAXBINS) {
                xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                             params.binsize +
                         params.binsize / 2.0 +
                         cMbase; // Use central bin mark
              if (xj3>Matbinmax[xii]){
                  Matbinmax[xii] = xj3;
              }
          // All these updates are atomic:
                #pragma omp atomic
                ++nbinroh[xj3];
                #pragma omp atomic
                cMbinroh[xj3] += xcenmor;
              }
            }
          }
        }
      }
    }
    totHoHo=0;
    totHoHet=0;
    for (i = 0; i < eneind - 1; ++i) {
        totHoHo += contaHoHo[i];
        totHoHet += contaHoHet[i];
        if (Matbinmax[i]>binmax){
          binmax=Matbinmax[i];
        }
    }
    HoHet = totHoHet/(totHoHet+totHoHo);
  }


  else if ((params.flagX) && (params.diphap == 2)) { // CROMOSOMA X PHASED
    binmax = 0;  // CROMOSOMA X PHASED
    // Separate the two genomes of females by copying the second to another position:
    conta=eneind;
    for (i=0;i<eneind;++i){
      if (sexo[i]){
        for (j=0;j<popInfo.numLoci;++j){
          indi[0][conta][p[j]] = indi[1][i][p[j]];
        }
        ++conta;
      }
    }
    eneindX=conta; // Number of haploid genomes
    eneunits = eneindX*(eneindX-1)/2;
    for (i = 0; i < eneindX - 1; ++i) {
      #pragma omp parallel for private(crec, Hoant, Honow, xj, xii, xj3, xconta, xconta2, xroh, xiniciocrom, xcenmor)
      for (xii = i + 1; xii < eneindX; ++xii) {
        Hoant=false;
        for (xj = 0; xj < popInfo.numcromo; ++xj) {
          xconta = rangocromo[xj];
          xconta2 = xconta;
          xroh = false;
          xiniciocrom = true;
          while (xconta2 < rangocromo[xj + 1]) {
            Honow=(indi[0][i][p[xconta2]] == indi[0][xii][p[xconta2]]);
            if (Hoant){
              if (Honow){++contaHoHo[xii];}
              else{++contaHoHet[xii];}
            }
            Hoant=Honow;
            if (Honow) {
              if (!xroh) {
                xroh = true;
                xconta = xconta2;
                if (xconta > rangocromo[xj]) {
                  --xconta;
                } // New criterion: distance between two het sites
              }
            } else {
              if (xroh) {
                xcenmor = fabs(posicM[p[xconta2]] - posicM[p[xconta]]);

                // // Recombination reduction due to males:
                //   xcenmor /= 0.6666666;
                // crec=(1.0-exp(-xcenmor/50.0))/2.0*0.6666666;
                // xcenmor=-log(1.0-2.0*crec)*50.0;

                if (xiniciocrom) {
                  xcenmor += cortel;
                } // Telomere correction at start
                if (xcenmor >= cMbase && xcenmor <= cMtope) {
                  xj3 = floor(float((xcenmor - cMbase) / params.binsize));
                  if (xj3 < MAXBINS) {
                    xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                                 params.binsize +
                             params.binsize / 2.0 +
                             cMbase; // Use central bin mark
                    if (xj3>Matbinmax[xii]){
                        Matbinmax[xii] = xj3;
                    }
                // All these updates are atomic:
                      #pragma omp atomic
                      ++nbinroh[xj3];
                      #pragma omp atomic
                      cMbinroh[xj3] += xcenmor;
                  }
                }
                xroh = false;
              }
              xiniciocrom = false;
            }
            ++xconta2;
          }
          if (xroh) {
            xcenmor = fabs(posicM[p[xconta2 - 1]] - posicM[p[xconta - 1]]);

            // // Recombination reduction due to males:
            //       xcenmor /= 0.6666666;
            // crec=(1.0-exp(-xcenmor/50.0))/2.0*0.6666666;
            // xcenmor=-log(1.0-2.0*crec)*50.0;

            xcenmor += cortel + sitedist; // correcc. Telom. Final
            if (xcenmor >= cMbase && xcenmor <= cMtope) {
              xj3 = floor(float((xcenmor - cMbase) / params.binsize));
              if (xj3 < MAXBINS) {
                xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                             params.binsize +
                         params.binsize / 2.0 +
                         cMbase; // Use central bin mark
              if (xj3>Matbinmax[xii]){
                  Matbinmax[xii] = xj3;
              }
          // All these updates are atomic:
                #pragma omp atomic
                ++nbinroh[xj3];
                #pragma omp atomic
                cMbinroh[xj3] += xcenmor;
              }
            }
          }
        }
      }
    }
    totHoHo=0;
    totHoHet=0;
    for (i = 0; i < eneindX - 1; ++i) {
        totHoHo += contaHoHo[i];
        totHoHet += contaHoHet[i];
        if (Matbinmax[i]>binmax){
          binmax=Matbinmax[i];
        }
    }
    HoHet = totHoHet/(totHoHet+totHoHo);
  }


  else if ((params.flagX) && (params.diphap == 1)) { // X CHROMOSOME MALES ONLY
    binmax = 0; // X CHROMOSOME MALES ONLY
    if (enemachos>1) {eneunits = enemachos*(enemachos-1)/2;}// HAPLOIDS (MALES)
    for (i = 0; i < eneind - 1; ++i) {
      if (!sexo[i]){
          #pragma omp parallel for private(crec, Hoant, Honow, xj, xii, xj3, xconta, xconta2, xroh, xiniciocrom, xcenmor)
          for (xii = i + 1; xii < eneind; ++xii) {
            if (!sexo[xii]){
              Hoant=false;
              for (xj = 0; xj < popInfo.numcromo; ++xj) {
                  xconta = rangocromo[xj];
                  xconta2 = xconta;
                  xroh = false;
                  xiniciocrom = true;
                  while (xconta2 < rangocromo[xj + 1]) {
                    Honow=(indi[0][i][p[xconta2]] == indi[0][xii][p[xconta2]]);
                    if (Hoant){
                      if (Honow){++contaHoHo[xii];}
                      else{++contaHoHet[xii];}
                    }
                    Hoant=Honow;
                    if (Honow) {
                      if (!xroh) {
                        xroh = true;
                        xconta = xconta2;
                        if (xconta > rangocromo[xj]) {
                          --xconta;
                        } // New criterion: distance between two het sites
                      }
                    } else {
                      if (xroh) {
                        xcenmor = fabs(posicM[p[xconta2]] - posicM[p[xconta]]);

                  //       // Recombination reduction due to males:
                  // xcenmor /= 0.6666666;
                  //       crec=(1.0-exp(-xcenmor/50.0))/2.0*0.6666666;
                  //       xcenmor=-log(1.0-2.0*crec)*50.0;

                        if (xiniciocrom) {
                          xcenmor += cortel;
                        } // Telomere correction at start
                        if (xcenmor >= cMbase && xcenmor <= cMtope) {
                          xj3 = floor(float((xcenmor - cMbase) / params.binsize));
                          if (xj3 < MAXBINS) {
                            xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                                        params.binsize +
                                    params.binsize / 2.0 +
                                    cMbase; // Use central bin mark
                            if (xj3>Matbinmax[xii]){
                                Matbinmax[xii] = xj3;
                            }
                        // All these updates are atomic:
                              #pragma omp atomic
                              ++nbinroh[xj3];
                              #pragma omp atomic
                              cMbinroh[xj3] += xcenmor;
                          }
                        }
                        xroh = false;
                      }
                      xiniciocrom = false;
                    }
                    ++xconta2;
                  }
                  if (xroh) {
                    xcenmor = fabs(posicM[p[xconta2 - 1]] - posicM[p[xconta - 1]]);

                  //   // Recombination reduction due to males:
                  // xcenmor /= 0.6666666;
                  //   crec=(1.0-exp(-xcenmor/50.0))/2.0*0.6666666;
                  //   xcenmor=-log(1.0-2.0*crec)*50.0;

                    xcenmor += cortel + sitedist; // correcc. Telom. Final
                    if (xcenmor >= cMbase && xcenmor <= cMtope) {
                      xj3 = floor(float((xcenmor - cMbase) / params.binsize));
                      if (xj3 < MAXBINS) {
                        xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                                    params.binsize +
                                params.binsize / 2.0 +
                                cMbase; // Use central bin mark
                      if (xj3>Matbinmax[xii]){
                          Matbinmax[xii] = xj3;
                      }
                  // All these updates are atomic:
                        #pragma omp atomic
                        ++nbinroh[xj3];
                        #pragma omp atomic
                        cMbinroh[xj3] += xcenmor;
                      }
                    }
                  }
              }
            }
          }
      }
    }
    totHoHo=0;
    totHoHet=0;
    for (i = 0; i < eneind; ++i) {
        totHoHo += contaHoHo[i];
        totHoHet += contaHoHet[i];
        if (Matbinmax[i]>binmax){
          binmax=Matbinmax[i];
        }
    }
    HoHet = totHoHet/(totHoHet+totHoHo);
  }


  else if ((params.flagX) && (params.diphap == 0)) { // X CHROMOSOME UNPHASED
    binmax = 0; // PART 1: TWO CHROMOSOMES X UNPHASED (FEMALES)
    eneunits = enehembras;
    #pragma omp parallel for private(crec, Hoant, Honow, xj, xj3, xconta, xconta2, xroh, xiniciocrom, xcenmor)
    for (i = 0; i < eneind; ++i) {
      if (sexo[i]){
          Hoant=false;
          for (xj = 0; xj < popInfo.numcromo; ++xj) {
            xconta = rangocromo[xj];
            xconta2 = xconta;
            xroh = false;
            xiniciocrom = true;
            while (xconta2 < rangocromo[xj + 1]) {
              Honow=(indi[0][i][p[xconta2]] == indi[1][i][p[xconta2]]);
              if (Hoant){
                if (Honow){++contaHoHo[i];}
                else{++contaHoHet[i];}
              }
              Hoant=Honow;
              if (Honow) {
                if (!xroh) {
                  xroh = true;
                  xconta = xconta2;
                  if (xconta > rangocromo[xj]) {
                    --xconta;
                  } // New criterion: distance between two het sites
                }
              }
              else {
                if (xroh) {
                  xcenmor = fabs(posicM[p[xconta2]] - posicM[p[xconta]]);

                  // // Recombination reduction due to males:
                  // xcenmor /= 0.6666666;
                  // crec=(1.0-exp(-xcenmor/50.0))/2.0*0.6666666;
                  // xcenmor=-log(1.0-2.0*crec)*50.0;

                  if (xiniciocrom) {
                    xcenmor += cortel;
                  } // Telomere correction at start
                  if ((xcenmor >= cMbase) && (xcenmor <= cMtope)) {
                    xj3 = floor(float((xcenmor - cMbase) / params.binsize));
                    if (xj3 < MAXBINS) {
                      xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                                  params.binsize +
                              params.binsize / 2.0 +
                              cMbase; // Use central bin mark
                      if (xj3>Matbinmax[i]){
                          Matbinmax[i] = xj3;
                      }
                  // All these updates are atomic:
                        #pragma omp atomic
                        ++nbinroh[xj3];
                        #pragma omp atomic
                        cMbinroh[xj3] += xcenmor;
                    }
                  }
                  xroh = false;
                }
                xiniciocrom = false;
              }
              ++xconta2;
            }
            if (xroh) {
              xcenmor = fabs(posicM[p[xconta2 - 1]] - posicM[p[xconta - 1]]);

              // // Recombination reduction due to males:
              //     xcenmor /= 0.6666666;
              // crec=(1.0-exp(-xcenmor/50.0))/2.0*0.6666666;
              // xcenmor=-log(1.0-2.0*crec)*50.0;

              xcenmor += cortel + sitedist; // correcc. Telom. Final
              if (xcenmor >= cMbase && xcenmor <= cMtope) {
                xj3 = floor(float((xcenmor - cMbase) / params.binsize));
                if (xj3 < MAXBINS) {
                  xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                              params.binsize +
                          params.binsize / 2.0 +
                          cMbase; // Use central bin mark
                  if (xj3>Matbinmax[i]){
                      Matbinmax[i] = xj3;
                  }
              // All these updates are atomic:
                    #pragma omp atomic
                    ++nbinroh[xj3];
                    #pragma omp atomic
                    cMbinroh[xj3] += xcenmor;
                }
              }
            }
          }
      }

    }
    totHoHo=0;
    totHoHet=0;
    for (i = 0; i < eneind; ++i) {
        totHoHo += contaHoHo[i];
        contaHoHo[i]=0;
        totHoHet += contaHoHet[i];
        contaHoHet[i]=0;
        if (Matbinmax[i]>binmax){
          binmax=Matbinmax[i];
        }
    }
    if (enemachos>1) {eneunits += enemachos*(enemachos-1)/2;}// PART 2: HAPLOIDS (MALES)
    for (i = 0; i < eneind - 1; ++i) {
      if (!sexo[i]){
          #pragma omp parallel for private(crec, Hoant, Honow, xj, xii, xj3, xconta, xconta2, xroh, xiniciocrom, xcenmor)
          for (xii = i + 1; xii < eneind; ++xii) {
            if (!sexo[xii]){
              Hoant=false;
              for (xj = 0; xj < popInfo.numcromo; ++xj) {
                  xconta = rangocromo[xj];
                  xconta2 = xconta;
                  xroh = false;
                  xiniciocrom = true;
                  while (xconta2 < rangocromo[xj + 1]) {
                    Honow=(indi[0][i][p[xconta2]] == indi[0][xii][p[xconta2]]);
                    if (Hoant){
                      if (Honow){++contaHoHo[xii];}
                      else{++contaHoHet[xii];}
                    }
                    Hoant=Honow;
                    if (Honow) {
                      if (!xroh) {
                        xroh = true;
                        xconta = xconta2;
                        if (xconta > rangocromo[xj]) {
                          --xconta;
                        } // New criterion: distance between two het sites
                      }
                    } else {
                      if (xroh) {
                        xcenmor = fabs(posicM[p[xconta2]] - posicM[p[xconta]]);

                  //       // Recombination reduction due to males:
                  // xcenmor /= 0.6666666;
                  //       crec=(1.0-exp(-xcenmor/50.0))/2.0*0.6666666;
                  //       xcenmor=-log(1.0-2.0*crec)*50.0;

                        if (xiniciocrom) {
                          xcenmor += cortel;
                        } // Telomere correction at start
                        if (xcenmor >= cMbase && xcenmor <= cMtope) {
                          xj3 = floor(float((xcenmor - cMbase) / params.binsize));
                          if (xj3 < MAXBINS) {
                            xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                                        params.binsize +
                                    params.binsize / 2.0 +
                                    cMbase; // Use central bin mark
                            if (xj3>Matbinmax[xii]){
                                Matbinmax[xii] = xj3;
                            }
                        // All these updates are atomic:
                              #pragma omp atomic
                              ++nbinroh[xj3];
                              #pragma omp atomic
                              cMbinroh[xj3] += xcenmor;
                          }
                        }
                        xroh = false;
                      }
                      xiniciocrom = false;
                    }
                    ++xconta2;
                  }
                  if (xroh) {
                    xcenmor = fabs(posicM[p[xconta2 - 1]] - posicM[p[xconta - 1]]);

                  //   // Recombination reduction due to males:
                  // xcenmor /= 0.6666666;
                  //   crec=(1.0-exp(-xcenmor/50.0))/2.0*0.6666666;
                  //   xcenmor=-log(1.0-2.0*crec)*50.0;

                    xcenmor += cortel + sitedist; // correcc. Telom. Final
                    if (xcenmor >= cMbase && xcenmor <= cMtope) {
                      xj3 = floor(float((xcenmor - cMbase) / params.binsize));
                      if (xj3 < MAXBINS) {
                        xcenmor = floor(float((xcenmor - cMbase) / params.binsize)) *
                                    params.binsize +
                                params.binsize / 2.0 +
                                cMbase; // Use central bin mark
                      if (xj3>Matbinmax[xii]){
                          Matbinmax[xii] = xj3;
                      }
                  // All these updates are atomic:
                        #pragma omp atomic
                        ++nbinroh[xj3];
                        #pragma omp atomic
                        cMbinroh[xj3] += xcenmor;
                      }
                    }
                  }
              }
            }
          }
      }
    }
    for (i = 0; i < eneind; ++i) {
        totHoHo += contaHoHo[i];
        totHoHet += contaHoHet[i];
        if (Matbinmax[i]>binmax){
          binmax=Matbinmax[i];
        }
    }
    HoHet = totHoHet/(totHoHet+totHoHo);
  }

  if (params.HoHet){
        FactorHoHetH = totHoHet/(totHoHet+totHoHo);
    }
  else {
      FactorHoHetH=1;
    }

  params.progress.SetTaskProgress(7);
  // params.progress.PrintProgress();
  fp = (1.0 + fs * (2.0 * eneind - 1.0)) / (2.0 * eneind - 1.0 + fs);

  std::stringstream salida2;
  salida2 << std::fixed << std::setprecision(0);
  salida2 << "# Sample size:"
          << "\n"
          << eneind << "\n";
          sample_size=eneind;
//  salida2 << "# Genome type (0: dip unphased, 1: hap, 2: dip phased):\n"<< params.diphap << "\n";
//          diphap=params.diphap;
  if (params.flagX){
  salida2 << "# X Chromosome (0: females unphased and males, 1: only males, 2: females phased and males):\n"<< params.diphap << "\n";
  }
  else{
  salida2 << "# Autosomes (0: diploids unphased, 1: haploids, 2: diploid phased):\n"<< params.diphap << "\n";
  }
  salida2 << std::fixed << std::setprecision(6);
  salida2 << "# Genome size in cM:"
          << "\n"
          << popInfo.Mtot * 100 << "\n";
          genomesize=popInfo.Mtot*100;
  salida2 << std::fixed << std::setprecision(0);
  salida2 << "# No. of chromosomes:"
          << "\n"
          << numcromo1 << "\n";
          ncrom=numcromo1;
  salida2 << "# Number of segments with continuous data over 40 cM or more:\n";
  salida2 << std::fixed << std::setprecision(0);
  salida2 << numcromo2 << "\n";
  salida2 << std::fixed << std::setprecision(6);
  salida2 << "# Average distance between consecutive sites in cM:"
          << "\n"
          << sitedist << "\n";
          dcM=sitedist;
  salida2 << "# Bin width in cM:"
          << "\n"
          << params.binsize << "\n";
          binwidthcM=params.binsize;
  salida2 << "# Mutation rate per cM:"
          << "\n"
          << mucM << "\n";
  salida2 << "# Genotyping error rate per cM:"
          << "\n"
          << errorcM << "\n";
  salida2 << "# Heterozygosity (>MAF sites):"
          << "\n"
          << Het_total << "\n";
  salida2 << "# Heterozygosity after a homozygous site (>MAF sites):"
          << "\n"
          << HoHet << "\n";
          // H=Het_total;
          H=HoHet;
  // salida2 << "# d(cM)/H ranked values:"
  //         << "\n"
  //         << dH_vals[0] << "\n"
  //         << dH_vals[1] << "\n"
  //         << dH_vals[2] << "\n"
  //         << dH_vals[3] << "\n"
  //         << dH_vals[4] << "\n";
  salida2 << "# Bin_class_in_cM:\t"
          << "Total_No._of_segments:\t"
          << "Prop_ofGenome_IBD_due_to_this_class:\n";

  basebinsize = params.lc;
  for (j3 = 0; j3 <= binmax; ++j3) {
    if (nbinroh[j3] > 0) {
      cMbinroh[j3] /= nbinroh[j3];
      // Size of the roh in cM (that is the bin class)
      salida2 << std::fixed << std::setprecision(6);
      salida2 << cMbinroh[j3];
      xrohval[j3]=cMbinroh[j3];
      // // Total number of rohs found over all the genomes and over all the
      // individuals for this bin
      salida2 << std::fixed << std::setprecision(0);
      salida2 << "\t" << nbinroh[j3];
      xnbin[j3]=nbinroh[j3];
      ynbin[j3]=nbinroh[j3];
      // Average Proportion of the genome that is IBD due to this bin (<< THIS
      // IS THE GOOD ONE)
      if (params.flagX){
        // xrohobs[j3]=(cMbinroh[j3] * nbinroh[j3] / eneunits) / (popInfo.Mtot * 0.66666 * 100);
        xrohobs[j3]=(cMbinroh[j3] * nbinroh[j3] / eneunits) / (popInfo.Mtot * 100);
      }
      else{
        xrohobs[j3]=(cMbinroh[j3] * nbinroh[j3] / eneunits) / (popInfo.Mtot * 100);
      }
      yrohobs[j3]=xrohobs[j3];
      salida2 << std::fixed << std::setprecision(6);
      salida2 << "\t"<< xrohobs[j3]<< "\n";
    }
    else {
      xrohval[j3]=basebinsize;
      xnbin[j3]=0;
      ynbin[j3]=0;
      xrohobs[j3]=0;
      yrohobs[j3]=0;
      salida2 << basebinsize << "\t0\t0\n";
    }
    basebinsize += params.binsize;
//    std::cout<<j3<<":"<<xnbin[j3]<<std::endl;
  }
  maxnlin=binmax+1;
  const std::string fichsal2 = fich + "_NeROH_STATS.txt";
  std::ofstream outputFile2;
  outputFile2.open(fichsal2);
  outputFile2 << salida2.str();
  outputFile2.close();

  params.progress.SetCurrentTask(1, "Estimating Ne");
  // If maxgen changes, this should change too
  params.progress.InitCurrentTask(NREPETICIONES * 330);
  params.progress.SetTaskProgress(0);

  const std::string fichsal4 = fich + "_NeROH";
  neroh(fichsal4);

  tpas = (omp_get_wtime() - tini);
  std::stringstream salida;
  salida << "# (NeROH v1.0)\n";
  salida << "# Command:";
  for (i = 0; i < argc; ++i) {
    salida << " " << argv[i];
  }
  salida << "\n";
  salida << "# Running time: ";
  salida << (float(tpas)) << "sec\n";
  salida << "#\n";
  salida << "# INPUT PARAMETERS:\n";
  salida << std::fixed << std::setprecision(0);
  if (params.flagX){
    salida << "# X Chromosome (0: females unphased and males, 1: only males, 2: females phased and males):\n";
    salida<< params.diphap << "\n";
  }
  else{
    salida << "# Autosomes (0: diploids unphased, 1: haploids, 2: diploid phased):\n";
    salida<< params.diphap << "\n";
  }
  salida << "# Number of chromosomes:\n";
  salida << std::fixed << std::setprecision(0);
  salida << popInfo.numcromo << "\n";
  salida << "# Number of segments with continuous data over 40 cM or more:\n";
  salida << std::fixed << std::setprecision(0);
  salida << numcromo2 << "\n";
  salida << "# Genome size in Morgans:\n";
  salida << std::fixed << std::setprecision(4);
  salida << popInfo.Mtot << "\n";
  salida << "# Genome size in Mb:\n";
  salida << std::fixed << std::setprecision(2);
  salida << popInfo.Mbtot << "\n";
  salida << "# Number of individuals in the input file:\n";
  salida << eneind << "\n";
  salida << "# Total number of markers in the input file:\n";
  salida << std::fixed << std::setprecision(0);
  salida << popInfo.fichLoci << "\n";
  salida << "# Number of markers with frequencies>MAF in the input file:\n";
  salida << eneloc << "\n";
  salida << "# Mutation rate per cM:\n";
  salida << std::fixed << std::setprecision(6);
  salida << mucM << "\n";
  salida << "# Genotyping error rate per cM:\n";
  salida << std::fixed << std::setprecision(6);
  salida << errorcM << "\n";
  salida << "# Estimated Fis value of the population (deviation from H-W proportions):\n";
  if (params.diphap == 1) {
    salida << "Not applicable \n";
  } else {
    salida << fp << "\n";
  }
  salida << "# Heterozygosity (only > MAF sites):\n";
  salida << Het_total << "\n";
  // salida << "# Heterozygosity under H-W equilibrium (all sites):\n";
  // salida << Het_total*eneloc/popInfo.numLoci<<"\n";

  if (params.printToStdOut) {
    std::string output = salida.str();
    std::cout << output << std::endl;
  } else {
    std::string fichsal = fich + "_NeROH_INPUT.txt";
    std::ofstream outputFile;
    outputFile.open(fichsal);
    outputFile << salida.str();
    outputFile.close();
  //      std::cout << " Finished preprocessing. Output file " << fichsal << " generated\n";
  }
  std::cout << "  End of processing.\n";
  std::cout << "  Total run time in seconds: " << (float(tpas)) << "sec.\n\n";

  std::remove(fichProgress.c_str());
  return 0;
}

int neroh(const std::string &fichentra) {
  double a, aa, b, bb, d, dd, maxc, Ne0, increm, producto = 1;
  double rohlow, rohhigh, minn, sumnbins = 0, dpob, dist1, dist2;
  double s, BB, ll, sk, sm, QQ;
  double h1, h2, DD, R, Rant;
  int i, ii, jj, j, imaxc, incre1 = 0, incre2 = 0, incre3 = 0;
  int mini, sizebins = 1;
  bool flagdebug = false, flagR, hayalgoquecomprimir;
  struct tm *loctime;
  double tinicio, tfinal;
  clock_t tini, tiniabs, tpas;

  std::string fichero = "";

  fichero = fichentra + "";
  std::string fichero_sal_NeH = fichentra + "_Ne.txt";
  std::string fichero_sal_roh = fichentra + "_DISTRIB.txt";
  std::string fichero_sal_evol = fichentra + "_evol"; // --debug
  std::string fichero_caca = fichentra + "_caca";     // --debug
  std::ofstream salida;

  tinicio = omp_get_wtime();

  for (i=0;i<5;++i){ // Ahora va en morgans
      dH_vals[i] /= 100;
  }

  for (i=0;i<5;++i){ // Apply the Ho/Het correction factor here
      dH_vals[i] /= FactorHoHetH;
  }
  d = dcM / 100; // Distance in Morgans
  wl = binwidthcM / 100; // Width in Morgans
  mu = mucM * 100;
  err = errorcM * 100;
  conta0=0;
  for (nlin=0;nlin<maxnlin;++nlin){
      xpeso[nlin] = 1;
//  std::cout<<nlin<<":"<<xnbin[nlin]<<std::endl;
      if (xnbin[nlin] > 0) {
        conta0 = 0;
      }
      else {
        ++conta0;
        if (conta0 > 1) { // Stop reading at two or more consecutive lines with zeros
          break;
        }
      }
  }
  nlin -= conta0;// Keep lines up to one before the series of zeros
  ++nlin;
//  std::cout<<maxnlin<<std::endl;
//  std::cout<<nlin<<std::endl;

  // Moving averages with +-2 bin grouping
  // First find the first bin with 10 or fewer cases
//   std::cout<<":1:"<<nlin<<std::endl;
  ii = nlin;
  for (i = 0; i < nlin; ++i) {
    if (ynbin[i] <= 10) {
      ii = i;
      break;
    }
  }
  // Then apply 5-point moving averages
  for (i = ii; i < nlin - 2; ++i) {
    a = 0;
    suma = 0;
    for (j = i - 2; j <= i + 2; ++j) {
      suma += yrohobs[j];
      a += ynbin[j];
    }
    xrohobs[i] = suma / (5);
    xnbin[i] = a / (5);
  }
  // For the penultimate bin, always use a 3-point moving average
  xrohobs[nlin - 2] =
      (yrohobs[nlin - 3] + yrohobs[nlin - 2] + yrohobs[nlin - 1]) / 3;
  xnbin[nlin - 2] = (ynbin[nlin - 3] + ynbin[nlin - 2] + ynbin[nlin - 1]) / 3;
  --nlin; // Remove the last bin

  // Sample the first bins to get an initial Ne estimate
  // Note: here Ne is diploid, unlike in GONE
  Nemed = 0;
  conta = 0;
  for (i = 5; i < 11; ++i) {
    if (i >= nlin) {
      break;
    }
    Ne0 = wl / (2 * xrohobs[i] * xrohval[i] * xrohval[i] / 10000);
    if (Ne0 > 10) {
      conta += 1;
      Nemed += (Ne0);
    }
  }
  if (conta > 0) {
    Nemed = (Nemed / conta);
  } else {
    Nemed = 1000;
  }
  if (Nemed < 10) {
    Nemed = 10;
  }

  for (i = 0; i < nlin; ++i) { // Compute average coalescence age for each bin
    ll = xrohval[i] / 100;
    a = 2 * ll + 0.5 / Nemed - 4 * d / (H - d);
    if (a > 0) {
      edadrohs[i] = 3 / (a);
    } else {
      edadrohs[i] = 3 / (2 * ll + .5 / Nemed);
    }
  }

  // Readjust class marks to correspond to class means (see theory):
  flagR = true;
  R = 0;
  Rant = 0;
  for (i = 0; i < nlin; ++i) {
    if (flagR && (i < (nlin - 1))) {
      h1 = xnbin[i];
      h2 = xnbin[i + 1];
      if (h1 >= h2) {
        DD = (wl / 3 * (h1 - h2) / 2 + wl / 2 * h2) / ((h1 - h2) / 2 + h2);
        R = (wl / 2 - DD) * 100; // Now in cM
      } else {
        flagR = false;
        R = Rant;
      }
      if ((Rant < R) && (i > 0)) {
        R = (Rant + R) / 2;
      }
      xrohval[i] -= R;
      Rant = R;
    } else {
      xrohval[i] -= R;
    }
  }

  // Merge consecutive lines whose age difference is less than 1 generation
  hayalgoquecomprimir = true;
  while (hayalgoquecomprimir) {
    hayalgoquecomprimir = false;
    for (i = nlin - 2; i > 0; (--i)) {
      if ((edadrohs[i] - edadrohs[i + 1]) < 1) {
        xrohval[i] = (xpeso[i] * xrohval[i] + xpeso[i + 1] * xrohval[i + 1]) /
                     (xpeso[i] + xpeso[i + 1]);
        xrohobs[i] = (xpeso[i] * xrohobs[i] + xpeso[i + 1] * xrohobs[i + 1]) /
                     (xpeso[i] + xpeso[i + 1]);
        edadrohs[i] =
            (xpeso[i] * edadrohs[i] + xpeso[i + 1] * edadrohs[i + 1]) /
            (xpeso[i] + xpeso[i + 1]);
        xnbin[i] = xnbin[i] + xnbin[i + 1];
        xpeso[i] = xpeso[i] + xpeso[i + 1];
        // Compress the list
        --nlin;
        for (j = i + 1; j < nlin; ++j) {
          xrohval[j] = xrohval[j + 1];
          xrohobs[j] = xrohobs[j + 1];
          xnbin[j] = xnbin[j + 1];
          xpeso[j] = xpeso[j + 1];
          edadrohs[j] = edadrohs[j + 1];
        }
        hayalgoquecomprimir = true;
        --i;
      }
      // else{
      //     break;
      // }
    }
  }
//  std::cout<<":2:"<<nlin<<std::endl;

  // Examine trailing lines and remove those with too few representatives
  j=nlin;
  for (i = nlin - 1; i > 0; (--i)) {
    if ((params.flaghc && (xnbin[i] < 5)) || (!params.flaghc && (xnbin[i] < 25))){
        --j;
    }
    else{
        break;
    }
  }
  nlin=j;

  minedadrohs = int(edadrohs[nlin - 1]);
  if (nlin < 6) {  // (was 8)
    std::cerr << " Too few bins." << std::endl;
    return 1;
  }

  if (nlin > 45) {
    incre1= 1;
  } // era 40
  if (nlin > 30) {
    incre2 = 1;
  } // era 20
  if (nlin > 15) {
    incre3 = 1;
  } // era 10

  sumnbins = 0;
  for (i = 0; i < nlin; ++i) {
    nbin[i] = xnbin[i];
    sumnbins += nbin[i];
    rohval[i] = xrohval[i];
    rohobs[i] = xrohobs[i];
  }

  tiniabs = clock();
  tpas = clock() - tiniabs;

  tinicio = omp_get_wtime();
  tini = clock();
  // std::cout << " Start of processing: " << fichero << std::endl;

  // Start of multithreading
  int maxThreads = std::min(omp_get_max_threads(), NREPETICIONES);
  int numThreads = std::min(maxThreads, 16);
  omp_set_num_threads(numThreads);

  std::vector<std::vector<double>> t_acusumNe(NREPETICIONES);
  std::vector<std::vector<double>> t_acurohprd(NREPETICIONES);
  std::vector<int> tcounter(NREPETICIONES);
  int progressCounter = 0;

  // Grow the first dimension
  // t_acusumNe.resize(numThreads);
  // t_acurohprd.resize(numThreads);
  // tcounter.resize(numThreads);

//  std::cout<<":3:"<<nlin<<std::endl;

  // Initialize the vectors
  // for (int icounter=0; icounter < numThreads; ++icounter) {
  for (int icounter=0; icounter < NREPETICIONES; ++icounter) {
    t_acusumNe[icounter].resize(ngenmax);
    t_acurohprd[icounter].resize(nlinmax);
    tcounter[icounter] = 0;
  }
  gmax = 300;  // (was 200)
  gmax2 = 100;  // (was 100)
  gmax3 = 200;  // (was 150)
  //                              OPTIMIZATION PARAMETERS
  //**************************************************************************************
  flagsolape = true;
  flagporbloques = true;
  resolucion = 4;  // (was 5)
  minedadrohs -= resolucion + 0;  // (was 2)
  if (minedadrohs < 0) {
    minedadrohs = 0;
  }
  topeposigen = gmax2 - 2 * resolucion;
  topeposiinicio1 = 35;
  topeposiinicio2 = 70;
  topeposiinicio3 = 100; // era 95
  topeposiinicio4 = 130; // era 120
  if (minedadrohs > topeposiinicio1) {
    minedadrohs = topeposiinicio1 - 15;  // (was 3)
  }
  if (topeposiinicio1 > topeposigen) {
    topeposiinicio1 = topeposigen;
  }
  if (topeposiinicio2 > topeposigen) {
    topeposiinicio2 = topeposigen;
  }
  if (topeposiinicio3 > topeposigen) {
    topeposiinicio3 = topeposigen;
  }
  if (topeposiinicio4 > topeposigen) {
    topeposiinicio4 = topeposigen;
  }
  topeposiinicio = topeposigen;
  topesalto = 20.0; // era 9
  invtopesalto = 1.0 / topesalto;
  topesalto2 = 100; // era 50
  invtopesalto2 = 1.0 / topesalto2;
  ndesini = 3000;
  ndes = 1000;
  efectomutsuave = 0.02;
  frecinversion = 0.3; // Inversion probability for last two segments
  tercio1 = 10;        // Top-tier parents (elite, not modified)
  tercio2 =90; // Second-tier parents (replaced by offspring if offspring are better)
  tercio12 = tercio1 + tercio2;
  nhijos = tercio12;

  omp_lock_t lock;
  omp_init_lock(&lock);

#pragma omp parallel for private(BB, Ne, Neblock, Necons, QQ, SC, SCbest, SCmed, SCmedanterior, aa, ancho1, ancho2, bb, bicho, bichoH, bichoP, conta, conta2, contagen, des, efecto, efectomut, efectomutlateral, frecmut, frecmutlateral, frecnorec, frecrec, gen, i, ii, ind1, ind2, ind3, indexmaxSC, indexminSC, j, jj, maxgen, maxsegmentos, minSC, nsegmed, nsegbest, nsegmentos, posi1, posi2, posi3, posiblock, posigen, rohprd, s, segdesdehasta, sk, sm, sumNe, tope)
  for (REPE = 0; REPE < NREPETICIONES; ++REPE) {
    // int tid = omp_get_thread_num();
    int tid = REPE;

    // Random pre-loading with selection
    bicho.nseg = 4;         // ++
    bicho.segbl[0] = 0;     // ++
    bicho.segbl[3] = gmax3; // ++
    bicho.segbl[4] = gmax;  // ++
    for (ii = 0; ii < ndesini; ++ii) {
      for (;;) {
        aa = uniforme01(genera);
        // if (aa < 0.2) {
        //   posi1 = int(uniforme01(genera) *
        //               (topeposiinicio1 - resolucion - minedadrohs)) +
        //           resolucion + minedadrohs;
        //   posi2 = int(uniforme01(genera) *
        //               (topeposiinicio1 - resolucion - minedadrohs)) +
        //           resolucion + minedadrohs;
        // } else if (aa < 0.5) {
        //   posi1 = int(uniforme01(genera) *
        //               (topeposiinicio2 - resolucion - minedadrohs)) +
        //           resolucion + minedadrohs;
        //   posi2 = int(uniforme01(genera) *
        //               (topeposiinicio2 - resolucion - minedadrohs)) +
        //           resolucion + minedadrohs;
        // } else if (aa < 0.8) {
        //   posi1 = int(uniforme01(genera) *
        //               (topeposiinicio3 - resolucion - minedadrohs)) +
        //           resolucion + minedadrohs;
        //   posi2 = int(uniforme01(genera) *
        //               (topeposiinicio3 - resolucion - minedadrohs)) +
        //           resolucion + minedadrohs;
        // } else {
        //   posi1 = int(uniforme01(genera) *
        //               (topeposiinicio4 - resolucion - minedadrohs)) +
        //           resolucion + minedadrohs;
        //   posi2 = int(uniforme01(genera) *
        //               (topeposiinicio4 - resolucion - minedadrohs)) +
        //           resolucion + minedadrohs;
        // }

        if (aa < 0.3) {  // (was 0.2)
            posi1 = int(uniforme01(genera) * (topeposiinicio1))+resolucion;
            posi2 = int(uniforme01(genera) * (topeposiinicio1))+resolucion;
        } else if (aa < 0.7) {  // (was 0.5)
            posi1 = int(uniforme01(genera) * (topeposiinicio2))+resolucion;
            posi2 = int(uniforme01(genera) * (topeposiinicio2))+resolucion;
        } else if (aa < 0.9) {  // (was 0.8)
            posi1 = int(uniforme01(genera) * (topeposiinicio3))+resolucion;
            posi2 = int(uniforme01(genera) * (topeposiinicio3))+resolucion;
        } else {
            posi1 = int(uniforme01(genera) * (topeposiinicio4))+resolucion;
            posi2 = int(uniforme01(genera) * (topeposiinicio4))+resolucion;
        }


        if (posi1 > posi2) {
          posi3 = posi1;
          posi1 = posi2;
          posi2 = posi3;
        }
        if ((posi2 - posi1) > 5) {
          break;
        }
      }
      bicho.segbl[1] = posi1; // ++
      bicho.segbl[2] = posi2; // ++
      bicho.Nebl[0] = Nemed;
      tope = uniforme01(genera) * 2 * topesalto;
      //        tope=topesalto;  // Alternative to the previous line
      for (jj = 1; jj < 3; ++jj) {
        aa = 1 + uniforme01(genera) * (tope); // ERA2
        if (uniforme01(genera) < 0.66) {
          aa = 1 / aa;
        } // ERA0.5
        aa = Nemed * aa;
        if (aa > 10000000) {
          aa = 10000000;
        }
        if (aa < 5) {
          aa = 5;
        }
        bicho.Nebl[jj] = aa; // ++
      }
      aa = 1.0 + uniforme01(genera) * tope / 2; // era 0.5
      aa = Nemed * aa;
      if (aa > 10000000) {
        aa = 10000000;
      }
      if (aa < 5) {
        aa = 5;
      }
      bicho.Nebl[3] = aa; // ++

      SC = CalculaSC(bicho, rohprd);
      bicho.SCval = SC; // ++
      // Select parents
      if (ii < tercio12) {
        bichoP[ii] = bicho;
      } else {
        maxSC = 0; // Find worst SC among future parents
        indexmaxSC = 0;
        for (jj = 0; jj < tercio12; ++jj) {
          if (bichoP[jj].SCval > maxSC) {
            indexmaxSC = jj;
            maxSC = bichoP[jj].SCval;
          }
        }
        if (bicho.SCval < bichoP[indexmaxSC].SCval) {
          bichoP[indexmaxSC] = bicho;
        }
      }
    }
    // Sort the resulting parents
    for (ii = 0; ii < tercio12 - 1; ++ii) {
      indexminSC = ii;
      for (jj = ii + 1; jj < tercio12; ++jj) {
        if (bichoP[jj].SCval < bichoP[indexminSC].SCval) {
          indexminSC = jj;
        }
      }
      if (indexminSC != ii) { // swap
        bicho = bichoP[ii];
        bichoP[ii] = bichoP[indexminSC];
        bichoP[indexminSC] = bicho;
      }
    }
    // Begin evolutionary cycles
    contagen = 0;
    SCmedanterior = SC;
    // if (flagdebug) {
    //  salida.open(fichero_caca, std::ios::app);
    //  for (j = 0; j < bichoP[0].nseg; ++j) {
    //    salida << j << "\t" << bichoP[0].segbl[j] << "\t"
    //           << bichoP[0].Nebl[j] / 2.0 << "\n";
    //  }
    //  salida << "\n";
    //  salida.close();
    // }
    maxgen = 330;
    for (gen = 0; gen < maxgen; ++gen) {
      // if (flagdebug) {
      //  if (gen == int(double(gen / 100.0)) * 100) {
      //    salida.open(fichero_caca, std::ios::app);
      //    for (j = 0; j < bichoP[0].nseg; ++j) {
      //      salida << "-> " << gen << "\t" << j << "\t" << bichoP[0].segbl[j]
      //             << "\t" << bichoP[0].Nebl[j] / 2.0 << "\n";
      //    }
      //    salida << "\n";
      //    salida.close();
      //  }
      // }

      switch (gen) {
      case 0: {
        frecmut = 0.5;
        efectomut = 0.3;
        frecmutlateral = 0.8;
        efectomutlateral = 9;
        frecrec = 0.6;
        frecnorec = 0.3;
        maxsegmentos = 3;
        flagsolape = true;
        break;
      }
      case 30: {
        frecmut = 0.5;
        efectomut = 0.3;
        frecmutlateral = 0.8;
        efectomutlateral = 9;
        frecrec = 0.5;
        frecnorec = 0.3;
        maxsegmentos += incre2;
        flagsolape = true;
        break;
      }
      case 90: {
        frecmut = 0.5;
        efectomut = 0.3;
        frecmutlateral = 0.8;
        efectomutlateral = 9;
        frecrec = 0.4;
        frecnorec = 0.5;
        maxsegmentos += incre3;
        flagsolape = true;
        break;
      }
      case 180: {
        frecmut = 0.5;
        efectomut = 0.3;
        frecmutlateral = 0.8;
        efectomutlateral = 9;
        frecrec = 0.3;
        frecnorec = 0;
        maxsegmentos += incre2;
        flagsolape = true;
        break;
      }
      case 250: {
        frecmut = 0.5;
        efectomut = 0.2;
        frecmutlateral = 0.8;
        efectomutlateral = 9;
        frecrec = 0.2;
        frecnorec = 0.5;
        maxsegmentos += incre1;
        flagsolape = false;
        break;
      }
      case 290: {
        frecmut = 0.2;
        efectomut = 0.1;
        frecmutlateral = 0.2;
        efectomutlateral = 2;
        frecrec = 0.1;
        frecnorec = 0.6;
        maxsegmentos += incre1 + incre2;
        flagsolape = false;
        break;
      }
      case 310: {
        frecmut = 0.5;
        efectomut = 0.03;
        frecmutlateral = 0.2;
        efectomutlateral = 2;
        frecrec = 0.05;
        frecnorec = 0.7;
        maxsegmentos += incre1 + incre2 + incre3;
        flagsolape = false;
        break;
      }
      }
      for (des = 0; des < ndes;
           ++des) { // Main offspring generation loop
        for (;;) {
          for (int k = 0; k < 1000; ++k) {
            flag = true; // Boolean flag (should be made thread-private in OMP)
            ind1 = int(uniforme01(genera) *
                       tercio1); // Select from elite parents
            ind2 = int(uniforme01(genera) *
                       tercio12); // Select from all parents
            if (uniforme01(genera) < 0.5) {
              ind3 = ind1;
              ind1 = ind2;
              ind2 = ind3;
            } // Swap ind1 and ind2 with 50% probability
            // Recombination between and within parents with equal probability
            if (uniforme01(genera) < frecrec) {
              aa = uniforme01(genera);
              if (aa < 0.1) {
                topeposiinicio = topeposiinicio1;
              }
              else if (aa < 0.5) {
                topeposiinicio = topeposiinicio2;
              }
              else if (aa < 0.9) {
                topeposiinicio = topeposiinicio3;
              }
              else {
                topeposiinicio = topeposiinicio4;
              }
              if (flagporbloques) {
                do {
                  posiblock = int(uniforme01(genera) *
                                  (bichoP[ind1].nseg -
                                   1));
                  posigen =
                      bichoP[ind1].segbl[posiblock] +
                      int(uniforme01(genera) *
                          (bichoP[ind1].segbl[posiblock + 1] -
                           bichoP[ind1].segbl[posiblock]));
                } while (posigen >
                         topeposiinicio);
              }
              else {
                posigen =
                    resolucion + minedadrohs +
                    int(uniforme01(genera) *
                        (topeposiinicio - minedadrohs));
              }
              // Within: equalize individuals and center posigen

              if (uniforme01(genera) < .5) {
                ind2 = ind1;
              }
              for (jj = 0; jj <= bichoP[ind1].nseg; ++jj) {
                if (abs(bichoP[ind1].segbl[jj] - posigen) <= resolucion) {
                  flag = false;
                  break;
                }
              }
              if (flag) {
                for (jj = 0; jj <= bichoP[ind2].nseg; ++jj) {
                  if (abs(bichoP[ind2].segbl[jj] - posigen) <= resolucion) {
                    flag = false;
                    break;
                  }
                }
              }
              if (flag) {
                ii = 0;
                for (jj = 0; jj < bichoP[ind1].nseg; ++jj) {
                  if (bichoP[ind1].segbl[jj] < posigen) {
                    bicho.segbl[ii] = bichoP[ind1].segbl[jj];
                    bicho.Nebl[ii] = bichoP[ind1].Nebl[jj];
                    ++ii;
                  } else {
                    break;
                  }
                }
                for (jj = 1; jj < bichoP[ind2].nseg + 1; ++jj) {
                  if (bichoP[ind2].segbl[jj] > posigen) {
                    bicho.segbl[ii] = posigen;
                    bicho.Nebl[ii] = bichoP[ind2].Nebl[jj - 1];
                    ++ii;
                    break;
                  }
                }
                for (jj = 1; jj < bichoP[ind2].nseg; ++jj) {
                  if (bichoP[ind2].segbl[jj] > posigen) {
                    bicho.segbl[ii] = bichoP[ind2].segbl[jj];
                    bicho.Nebl[ii] = bichoP[ind2].Nebl[jj];
                    ++ii;
                  }
                }
                bicho.nseg = ii;
                bicho.segbl[ii] = bichoP[ind2].segbl[bichoP[ind2].nseg];
              }
            } else {
              bicho = bichoP[ind1];
            }
            if (flag) {
              break;
            }
          }
          // Merge two blocks with probability frecnorec using harmonic mean of Ne
          if (bicho.nseg > 3) {
            if ((uniforme01(genera) < frecnorec) ||
                (bicho.nseg > maxsegmentos)) {
              posiblock =
                  int(uniforme01(genera) * (bicho.nseg - 2)); // Merge randomly
              ancho1 = bicho.segbl[posiblock + 1] -
                       bicho.segbl[posiblock];
              ancho2 = bicho.segbl[posiblock + 2] -
                       bicho.segbl[posiblock + 1];
              aa = (ancho1 + ancho2) /
                   (ancho1 / bicho.Nebl[posiblock] +
                    ancho2 / bicho.Nebl[posiblock + 1]);
              flag = true;
              if (posiblock > 0) {
                bb = aa / bicho.Nebl[posiblock - 1];
                if (bb < invtopesalto2) {
                  flag = false;
                }
              }
              if (posiblock < bicho.nseg - 2) {
                bb = aa / bicho.Nebl[posiblock + 2];
                if (bb > topesalto2) {
                  flag = false;
                }
              }
              if (flag) {
                bicho.Nebl[posiblock] = aa;
                for (jj = posiblock + 1; jj <= bicho.nseg; ++jj) {
                  bicho.Nebl[jj] = bicho.Nebl[jj + 1];
                  bicho.segbl[jj] = bicho.segbl[jj + 1];
                }
                --bicho.nseg;
              }
            }
          }
          // Check segment count limit
          flag = true;
          if ((bicho.nseg > maxsegmentos) || (bicho.nseg < 3)) {
            flag = false;
          }
          if (flag) {
            break;
          }
        }
        // Ne mutation: randomly change Ne values of segments
        for (posiblock = 0; posiblock < bicho.nseg; ++posiblock) {
          if (uniforme01(genera) < frecmut) {
            efecto = bicho.Nebl[posiblock] * uniforme01(genera) * efectomut;
          }
          else {
            efecto =
                bicho.Nebl[posiblock] * uniforme01(genera) * efectomutsuave;
          }
          if (uniforme01(genera) < 0.5) {
            efecto = -efecto;
          }
          aa = bicho.Nebl[posiblock] + efecto;
          flag = true;
          if (posiblock > 0) {
            bb = aa / bicho.Nebl[posiblock - 1];
            if (bb < invtopesalto2) {
              flag = false;
            }
          }
          if (posiblock < bicho.nseg - 1) {
            bb = aa /
                 bicho.Nebl[posiblock + 1];
            if (bb > topesalto2) {
              flag = false;
            }
          }
          if (flag) {
            if (aa > 10000000) {
              aa = 10000000;
            }
            if (aa < 5) {
              aa = 5;
            }
            bicho.Nebl[posiblock] = aa;
          }
        }
        // Constrain last two segments to be similar:
        aa = bicho.Nebl[bicho.nseg - 1];
        bb = bicho.Nebl[bicho.nseg - 2];
        if (aa > (bb * 1.2)) {
          bicho.Nebl[bicho.nseg - 1] = bb * 1.2;
        }
        if (aa < (bb / 1.2)) {
          bicho.Nebl[bicho.nseg - 1] = bb / 1.2;
        }

        for (posiblock = 0; posiblock < bicho.nseg; ++posiblock) {
          if (bicho.Nebl[posiblock] < 5) {
            bicho.Nebl[posiblock] = 5;
          }
        }

        // Lateral mutation: randomly shift segment boundaries
        for (posiblock = 1; posiblock < bicho.nseg - 1; ++posiblock) {
          if (uniforme01(genera) < frecmutlateral) {
            efecto = int(uniforme01(genera) * efectomutlateral) + 1;
          }
          else {
            efecto = 0;
            //                    if (uniforme01(genera)<0.5){efecto=0;}
          }
          if (uniforme01(genera) < 0.5) { // ERA0.6
            if ((bicho.segbl[posiblock] - bicho.segbl[posiblock - 1]) >
                resolucion + efecto) {
              bicho.segbl[posiblock] -= efecto;
            }
          } else {
            if ((bicho.segbl[posiblock + 1] - bicho.segbl[posiblock]) >
                resolucion + efecto) {
              bicho.segbl[posiblock] += efecto;
            }
          }
        }


        // Additional lateral mutation: shift segment boundaries in block
        if (uniforme01(genera) < 0.05) {
          efecto = 1;
          jj=int(uniforme01(genera) * (bicho.nseg-1))+2;
          for (posiblock = jj; posiblock < bicho.nseg - 1; ++posiblock) {
              bicho.segbl[posiblock] += efecto;
              if ((bicho.segbl[posiblock + 1] - bicho.segbl[posiblock]) >resolucion) {
                (++(++efecto));
              }
          }
        }

        // // Additional lateral mutation: shift segment boundaries in block
        // if (uniforme01(genera) < 0.05) {
        //   efecto = 1;
        //   if (uniforme01(genera) < 0.5){efecto=-efecto;}
        //   jj=int(uniforme01(genera) * (bicho.nseg-1))+3;
        //   for (posiblock = jj; posiblock < bicho.nseg - 1; ++posiblock) {
        //     if (efecto<0) {
        //       if ((bicho.segbl[posiblock] - bicho.segbl[posiblock - 1]) > resolucion + efecto) {
        //         bicho.segbl[posiblock] += efecto;
        //         (--efecto);
        //       }
        //     } else {
        //       bicho.segbl[posiblock] += efecto;
        //       if ((bicho.segbl[posiblock + 1] - bicho.segbl[posiblock]) >resolucion) {
        //         (++(++efecto));
        //       }
        //     }
        //   }
        // }

        // End of OMP section
        SC = CalculaSC(bicho, rohprd);
        bicho.SCval = SC;
        if (des < nhijos) {
          bichoH[des] = bicho;
        } else {
          maxSC = -1;
          indexmaxSC = 0;
          for (i = 0; i < nhijos; ++i) {
            if (bichoH[i].SCval > maxSC) {
              maxSC = bichoH[i].SCval;
              indexmaxSC = i;
            }
          }
          if (SC < maxSC) {
            bichoH[indexmaxSC] = bicho;
          }
        }
      }
      // The top tier of bichoP parents is not modified.

      if (flagsolape) {
        // Insert the best offspring into the second tier of parents if they are better
        for (ii = tercio1; ii < tercio12; ++ii) {
          flag = true;
          minSC = maxdouble; // Find minimum SC among offspring
          indexminSC = 0;
          for (i = 0; i < nhijos; ++i) {
            if (bichoH[i].SCval < minSC) {
              minSC = bichoH[i].SCval;
              indexminSC = i;
            }
          }
          for (i = tercio1; i < tercio12;
               ++i) { // Check if this minimum is better than any parent in the second tier
            if (minSC < bichoP[i].SCval) {
              bichoP[i] = bichoH[indexminSC];
              bichoH[indexminSC].SCval = maxdouble;
              flag = false;
              break;
            }
          }
          if (flag) {
            break;
          }
        }
      } else {
        // Copy offspring to parents:
        for (ii = 0; ii < tercio12; ++ii) {
          bichoP[ii] = bichoH[ii];
        }
      }

      // Inversion: swap two last adjacent Ne values with probability frecinversion
      if (uniforme01(genera) < frecinversion) {
        for (ii = 0; ii < tercio12; ++ii) {
          bicho = bichoP[ii];
          if (bicho.nseg >= 4) {
            posiblock = bicho.nseg - 2;
            aa = bicho.Nebl[posiblock + 1];
            flag = true;
            bb = aa / bicho.Nebl[posiblock - 1];
            if (bb > invtopesalto && bb < topesalto) {
              bicho.Nebl[posiblock + 1] = bicho.Nebl[posiblock];
              bicho.Nebl[posiblock] = aa;
              SC = CalculaSC(bicho, rohprd);
              bicho.SCval = SC; // ++
              bichoP[ii] = bicho;
            }
          }
        }
      }

      // Sort parents
      for (ii = 0; ii < tercio12 - 1; ++ii) {
        indexminSC = ii;
        for (jj = ii + 1; jj < tercio12; ++jj) {
          if (bichoP[jj].SCval < bichoP[indexminSC].SCval) {
            indexminSC = jj;
          }
        }
        if (indexminSC != ii) { // swap
          bicho = bichoP[ii];
          bichoP[ii] = bichoP[indexminSC];
          bichoP[indexminSC] = bicho;
        }
      }

      SCmed = nsegmed = 0;
      for (ii = 0; ii < tercio12; ++ii) {
        SCmed += bichoP[ii].SCval;
        nsegmed += bichoP[ii].nseg;
      }
      SCmed /= tercio12;
      nsegmed /= tercio12;
      SCbest = bichoP[0].SCval;
      nsegbest = bichoP[0].nseg;
      // if (flagdebug) {
      //  salida.open(fichero_sal_evol, std::ios::app);
      //  salida << gen << "\t" << SCbest << "\t" << SCmed << "\t" << nsegbest
      //         << "\t" << nsegmed << "\n";
      //  salida.close();
      // }
      if (gen % 50 == 0 && gen > 0) {
        omp_set_lock(&lock);
        progressCounter += 50;
        params.progress.SetTaskProgress(progressCounter);
        // params.progress.PrintProgress();
        omp_unset_lock(&lock);
      }
    }
    // if (flagdebug) {
    //  salida.open(fichero_caca, std::ios::app);
    //  for (j = 0; j < bichoP[0].nseg; ++j) {
    //    salida << j << "\t" << bichoP[0].segbl[j] << "\t"
    //           << bichoP[0].Nebl[j] / 2.0 << "\n";
    //  }
    //  salida.close();
    // }

    bicho = bichoP[0]; // Best individual
    SC = CalculaSC(bicho, rohprd);       // This loads Neblock, segdesdehasta, SC, ...

    for (i = 0; i <= bicho.nseg; ++i) {
      Neblock[i] = bicho.Nebl[i];
      segdesdehasta[i] = bicho.segbl[i];
    }

    conta = 0;
    for (i = 0; i < bicho.nseg; ++i) {
      for (j = segdesdehasta[i]; j < segdesdehasta[i + 1]; ++j) {
        Ne[conta] = Neblock[i];
        conta = conta + 1;
      }
    }
    conta2 = conta;
    // Geometric mean of the top 'muestrasalida' solutions
    for (i = 0; i < conta2; ++i) {
      sumNe[i] = 1;
    } //
    for (ii = 0; ii < muestrasalida; ++ii) {
      conta = 0;
      for (i = 0; i < bichoP[ii].nseg; ++i) {
        for (j = bichoP[ii].segbl[i]; j < bichoP[ii].segbl[i + 1]; ++j) {
          sumNe[conta] *= pow(bichoP[ii].Nebl[i], 1.0 / muestrasalida);
          conta = conta + 1;
          if (conta > gmax) {
            conta = gmax;
            break;
          }
        }
        if (conta > gmax) {
          conta = gmax;
          break;
        }
      }
    }

    tfinal = time(NULL);
    tpas = (clock() - tiniabs);
    double mu1=1+mu;
    for (i = 0; i < nlin; ++i) {
      ll = xrohval[i] / 100;
      //            edadrohs[i]=3/(2*ll);
      Nerohs[i] = wl / (2 * xrohobs[i] * ll * ll);
      //            BB=ll*(1-d/H)-2*d/H  +  (ll*(1-d/H)+2*d/H)*exp(-ll*(H/d-1));
      BB = ll;
      s = ll / (d / H);
      // TODO(Enrique): this if-block is unused
      // if (s <= 2) {
      //  sk = 7.806;
      //  sm = 0.878;
      // } else if (s <= 3.8) {
      //  sk = 3.62;
      //  sm = 1.3212;
      // } else {
      //  sk = 1.108;
      //  sm = 2.234;
      // }
      //            QQ=exp(-s/sm)*sk;
      QQ = 0;
      double newNe;
      for (j = 1; j < 1000; ++j) {
        newNe = (4 * wl * BB * mu1*mu1) / (xrohobs[i] * pow((2 + QQ) * ll * mu1 + 1 / (2 * Nerohs[i]) - 4 * d / H, 3));
        // newNe = (4 * wl * BB ) / (xrohobs[i] * pow((2 + QQ) * ll + 1 / (2 * Nerohs[i]) - 4 * d / H, 3));
        Nerohs[i] = 0.1 * newNe + 0.9 * Nerohs[i];
        if (abs(newNe-Nerohs[i])/newNe < 0.0001){break;}
      }
    }
    // std::cout<< std::fixed << std::setprecision(8)<< "  d:"<<d<<" H:"<<H<<std::endl;
    for (ii = 0; ii < conta; ++ii) {
      t_acusumNe[tid][ii] += sumNe[ii];
    }
    for (ii = 0; ii < nlin; ++ii) {
      t_acurohprd[tid][ii] += rohprd[ii];
    }
    tcounter[tid]++;
  }
  // for (int tidx = 0; tidx < numThreads; ++tidx) {
  for (int tidx = 0; tidx < NREPETICIONES; ++tidx) {
    if (tcounter[tidx] > 0) {
      for (ii = 0; ii < ngenmax; ++ii) {
        // t_acusumNe[tidx][ii] /= static_cast<double>(tcounter[tidx]);
        acusumNe[ii] += t_acusumNe[tidx][ii];
      }

      for (ii = 0; ii < nlinmax; ++ii) {
        // t_acurohprd[tidx][ii] /= static_cast<double>(tcounter[tidx]);
        acurohprd[ii] += t_acurohprd[tidx][ii];
      }
    }
  }

  for (int zz = 0; zz < ngenmax; ++zz) {
    // acusumNe[zz] /= static_cast<double>(numThreads);
    acusumNe[zz] /= static_cast<double>(NREPETICIONES);
  }
  for (int zz = 0; zz < nlinmax; ++zz) {
    // acurohprd[zz] /= static_cast<double>(numThreads);
    acurohprd[zz] /= static_cast<double>(NREPETICIONES);
  }

  int maxedad;
  int minedad = floor(edadrohs[nlin - 1] + 1);
  // int minedad = floor((edadrohs[nlin - 1] + edadrohs[nlin - 2]) / 2);
  // minedad=std::max(minedad,0);
  double cortecMmax=dH5_total * 4.5;
  cortecMmax = std::min(cortecMmax, 10.0);
  for (i = 0; i < nlin; ++i) {
    if (Nerohs[i]>1){
      if (cortecMmax<rohval[i]){
        if (i==0){
          maxedad=edadrohs[i];
          break;
        }
        else{
          double factor=(rohval[i]-cortecMmax)/(rohval[i]-rohval[i-1]);
          maxedad=int(edadrohs[i]+(edadrohs[i-1]-edadrohs[i])*factor);
          break;
        }
      }
    }
  }
  maxedad = std::min(maxedad, ngenmax);

  salida.open(fichero_sal_NeH, std::ios::out);
  if ((params.diphap == 1) && (!params.flagX)) {
    salida << "Generation"
           << "\t"
           << "Ne_hap"
           << "\n";
  } else {
    salida << "Generation"
           << "\t"
           << "Ne_dip"
           << "\n";
  }
  for (j = minedad; j < maxedad; ++j) {
      if ((params.diphap == 1) && (!params.flagX)) {
        salida << j << "\t" << acusumNe[j] * 2 << "\n";
      } else {
        salida << j << "\t" << acusumNe[j] << "\n";
      }
      // if (j == minedad) {
      //   salida << "#^^^^^-low_reliability-^^^^^\n";
      // }
  }
  salida.close();

  salida.open(fichero_sal_roh, std::ios::out);
  if (params.diphap == 1) {
    salida
        << "Num_of_segments\tBin_class_in_cM\tObserved_Prop_ofGenome_IBD_due_"
           "to_this_class\tExpected_Prop_ofGenome_IBD_due_to_this_class\tNe_"
           "hap_for_this_class\tAverage_coal_time_of_this_class\n";
    for (i = 0; i < nlin; ++i) {
      if (Nerohs[i]>1){
        salida << nbin[i] << "\t" << rohval[i] << "\t" << rohobs[i] << "\t"
             << acurohprd[i] << "\t" << Nerohs[i] * 2 << "\t" << edadrohs[i]
             << "\n";
      }
    }
  } else {
    salida
        << "Num_of_segments\tBin_class_in_cM\tObserved_Prop_ofGenome_IBD_due_"
           "to_this_class\tExpected_Prop_ofGenome_IBD_due_to_this_class\tNe_"
           "dip_for_this_class\tAverage_coal_time_of_this_class\n";
    for (i = 0; i < nlin; ++i) {
      if (Nerohs[i]>1){
        salida << nbin[i] << "\t" << rohval[i] << "\t" << rohobs[i] << "\t"
              << acurohprd[i] << "\t" << Nerohs[i] << "\t" << edadrohs[i]
              << "\n";
      }
    }
  }
  salida.close();

  tfinal = omp_get_wtime() - tinicio;
  // TODO(Enrique): esto no se usa
  tpas = (clock() - tiniabs);

  // std::cout << " End of processing: " << fichero << "\n";
  // std::cout << " Total run time in seconds: " << tfinal << "\n";

  return 0;
}

// Core of computation:
// This subroutine calculates the values of SD2_c, SW_c and the predicted d2_c
// for a particular set of Ne_t values stored in Ne[].

double CalculaSC(serie& bicho, double* rohprd) {
  // TODO(Enrique): exponente, hastasegmento,
  //                sumC, Hd, d2, d1 no se usan
  int i, t, ii, jj, hastat;
  double mu1, Hd4, Hde, factorC, factorX, sumCX, e2lHd4, corr2, l;
  double A, N, N21, N2, C, KK, addsumCX, Q, s, sk, sm, eQ[nlinmax];
  double Net[ngenmax] = {0};
  bool negativos;

  nsegmentos = bicho.nseg;
  for (i = 0; i <= nsegmentos; ++i) {
    Neblock[i] = bicho.Nebl[i];
    segdesdehasta[i] = bicho.segbl[i];
  }
  Necons = Neblock[nsegmentos - 1];
  hastat = segdesdehasta[nsegmentos - 1];
  for (i = 0; i < nsegmentos; ++i) {
    for (ii = segdesdehasta[i]; ii < segdesdehasta[i + 1]; ++ii) {
      Net[ii] = Neblock[i];
    }
  }
  mu1 = 1 + mu;
  SC = 0;
  for (ii = 0; ii < nlin; ++ii) {
    rohprd[ii] = 0;
  }
  for (jj = 0; jj < 5; jj++) {
    Hd4 = pow((1 - dH_vals[jj] * mu) / (1 - dH_vals[jj] * mu1), 4);
    for (ii = 0; ii < nlin; ++ii) {
      l = (rohval[ii]) / 100;
      s = l / dH_vals[jj];
      if (s <= 2) {
        sk = 7.806;
        sm = 0.878;
      } else if (s <= 3.8) {
        sk = 3.62;
        sm = 1.3212;
      } else {
        sk = 1.108;
        sm = 2.234;
      }
      Q = exp(-l / dH_vals[jj] / sm) * sk;
      // Q=0.1; // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      eQ[ii] = exp(-2 * l * mu1 - Q * l);
    }
    for (ii = 0; ii < nlin; ++ii) {
      l = (rohval[ii]) / 100;
      // Hde=exp(-2*l*mu1)*Hd4;
      Hde = eQ[ii] * Hd4;
      factorC = 1;
      // TODO(Enrique): factorX no se usa
      factorX = 1;
      sumCX = 0;
      e2lHd4 = 1;
      corr2 = 4 * exp(-2 * l * err) *
              ((l * mu1 * (1 - dH_vals[jj]) - 2 * dH_vals[jj]) * mu1 +
               (l * (1 - dH_vals[jj]) + 2 * dH_vals[jj] * mu1) *
                   exp(-l * (1 / dH_vals[jj] - 1)));
      double tt = 1;
      for (t = 1; t < hastat; ++t) {
        N2 = 0.5 / Net[t];
        N21 = 1.0 - N2;
        e2lHd4 *= Hde;
        C = factorC * N2;
        factorX = corr2 * tt * tt * e2lHd4;
        sumCX += factorX * C;
        factorC *= N21;
        tt++;
      }
      N = Necons;
      N2 = 0.5 / N;
      N21 = 1.0 - N2;
      KK = Hde * N21;
      addsumCX = corr2 * factorC * e2lHd4 * N2;
      addsumCX *=
          ((hastat - 1) * (hastat - 1) * KK * KK -
           (2 * hastat * hastat - 2 * hastat - 1) * KK + hastat * hastat) /
          pow(1 - KK, 3);
      sumCX += addsumCX;

      rohprd[ii] += sumCX * wl;
    }
  }
  negativos=false;
  for (ii = 0; ii < nlin; ++ii) {
    rohprd[ii] /= 5.0;
    if (rohprd[ii] <0){
      negativos=true;
    }
    A = log10(rohobs[ii]) - log10(rohprd[ii]);
    //  A=(rohobs[ii])-(rohprd[ii]);
    SC += (A *= A);
  }
  // if (negativos){SC=std::numeric_limits<double>::max();}
  return SC;
}
