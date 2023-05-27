#include "convenience.h"
#include "cell.h"
#include "tsc_block.h"

using namespace std;

string help_message()
{
  std::string t = "\n----------------------------------------------------------------------------\n";
  t.append("          These commands and arguments are known by NoSpherA2:\n");
  t.append("----------------------------------------------------------------------------\n\n");
  t.append("   -wfn            <FILENAME>.xxx         Read the following wavefunction file.\n");
  t.append("                                          Supported filetypes: .wfn/wfx/ffn; .molden; .xyz; .gbw; fch* (UNTESTED!)\n");
  t.append("   -fchk           <FILENAME>.fchk        Write a wavefunction to the given filename\n");
  t.append("   -b              <FILENAME>             Read this basis set\n");
  t.append("   -d              <PATH>                 Path to basis_sets directory with basis_sets in tonto style\n");
  t.append("   --help/-help/--h                       print this help\n");
  t.append("   -v                                     Turn on Verbose (debug) Mode (Slow and a LOT of output!)\n");
  t.append("   -v2                                    Even more stuff\n");
  t.append("   -mult           <NUMBER>               Input multiplicity of wavefunction (otherwise attempted to be read from the wfn)\n");
  t.append("   -method         <METHOD NAME>          Can be RKS or RHF to distinguish between DFT and HF\n");
  t.append("   -cif            <FILENAME>.cif         CIF to get labels of atoms to use for calculation of scatteriung factors\n");
  t.append("   -IAM                                   Make scattering factors based on Thakkar functions for atoms in CIF\n");
  t.append("   -xyz            <FILENAME>.xyz         Read atom positions from this xyz file for IAM\n");
  t.append("   -hkl            <FILENAME>.hkl         hkl file (ideally merged) to use for calculation of form factors.\n");
  t.append("   -group          <LIST OF INT NUMBERS>  Disorder groups to be read from the CIF for consideration as asym unit atoms (space separated).\n");
  t.append("   -acc            0,1,2,3,4...           Accuracy of numerical grids used, where the bumber indicates a pre-defined level. 4 should be considered maximum,\n");
  t.append("                                          anything above will most likely introduce numberical error and is just implemented for testing purposes.");
  t.append("   -gbw2wfn                               Only reads wavefucntion from .gbw specified by -wfn and prints it into .wfn format.\n");
  t.append("   -tscb           <FILENAME>.tsb         Convert binary tsc file to bigger, less accurate human-readable form.\n");
  t.append("   -twin     3x3 floating-point-matrix in the form -1 0 0 0 -1 0 0 0 -1 which contains the twin matrix to use.\n");
  t.append("             If there is more than a single twin law to be used, use the twin command multiple times.\n");
  t.append("   -merge          <List of .tsc files>   Names/Paths to .tsc files to be merged.\n");
  t.append("   -merge_nocheck  <List of .tsc files>   Names/Paths to .tsc files to be merged. They need to have identical hkl values.\n");
  t.append("   -mtc            <List of .wfns + parts>  Performs calculation for a list of wavefunctions (=Multi-Tsc-Calc), where asymmetric unit is.\n");
  t.append("                                            taken from given CIF. Also disorder groups are required per file as comma separated list\n");
  t.append("                                            without spaces.\n   Typical use Examples:\n");
  t.append("      Normal:       NoSpherA2.exe -cif A.cif -hkl A.hkl -wfn A.wfx -acc 1 -cpus 7\n");
  t.append("      thakkar-tsc:  NoSpherA2.exe -cif A.cif -hkl A.hkl -xyz A.xyz -acc 1 -cpus 7 -IAM\n");
  t.append("      Disorder:     NoSpherA2.exe -cif A.cif -hkl A.hkl -acc 1 -cpus 7 -mtc 1.wfn 0,1 2.wfn 0,2 3.wfn 0,3\n");
  t.append("      fragHAR:      NoSpherA2.exe -cif A.cif -hkl A.hkl -acc 1 -cpus 7 -cmtc 1.wfn 1.cif 0 2.wfn 2.cif 0 3_1.wfn 3_1.cif 0,1 3_2.wfn 3_2.cif 0,2\n");
  t.append("      merging tscs: NoSpherA2.exe -merge A.tsc B.tsc C.tsc\n");
  t.append("      merge tsc(2): NoSpherA2.exe -merge_nocheck A.tsc B.tsc C.tsc  (MAKE SURE THEY HAVE IDENTICAL HKL INIDCES!!)\n");
  return t;
}
string NoSpherA2_message()
{
  string t = "    _   __     _____       __              ___   ___\n";
  t.append("   / | / /___ / ___/____  / /_  ___  _____/   | |__ \\\n");
  t.append("  /  |/ / __ \\\\__ \\/ __ \\/ __ \\/ _ \\/ ___/ /| | __/ /\n");
  t.append(" / /|  / /_/ /__/ / /_/ / / / /  __/ /  / ___ |/ __/\n");
  t.append("/_/ |_/\\____/____/ .___/_/ /_/\\___/_/  /_/  |_/____/\n");
  t.append("                /_/\n");
  t.append("This software is part of the cuQCT software suite developed by Florian Kleemiss.\n");
  t.append("Please give credit and cite corresponding pieces!\n");
  t.append("NoSpherA2 was published at : Kleemiss et al. Chem.Sci., 2021, 12, 1675 - 1692\n");
  return t;
}

string build_date() {
  return ("This Executable was built on: " + string(__DATE__) + " " + string(__TIME__) + "\n");
}



bool is_similar_rel(const double& first, const double& second, const double& tolerance)
{
  double diff = abs(first - second);
  if (diff > abs((first + second + 0.01) * tolerance / 2))
    return false;
  else
    return true;
};

bool is_similar(const double& first, const double& second, const double& tolerance)
{
  double diff = abs(first - second);
  if (diff > pow(10, tolerance))
    return false;
  else
    return true;
};

bool is_similar_abs(const double& first, const double& second, const double& tolerance)
{
  double diff = abs(first - second);
  if (diff > abs(tolerance))
    return false;
  else
    return true;
};

void cls()
{
  //    cout << string( 100, '\n' );
#ifdef _WIN32
  if (system("CLS")) cout << "this should not happen...!" << endl;
#else
  if (system("clear")) cout << "this should not happen...!" << endl;
#endif
};

string atnr2letter(const int& nr)
{
  vector <string> Labels{ "DM","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca"
    ,"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"
    ,"Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe"
    ,"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn"
    ,"Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr"
  };
  if (nr == 0) {
    //Exception for Q peaks in residual maps
    return "Q";
  }
  if (nr > 103 || nr < 0) {
    if (nr == 119) {
      //Exception for Q in ECPs from ORCA
      return "Q";
    }
    cout << "Only yet implemented from H-Lr, ask Florian for improvements or give a reasonable number between 1-103!" << endl;
    return ("PROBLEM");
  }
  else return Labels[nr];
};

int get_Z_from_label(const char* tmp)
{
  if (strcmp(tmp, "H") == 0)  return 0;
  else  if (strcmp(tmp, "D") == 0)   return 0;
  else  if (strcmp(tmp, "T") == 0)   return 0;
  else 	if (strcmp(tmp, "He") == 0)  return 1;
  else 	if (strcmp(tmp, "Li") == 0)  return 2;
  else 	if (strcmp(tmp, "Be") == 0)  return 3;
  else 	if (strcmp(tmp, "B") == 0)   return 4;
  else 	if (strcmp(tmp, "C") == 0)   return 5;
  else 	if (strcmp(tmp, "N") == 0)   return 6;
  else 	if (strcmp(tmp, "O") == 0)   return 7;
  else 	if (strcmp(tmp, "F") == 0)   return 8;
  else 	if (strcmp(tmp, "Ne") == 0)  return 9;
  else 	if (strcmp(tmp, "Na") == 0)  return 10;
  else 	if (strcmp(tmp, "Mg") == 0)  return 11;
  else 	if (strcmp(tmp, "Al") == 0)  return 12;
  else 	if (strcmp(tmp, "Si") == 0)  return 13;
  else 	if (strcmp(tmp, "P") == 0)   return 14;
  else 	if (strcmp(tmp, "S") == 0)   return 15;
  else 	if (strcmp(tmp, "Cl") == 0)  return 16;
  else 	if (strcmp(tmp, "Ar") == 0)  return 17;
  else 	if (strcmp(tmp, "K") == 0)   return 18;
  else 	if (strcmp(tmp, "Ca") == 0)  return 19;
  else 	if (strcmp(tmp, "Sc") == 0)  return 20;
  else 	if (strcmp(tmp, "Ti") == 0)  return 21;
  else 	if (strcmp(tmp, "V") == 0)   return 22;
  else 	if (strcmp(tmp, "Cr") == 0)  return 23;
  else 	if (strcmp(tmp, "Mn") == 0)  return 24;
  else 	if (strcmp(tmp, "Fe") == 0)  return 25;
  else 	if (strcmp(tmp, "Co") == 0)  return 26;
  else 	if (strcmp(tmp, "Ni") == 0)  return 27;
  else 	if (strcmp(tmp, "Cu") == 0)  return 28;
  else 	if (strcmp(tmp, "Zn") == 0)  return 29;
  else 	if (strcmp(tmp, "Ga") == 0)  return 30;
  else 	if (strcmp(tmp, "Ge") == 0)  return 31;
  else 	if (strcmp(tmp, "As") == 0)  return 32;
  else 	if (strcmp(tmp, "Se") == 0)  return 33;
  else 	if (strcmp(tmp, "Br") == 0)  return 34;
  else 	if (strcmp(tmp, "Kr") == 0)  return 35;
  else 	if (strcmp(tmp, "Rb") == 0)  return 36;
  else 	if (strcmp(tmp, "Sr") == 0)  return 37;
  else 	if (strcmp(tmp, "Y") == 0)   return 38;
  else 	if (strcmp(tmp, "Zr") == 0)  return 39;
  else 	if (strcmp(tmp, "Nb") == 0)  return 40;
  else 	if (strcmp(tmp, "Mo") == 0)  return 41;
  else 	if (strcmp(tmp, "Tc") == 0)  return 42;
  else 	if (strcmp(tmp, "Ru") == 0)  return 43;
  else 	if (strcmp(tmp, "Rh") == 0)  return 44;
  else 	if (strcmp(tmp, "Pd") == 0)  return 45;
  else 	if (strcmp(tmp, "Ag") == 0)  return 46;
  else 	if (strcmp(tmp, "Cd") == 0)  return 47;
  else 	if (strcmp(tmp, "In") == 0)  return 48;
  else 	if (strcmp(tmp, "Sn") == 0)  return 49;
  else 	if (strcmp(tmp, "Sb") == 0)  return 50;
  else 	if (strcmp(tmp, "Te") == 0)  return 51;
  else 	if (strcmp(tmp, "I") == 0)   return 52;
  else 	if (strcmp(tmp, "Xe") == 0)  return 53;
  else 	if (strcmp(tmp, "Cs") == 0)  return 54;
  else 	if (strcmp(tmp, "Ba") == 0)  return 55;
  else 	if (strcmp(tmp, "La") == 0)  return 56;
  else 	if (strcmp(tmp, "Ce") == 0)  return 57;
  else 	if (strcmp(tmp, "Pr") == 0)  return 58;
  else 	if (strcmp(tmp, "Nd") == 0)  return 59;
  else 	if (strcmp(tmp, "Pm") == 0)  return 60;
  else 	if (strcmp(tmp, "Sm") == 0)  return 61;
  else 	if (strcmp(tmp, "Eu") == 0)  return 62;
  else 	if (strcmp(tmp, "Gd") == 0)  return 63;
  else 	if (strcmp(tmp, "Tb") == 0)  return 64;
  else 	if (strcmp(tmp, "Dy") == 0)  return 65;
  else 	if (strcmp(tmp, "Ho") == 0)  return 66;
  else 	if (strcmp(tmp, "Er") == 0)  return 67;
  else 	if (strcmp(tmp, "Tm") == 0)  return 68;
  else 	if (strcmp(tmp, "Yb") == 0)  return 69;
  else 	if (strcmp(tmp, "Lu") == 0)  return 70;
  else 	if (strcmp(tmp, "Hf") == 0)  return 71;
  else 	if (strcmp(tmp, "Ta") == 0)  return 72;
  else 	if (strcmp(tmp, "W") == 0)   return 73;
  else 	if (strcmp(tmp, "Re") == 0)  return 74;
  else 	if (strcmp(tmp, "Os") == 0)  return 75;
  else 	if (strcmp(tmp, "Ir") == 0)  return 76;
  else 	if (strcmp(tmp, "Pt") == 0)  return 77;
  else 	if (strcmp(tmp, "Au") == 0)  return 78;
  else 	if (strcmp(tmp, "Hg") == 0)  return 79;
  else 	if (strcmp(tmp, "Ti") == 0)  return 80;
  else 	if (strcmp(tmp, "Pb") == 0)  return 81;
  else 	if (strcmp(tmp, "Bi") == 0)  return 82;
  else 	if (strcmp(tmp, "Po") == 0)  return 83;
  else 	if (strcmp(tmp, "At") == 0)  return 84;
  else 	if (strcmp(tmp, "Rn") == 0)  return 85;

  else 	if (strcmp(tmp, "Fr") == 0)  return 86;
  else 	if (strcmp(tmp, "Ra") == 0)  return 87;
  else 	if (strcmp(tmp, "Ac") == 0)  return 88;
  else 	if (strcmp(tmp, "Th") == 0)  return 89;
  else 	if (strcmp(tmp, "Pa") == 0)  return 90;
  else 	if (strcmp(tmp, "U") == 0)   return 91;
  else 	if (strcmp(tmp, "Np") == 0)  return 92;
  else 	if (strcmp(tmp, "Pu") == 0)  return 93;
  else 	if (strcmp(tmp, "Am") == 0)  return 94;
  else 	if (strcmp(tmp, "Cm") == 0)  return 95;
  else 	if (strcmp(tmp, "Bk") == 0)  return 96;
  else 	if (strcmp(tmp, "Cf") == 0)  return 97;
  else 	if (strcmp(tmp, "Es") == 0)  return 98;
  else 	if (strcmp(tmp, "Fm") == 0)  return 99;
  else 	if (strcmp(tmp, "Md") == 0)  return 100;
  else 	if (strcmp(tmp, "No") == 0)  return 101;
  else 	if (strcmp(tmp, "Lr") == 0)  return 102;
  else 	if (strcmp(tmp, "Rf") == 0)  return 103;
  else 	if (strcmp(tmp, "Db") == 0)  return 104;
  else 	if (strcmp(tmp, "Sg") == 0)  return 105;
  else 	if (strcmp(tmp, "Bh") == 0)  return 106;
  else 	if (strcmp(tmp, "Hs") == 0)  return 107;
  else 	if (strcmp(tmp, "Mt") == 0)  return 108;
  else 	if (strcmp(tmp, "Ds") == 0)  return 109;
  else 	if (strcmp(tmp, "Rg") == 0)  return 110;
  else                               return -1;
}
/*
cosinus_annaeherung::cosinus_annaeherung() : mSize(0), mBase_values(nullptr), mStepwidth(1.0) {
  resize(100);
}

void cosinus_annaeherung::resize(size_t size)
{
  mSize = size;
  if (mBase_values) delete[] mBase_values;
  mBase_values = new double[mSize + 1];
#pragma omp parallel for
  for (auto i = 0; i < mSize + 1; i++)  // Fuer einen Werte mehr die Stueststellen speichern
  {
    double y = cos((MPI2 * i) / mSize);
    // cout << "resize: i="<<i<<" y=" << y << endl;
    mBase_values[i] = y;
  }
  mStepwidth = MPI2 / size;
}

double cosinus_annaeherung::calculate_error_at(double x) const
{
  return cos(x) - get(x);
}
*/
void copy_file(string& from, string& to)
{
  ifstream source(from.c_str(), ios::binary);
  ofstream dest(to.c_str(), ios::binary);

  dest << source.rdbuf();

  source.close();
  dest.close();
};

//---------------------------Configuration files ---------------------------------------------------

string get_home_path(void)
{
#ifdef _WIN32
  string temp1 = getenv("HOMEDRIVE");
  string temp2 = getenv("HOMEPATH");
  temp1.append(temp2);
  return temp1;
#else
  string home = getenv("HOME");
  return home;
#endif
}

void join_path(string& s1, string& s2)
{
#ifdef _WIN32
  s1.append("\\");
#else
  if (s1.substr(s1.length() - 1) == "/")
    s1.append("/");
#endif
  s1.append(s2);
}

string get_filename_from_path(const string& input)
{
#ifdef _WIN32
  return input.substr(input.rfind("\\") + 1);
#else
  return input.substr(input.rfind("/") + 1);
#endif
}

string get_foldername_from_path(const string& input)
{
#ifdef _WIN32
  return input.substr(0, input.rfind("\\") + 1);
#else
  return input.substr(0, input.rfind("/") + 1);
#endif
}

string get_basename_without_ending(const string& input)
{
  return input.substr(0, input.rfind("."));
}

void write_template_confi()
{
  string line;
  string programs = get_home_path();
  string filename = ".cuQCT.conf";
  join_path(programs, filename);
  if (exists(programs)) {
    cout << "File already exists! Aborting!" << endl;
    return;
  }
  ofstream conf(programs.c_str());
#ifdef _WIN32
  conf << "gaussian=\"D:\\g09\\g09\\\"" << endl;
  conf << "turbomole=\"D:\\turbomole\\dscf7.1\\\"" << endl;
  conf << "basis=\"D:\\tonto\\basis_sets\\\"" << endl;
#else
  conf << "gaussian=\"/usr/local/g09/g09\"" << endl;
  conf << "turbomole=\"/usr/local/bin/dscf7.1\"" << endl;
  conf << "basis=\"/basis_sets/\"" << endl;
#endif
  conf << "cpu=4" << endl;
  conf << "mem=4.0" << endl;
  conf << "rho=1" << endl;
  conf << "rdg=1" << endl;
  conf << "eli=0" << endl;
  conf << "elf=0" << endl;
  conf << "lap=0" << endl;
  conf << "esp=0" << endl;
  conf << "efv=0" << endl;
  conf << "def=0" << endl;
  conf << "hir=0" << endl;
  conf.flush();
  conf.close();
#ifdef _WIN32
  //	const wchar_t* fileLPCWSTR = programs.c_str();
  //	wstring stemp = wstring(programs.begin(), programs.end());
  //	int attr = GetFileAttributes(stemp.c_str());
  //	if ((attr & FILE_ATTRIBUTE_HIDDEN) == 0) {
  //		SetFileAttributes(stemp.c_str(), attr | FILE_ATTRIBUTE_HIDDEN);
  //	}
#endif
  return;
};

int program_confi(string& gaussian_path, string& turbomole_path, string& basis, int& ncpus, double& mem, bool debug, bool expert, unsigned int counter)
{
  counter++;
  if (counter == 3) {
    cout << "Too many iterations of tries to read config file, better abort..." << endl;
    return -1;
  }
  string programs = get_home_path();
  string filename = ".cuQCT.conf";
  join_path(programs, filename);
  ifstream conf(programs.c_str());
  if (debug) cout << programs << endl;
  string line;
  if (conf.good()) {
    if (debug) cout << "File is valid, continuing..." << endl;
  }
  else {
    if (expert) {
      cout << "couldn't open or find .cuQCT.conf, in your home folder: " << programs << ", writing a template for you!" << endl;
      write_template_confi();
      if (program_confi(gaussian_path, turbomole_path, basis, ncpus, mem, debug, expert, counter) != 1) return -1;
      cout << "Wrote a template for you, read default values!" << endl;
      return 0;
    }
    else {
      write_template_confi();
      program_confi(gaussian_path, turbomole_path, basis, ncpus, mem, debug);
      return 0;
    }
  }
  conf.seekg(0);
  getline(conf, line);
  size_t length;
  char* tempchar = new char[200];
  int run = 0;
  while (!conf.eof()) {
    switch (run) {
    case 0:
      length = line.copy(tempchar, line.size() - 11, 10);
      tempchar[length] = '\0';
      gaussian_path = tempchar;
      break;
    case 1:
      length = line.copy(tempchar, line.size() - 12, 11);
      tempchar[length] = '\0';
      turbomole_path = tempchar;
      break;
    case 2:
      length = line.copy(tempchar, line.size() - 8, 7);
      tempchar[length] = '\0';
      basis = tempchar;
      break;
    case 3: length = line.copy(tempchar, line.size() - 3, 4);
      tempchar[length] = '\0';
      ncpus = stoi(tempchar);
      break;
    case 4: length = line.copy(tempchar, line.size() - 3, 4);
      tempchar[length] = '\0';
      mem = stod(tempchar);
      break;
    default:
      if (debug) cout << "found everything i was looking for, if you miss something check the switch" << endl;
      break;
    }
    if (debug) cout << run << ". line: " << tempchar << endl;
    run++;
    getline(conf, line);
  }
  return 1;
};

int filetype_identifier(string& file, bool debug)
{
  /*
  List of filetypes and correpsonding values:
          -1: unreadable keyword
  -i/o *.wfn 		2: wfn
  -i/o *.ffn 		4: ffn
  -i *.out 		1: crystal output
  -c/o *.cub(e) 	3: cube file
  -g/o *.grd      6: XDGraph grid file
  -o *.(F)fc(C)hk 5: fchk
  */
  if (debug) {
    cout << "Testing WFN:  " << file.find(".wfn") << endl
      << "Testing out:  " << file.find(".out") << endl
      << "Testing FFN:  " << file.find(".ffn") << endl
      << "Testing CUB:  " << file.find(".cub") << endl
      << "Testing CUBE: " << file.find(".cube") << endl
      << "Testing Grid: " << file.find(".grd") << endl
      << "Testing fchk: " << file.find(".fchk") << endl
      << "Testing FChk: " << file.find(".FChk") << endl
      << "Testing Fchk: " << file.find(".Fchk") << endl;
    cout << "string::npos: " << string::npos << endl;
  }
  int temp_type = 0;
  size_t found, temp;
  temp = 0;
  if (debug) cout << "Temp before any checks: " << temp << endl;
  vector <string> types{ ".out",".wfn",".ffn",".cub",".cube",".grd",".fchk",".Fchk",".FChk" };
  if (file.find(".wfn") != string::npos) {
    if (debug) cout << "Checking for" << ".wfn" << endl;
    temp_type = 2;
    found = file.rfind(".wfn");
    if (debug) cout << "Found: " << found << endl;
    for (int i = 0; i < types.size(); i++)
      if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
        temp = file.rfind(types[i]);
    if (debug) cout << "Temp: " << temp << endl;
    if (temp == found) return temp_type;
    else {
      temp = 0;
      if (debug) cout << "Moving on!" << endl;
    }
  }
  if (file.find(".out") != string::npos) {
    if (debug) cout << "Checking for" << ".out" << endl;
    temp_type = 1;
    found = file.rfind(".out");
    for (int i = 0; i < types.size(); i++)
      if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
        temp = file.rfind(types[i]);
    if (temp == found) return temp_type;
    else temp = 0;
  }
  if (file.find(".ffn") != string::npos) {
    if (debug) cout << "Checking for" << ".ffn" << endl;
    temp_type = 4;
    found = file.rfind(".ffn");
    for (int i = 0; i < types.size(); i++)
      if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
        temp = file.rfind(types[i]);
    if (temp == found) return temp_type;
    else temp = 0;
  }
  if (file.find(".cub") != string::npos) {
    if (debug) cout << "Checking for" << ".cub" << endl;
    temp_type = 3;
    found = file.rfind(".cub");
    for (int i = 0; i < types.size(); i++)
      if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
        temp = file.rfind(types[i]);
    if (temp == found) return temp_type;
    else {
      temp = 0;
      if (debug) cout << "Moving on!" << endl;
    }
  }
  if (file.find(".cube") != string::npos) {
    temp_type = 3;
    found = file.rfind(".cube");
    for (int i = 0; i < types.size(); i++)
      if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
        temp = file.rfind(types[i]);
    if (temp == found) return temp_type;
    else temp = 0;
  }
  if (file.find(".grd") != string::npos) {
    temp_type = 6;
    found = file.rfind(".grd");
    for (int i = 0; i < types.size(); i++)
      if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
        temp = file.rfind(types[i]);
    if (temp == found) return temp_type;
    else temp = 0;
  }
  if (file.find(".fchk") != string::npos) {
    temp_type = 5;
    found = file.rfind(".fchk");
    for (int i = 0; i < types.size(); i++)
      if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
        temp = file.rfind(types[i]);
    if (temp == found) return temp_type;
    else temp = 0;
  }
  if (file.find(".FChk") != string::npos) {
    temp_type = 5;
    found = file.rfind(".FChk");
    for (int i = 0; i < types.size(); i++)
      if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
        temp = file.rfind(types[i]);
    if (temp == found) return temp_type;
    else temp = 0;
  }
  if (file.find(".Fchk") != string::npos) {
    temp_type = 5;
    found = file.rfind(".Fchk");
    for (int i = 0; i < types.size(); i++)
      if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
        temp = file.rfind(types[i]);
    if (temp == found) return temp_type;
    else temp = 0;
  }
  return -1;
}

string go_get_string(ifstream& file, string search, bool rewind)
{
  if (rewind) {
    file.clear();
    file.seekg(0, file.beg);
  }
  string line;
  while (line.find(search) == string::npos && !file.eof() && getline(file, line))
    continue;
  if (file.eof())
    return "";
  else
    return line;
}

string shrink_string(string& input)
{
  while (input.find(" ") != -1) { input.erase(input.find(" "), 1); }
  while (input.find("1") != -1) { input.erase(input.find("1"), 1); }
  while (input.find("2") != -1) { input.erase(input.find("2"), 1); }
  while (input.find("3") != -1) { input.erase(input.find("3"), 1); }
  while (input.find("4") != -1) { input.erase(input.find("4"), 1); }
  while (input.find("5") != -1) { input.erase(input.find("5"), 1); }
  while (input.find("6") != -1) { input.erase(input.find("6"), 1); }
  while (input.find("7") != -1) { input.erase(input.find("7"), 1); }
  while (input.find("8") != -1) { input.erase(input.find("8"), 1); }
  while (input.find("9") != -1) { input.erase(input.find("9"), 1); }
  while (input.find("0") != -1) { input.erase(input.find("0"), 1); }
  while (input.find("(") != -1) { input.erase(input.find("("), 1); }
  while (input.find(")") != -1) { input.erase(input.find(")"), 1); }
  return input;
};

string shrink_string_to_atom(string& input, const int& atom_number)
{
  while (input.find(" ") != -1) { input.erase(input.find(" "), 1); }
  while (input.find("1") != -1) { input.erase(input.find("1"), 1); }
  while (input.find("2") != -1) { input.erase(input.find("2"), 1); }
  while (input.find("3") != -1) { input.erase(input.find("3"), 1); }
  while (input.find("4") != -1) { input.erase(input.find("4"), 1); }
  while (input.find("5") != -1) { input.erase(input.find("5"), 1); }
  while (input.find("6") != -1) { input.erase(input.find("6"), 1); }
  while (input.find("7") != -1) { input.erase(input.find("7"), 1); }
  while (input.find("8") != -1) { input.erase(input.find("8"), 1); }
  while (input.find("9") != -1) { input.erase(input.find("9"), 1); }
  while (input.find("0") != -1) { input.erase(input.find("0"), 1); }
  while (input.find("(") != -1) { input.erase(input.find("("), 1); }
  while (input.find(")") != -1) { input.erase(input.find(")"), 1); }
  string temp = atnr2letter(atom_number);
  err_checkf(temp != "PROBLEM", "Problem identifying atoms!", std::cout);
  if (input.find(temp) != 1) return temp;
  if (temp != "PROBLEM")
    while (input.size() > temp.size())
      input.pop_back();
  return input;
};

void progress_bar::write(double fraction)
{
  // clamp fraction to valid range [0,1]
  if (fraction < 0)
    fraction = 0;
  else if (fraction > 1)
    fraction = 1;

  auto width = bar_width - message.size();
  auto offset = bar_width - static_cast<unsigned>(width * round(fraction / precision) * precision);

  os << '\r' << message;
  os.write(full_bar.data() + offset, width);
  os << " [" << std::setw(3) << static_cast<int>(100 * round(fraction / precision) * precision) << "%] " << std::flush;
}

void readxyzMinMax_fromCIF(
  string cif,
  double* CoordMinMax,
  int* NbSteps,
  vector < vector < double > >& cm,
  double Resolution,
  ofstream& file,
  bool debug
)
{
  if (debug)
    file << "starting to read cif!" << endl;
  if (!exists(cif)) {
    file << "CIF does not exists!" << endl;
    return;
  }
  ifstream cif_input(cif.c_str(), ios::in);
  vector<bool> found;
  string line;
  found.resize(7);
  for (int k = 0; k < 7; k++)
    found[k] = false;
  double a = 0.0, b = 0.0, c = 0.0, v = 0.0;
  double alpha = 0.0, beta = 0.0, gamma = 0.0;
  vector <string> cell_keywords;
  cell_keywords.push_back("_cell_length_a");
  cell_keywords.push_back("_cell_length_b");
  cell_keywords.push_back("_cell_length_c");
  cell_keywords.push_back("_cell_angle_alpha");
  cell_keywords.push_back("_cell_angle_beta");
  cell_keywords.push_back("_cell_angle_gamma");
  cell_keywords.push_back("_cell_volume");
  if (debug)
    file << "Starting while !.eof()" << endl;
  while (!cif_input.eof()) {
    if (debug)
      file << "While line! " << setw(80) << line
      << setw(10) << a << found[0]
      << setw(10) << b << found[1]
      << setw(10) << c << found[2]
      << setw(10) << alpha << found[3]
      << setw(10) << beta << found[4]
      << setw(10) << gamma << found[5]
      << setw(10) << v << found[6] << endl;
    getline(cif_input, line);
    for (int k = 0; k < cell_keywords.size(); k++) {
      if (line.find(cell_keywords[k]) != string::npos) {
        switch (k) {
        case 0:
          a = stod(line.substr(cell_keywords[k].length(), line.find("(")));
          break;
        case 1:
          b = stod(line.substr(cell_keywords[k].length(), line.find("(")));
          break;
        case 2:
          c = stod(line.substr(cell_keywords[k].length(), line.find("(")));
          break;
        case 3:
          alpha = stod(line.substr(cell_keywords[k].length(), line.find("(")));
          break;
        case 4:
          beta = stod(line.substr(cell_keywords[k].length(), line.find("(")));
          break;
        case 5:
          gamma = stod(line.substr(cell_keywords[k].length(), line.find("(")));
          break;
        case 6:
          v = stod(line.substr(cell_keywords[k].length(), line.find("(")));
          break;
        default:
          file << "This is weird... should never get here... aborting!" << endl;
          return;
        }
        found[k] = true;
      }
    }
    if (found[0] == true && found[1] == true && found[2] == true && found[3] == true && found[4] == true && found[5] == true && found[6] == true)
      break;
  }
  double ca = cos(PI_180 * alpha);
  double cb = cos(PI_180 * beta);
  double cg = cos(PI_180 * gamma);
  double sa = sin(PI_180 * alpha);
  double sb = sin(PI_180 * beta);
  double sg = sin(PI_180 * gamma);
  double Vp = sqrt((1 - ca * ca - cb * cb - cg * cg) + 2 * (ca * cb * cg));
  double V = a * b * c * Vp;

  if (debug)
    file << "Making cm" << endl
    << a << " " << b << " " << c << " " << ca << " " << cb << " " << cg << " " << sa << " " << sb << " " << sg << " " << V << " vs. " << v << endl;

  cm[0][0] = a / 0.529177249;
  cm[1][0] = 0.0;
  cm[2][0] = 0.0;

  cm[0][1] = b * cg / 0.529177249;
  cm[1][1] = b * sg / 0.529177249;
  cm[2][1] = 0.0;

  cm[0][2] = c * cb / 0.529177249;
  cm[1][2] = (c * (ca - cb * cg) / sg) / 0.529177249;
  cm[2][2] = V / (a * b * sg) / 0.529177249;

  CoordMinMax[0] = 0.0;
  CoordMinMax[1] = 0.0;
  CoordMinMax[2] = 0.0;

  CoordMinMax[3] = (a + b * cg + c * cb) / 0.529177249;
  CoordMinMax[4] = (b * sg + c * (ca - cb * cg) / sg) / 0.529177249;
  CoordMinMax[5] = V / (a * b * sg) / 0.529177249;

  NbSteps[0] = (int) ceil(a / Resolution);
  NbSteps[1] = (int) ceil(b / Resolution);
  NbSteps[2] = (int) ceil(c / Resolution);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      cm[i][j] /= NbSteps[j];

  cif_input.close();
}

const double normgauss(const int& type, const double& exp)
{
  int t[3];
  if (type > 0)
    type2vector(type, t);
  else
    t[0] = t[1] = t[2] = 0;
  int temp = ft[t[0]] * ft[t[1]] * ft[t[2]];
  int temp2 = ft[2 * t[0]] * ft[2 * t[1]] * ft[2 * t[2]];
  return pow(2 * exp / PI, 0.75) * sqrt(pow(8 * exp, t[0] + t[1] + t[2]) * temp / temp2);
};
bool generate_sph2cart_mat(vector<vector<double>>&p, vector<vector<double>>& d, vector<vector<double>>& f, vector<vector<double>>& g) {
  //                                                   
  //From 3P: P0 P1 P2
  //To 3P : Z X Y (4 2 3, as in ORCA format)
  // 
  p.resize(3);
#pragma omp parallel for
  for (int i = 0; i < 3; i++) {
    p[i].resize(3, 0.0);
  }
  p[0][2] = 1.0;
  p[1][0] = 1.0;
  p[2][1] = 1.0;
  
  //                                                   
  //From 5D: D 0, D + 1, D - 1, D + 2, D - 2           
  //To 6D : 5  6  7  8  9 10                           
  //XX, YY, ZZ, XY, XZ, YZ      
  // 
  d.resize(6);
#pragma omp parallel for
  for (int i = 0; i < 6; i++) {
    d[i].resize(5, 0.0);
  }
  //XX = -0.5/SQRT(3) * D0 + 0.5 * D2
  d[0][0] = -0.5 / sqrt(3);
  d[0][3] = 0.5;
  //YY = -0.5/SQRT(3) * D0 - 0.5 * D2
  d[1][0] = -0.5 / sqrt(3);
  d[1][3] = -0.5;
  //ZZ = SQRT(1/3) * D0
  d[2][0] = sqrt(1.0 / 3.0);
  //XY = D-2
  d[3][4] = 1.0;
  //XZ = D1
  d[4][1] = 1.0;
  //YZ = D-1
  d[5][2] = 1.0;

  //From 7F: F 0, F + 1, F - 1, F + 2, F - 2, F + 3, F - 3
  //To 10F : 11   12   13   14   15   16   17   18   19  20
  //XXX, YYY, ZZZ, XXY, XXZ, YYZ, XYY, XZZ, YZZ, XYZ (AIMALL order!)
  //
  f.resize(10);
#pragma omp parallel for
  for (int i = 0; i < 10; i++) {
    f[i].resize(7, 0.0);
  }
  f[0][1] = -sqrt(0.025);
  f[0][5] = -sqrt(1.0 / 24.0);

  f[1][2] = -sqrt(0.025);
  f[1][6] = sqrt(1.0 / 24.0);

  f[2][0] = sqrt(1.0 / 15.0);

  f[3][2] = -sqrt(0.025);
  f[3][6] = -sqrt(0.375);

  f[4][0] = -sqrt(0.15);
  f[4][3] = 0.5;

  f[5][0] = -sqrt(0.15);
  f[5][3] = -0.5;

  f[6][1] = -sqrt(0.025);
  f[6][5] = sqrt(0.375);

  f[7][1] = sqrt(0.4);

  f[8][2] = sqrt(0.4);

  f[9][4] = 1.0;

  g.resize(15);
#pragma omp parallel for
  for (int i = 0; i < 15; i++) {
    g[i].resize(9, 0.0);
  }
  g[0][0] = 0.375 * sqrt(1.0 / 35.0);
  g[0][3] = -0.25 / sqrt(7);
  g[0][7] = -0.125;

  g[1][0] = g[0][0];
  g[1][3] = -g[0][3];
  g[1][7] = g[0][7];

  g[2][0] = sqrt(1.0 / 35.0);

  g[3][4] = -sqrt(1.0 / 28.0);
  g[3][8] = -0.5;

  g[4][1] = -0.5 / sqrt(98.0 / 63.0);
  g[4][5] = -1.0 / sqrt(8.0);

  g[5][4] = g[3][4];
  g[5][8] = -g[3][8];

  g[6][2] = g[4][1];
  g[6][6] = -g[4][5];

  g[7][1] = sqrt(2.0 / 7.0);

  g[8][2] = g[7][1];

  g[9][0] = 3.0 / sqrt(560);
  g[9][7] = 0.75;

  g[10][0] = -3.0 / sqrt(35);
  g[10][3] = 1.5 / sqrt(7);

  g[11][0] = g[10][0];
  g[11][3] = -g[10][3];

  g[12][2] = g[4][1];
  g[12][6] = -0.75 * sqrt(2);

  g[13][1] = g[4][1];
  g[13][5] = -g[12][6];

  g[14][4] = 3.0 / sqrt(7);
  return true;
}
bool generate_cart2sph_mat(vector<vector<double>>& d, vector<vector<double>>& f, vector<vector<double>>& g, vector<vector<double>>& h)
{
  //                                                   
  //From 5D: D 0, D + 1, D - 1, D + 2, D - 2           
  //To 6D : 1  2  3  4  5  6                           
  //XX, YY, ZZ, XY, XZ, YZ      
  // 
  d.resize(6);
#pragma omp parallel for
  for (int i = 0; i < 6; i++) {
    d[i].resize(5, 0.0);
  }
  //D0 = -0.5 * XX - 0.5 * YY + ZZ
  d[0][0] = -0.5;
  d[1][0] = -0.5;
  d[2][0] = 1.0;
  //D + 1 = XZ
  d[4][1] = 1.0;
  //D - 1 = YZ
  d[5][2] = 1.0;
  //D + 2 = SQRT(3) / 2 * (XX - YY)
  d[0][3] = sqrt(3.0) / 2.0;
  d[1][3] = -sqrt(3.0) / 2.0;
  //D - 2 = XY
  d[3][4] = 1.0;

  //From 7F: F 0, F + 1, F - 1, F + 2, F - 2, F + 3, F - 3
  //To 10F : 1   2   3   4   5   6   7   8   9  10
  //XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ(Gaussian sequence, not identical to Multiwfn)
  //
  f.resize(10);
#pragma omp parallel for
  for (int i = 0; i < 10; i++) {
    f[i].resize(7, 0.0);
  }
  //F 0 = -3 / (2 * sqrt5) * (XXZ + YYZ) + ZZZ
  f[2][0] = 1.0;
  f[5][0] = -1.5 / sqrt(5.0);
  f[8][0] = -1.5 / sqrt(5.0);
  //F + 1 = -sqrt(3 / 8) * XXX - sqrt(3 / 40) * XYY + sqrt(6 / 5) * XZZ
  f[0][1] = -sqrt(3.0 / 8.0);
  f[3][1] = -sqrt(3.0 / 40.0);
  f[6][1] = sqrt(6.0 / 5.0);
  //F - 1 = -sqrt(3 / 40) * XXY - sqrt(3 / 8) * YYY + sqrt(6 / 5) * YZZ
  f[1][2] = -sqrt(3.0 / 8.0);
  f[4][2] = -sqrt(3.0 / 40.0);
  f[7][2] = sqrt(6.0 / 5.0);
  //F + 2 = sqrt3 / 2 * (XXZ - YYZ)
  f[5][3] = sqrt(3.0) / 2.0;
  f[8][3] = -sqrt(3.0) / 2.0;
  //F - 2 = XYZ
  f[9][4] = 1.0;
  //F + 3 = sqrt(5 / 8) * XXX - 3 / sqrt8 * XYY
  f[0][5] = sqrt(5.0 / 8.0);
  f[3][5] = -3.0 / sqrt(8.0);
  //F - 3 = 3 / sqrt8 * XXY - sqrt(5 / 8) * YYY
  f[1][6] = -sqrt(5.0 / 8.0);
  f[4][6] = 3.0 / sqrt(8.0);

  //From 9G: G 0, G + 1, G - 1, G + 2, G - 2, G + 3, G - 3, G + 4, G - 4
  //To 15G : 1    2    3    4    5    6    7    8
  //ZZZZ, YZZZ, YYZZ, YYYZ, YYYY, XZZZ, XYZZ, XYYZ
  //9   10   11   12   13   14   15
  //XYYY, XXZZ, XXYZ, XXYY, XXXZ, XXXY, XXXX
  //
  g.resize(15);
#pragma omp parallel for
  for (int i = 0; i < 15; i++) {
    g[i].resize(9, 0.0);
  }
  //G 0 = ZZZZ + 3 / 8 * (XXXX + YYYY) - 3 * sqrt(3 / 35) * (XXZZ + YYZZ - 1 / 4 * XXYY)
  g[0][0] = 1.0;
  g[2][0] = -3.0 * sqrt(3.0 / 35.0);
  g[4][0] = 3.0 / 8.0;
  g[9][0] = -3.0 * sqrt(3.0 / 35.0);
  g[11][0] = 3.0 / 4.0 * sqrt(3.0 / 35.0);
  g[14][0] = 3.0 / 8.0;
  //G + 1 = 2 * sqrt(5 / 14) * XZZZ - 3 / 2 * sqrt(5 / 14) * XXXZ - 3 / 2 / sqrt14 * XYYZ
  g[5][1] = 2.0 * sqrt(5.0 / 14.0);
  g[7][1] = -1.5 / sqrt(14.0);
  g[12][1] = -1.5 * sqrt(5.0 / 14.0);
  //G - 1 = 2 * sqrt(5 / 14) * YZZZ - 3 / 2 * sqrt(5 / 14) * YYYZ - 3 / 2 / sqrt14 * XXYZ
  g[1][2] = 2.0 * sqrt(5.0 / 14.0);
  g[3][2] = -1.5 * sqrt(5.0 / 14.0);
  g[10][2] = -1.5 / sqrt(14.0);
  //G + 2 = 3 * sqrt(3 / 28) * (XXZZ - YYZZ) - sqrt5 / 4 * (XXXX - YYYY)
  g[2][3] = -3.0 * sqrt(3.0 / 28.0);
  g[4][3] = sqrt(5.0) / 4.0;
  g[9][3] = 3.0 * sqrt(3.0 / 28.0);
  g[14][3] = -sqrt(5.0) / 4.0;
  //G - 2 = 3 / sqrt7 * XYZZ - sqrt(5 / 28) * (XXXY + XYYY)
  g[6][4] = 3.0 / sqrt(7.0);
  g[8][4] = -sqrt(5.0 / 28.0);
  g[13][4] = -sqrt(5.0 / 28.0);
  //G + 3 = sqrt(5 / 8) * XXXZ - 3 / sqrt8 * XYYZ
  g[7][5] = -3.0 / sqrt(8.0);
  g[12][5] = sqrt(5.0 / 8.0);
  //G - 3 = -sqrt(5 / 8) * YYYZ + 3 / sqrt8 * XXYZ
  g[3][6] = -sqrt(5.0 / 8.0);
  g[10][6] = 3.0 / sqrt(8.0);
  //G + 4 = sqrt35 / 8 * (XXXX + YYYY) - 3 / 4 * sqrt3 * XXYY
  g[4][7] = sqrt(35.0) / 8.0;
  g[11][7] = -3.0 / 4.0 * sqrt(3.0);
  g[14][7] = sqrt(35.0) / 8.0;
  //G - 4 = sqrt5 / 2 * (XXXY - XYYY)
  g[8][8] = -sqrt(5.0) / 2.0;
  g[13][8] = sqrt(5.0) / 2.0;

  //From 11H: H 0, H + 1, H - 1, H + 2, H - 2, H + 3, H - 3, H + 4, H - 4, H + 5, H - 5
  //To 21H : 1     2     3     4     5     6     7     8     9    10
  //ZZZZZ YZZZZ YYZZZ YYYZZ YYYYZ YYYYY XZZZZ XYZZZ XYYZZ XYYYZ
  //11    12    13    14    15    16    17    18    19    20    21
  //XYYYY XXZZZ XXYZZ XXYYZ XXYYY XXXZZ XXXYZ XXXYY XXXXZ XXXXY XXXXX
  //
  h.resize(21);
#pragma omp parallel for
  for (int i = 0; i < 21; i++) {
    h[i].resize(11);
    std::fill(h[i].begin(), h[i].end(), 0.0);
  }
  //H 0 = ZZZZZ - 5 / sqrt21 * (XXZZZ + YYZZZ) + 5 / 8 * (XXXXZ + YYYYZ) + sqrt(15 / 7) / 4 * XXYYZ
  h[0][0] = 1.0;
  h[11][0] = -5.0 / sqrt(21.0);
  h[2][0] = -5.0 / sqrt(21.0);
  h[18][0] = 5.0 / 8.0;
  h[4][0] = 5.0 / 8.0;
  h[13][0] = sqrt(15.0 / 7.0) / 4.0;
  //H + 1 = sqrt(5 / 3) * XZZZZ - 3 * sqrt(5 / 28) * XXXZZ - 3 / sqrt28 * XYYZZ + sqrt15 / 8 * XXXXX + sqrt(5 / 3) / 8 * XYYYY + sqrt(5 / 7) / 4 * XXXYY
  h[6][1] = sqrt(5.0 / 3.0);
  h[15][1] = -3.0 * sqrt(5.0 / 28.0);
  h[8][1] = -3.0 / sqrt(28.0);
  h[20][1] = sqrt(15.0) / 8.0;
  h[10][1] = sqrt(5.0 / 3.0) / 8.0;
  h[17][1] = sqrt(5.0 / 7.0) / 4.0;
  //H - 1 = sqrt(5 / 3) * YZZZZ - 3 * sqrt(5 / 28) * YYYZZ - 3 / sqrt28 * XXYZZ + sqrt15 / 8 * YYYYY + sqrt(5 / 3) / 8 * XXXXY + sqrt(5 / 7) / 4 * XXYYY
  h[1][2] = sqrt(5.0 / 3.0);
  h[3][2] = -3.0 * sqrt(5.0 / 28.0);
  h[12][2] = -3.0 / sqrt(28.0);
  h[5][2] = sqrt(15.0) / 8.0;
  h[19][2] = sqrt(5.0 / 3.0) / 8.0;
  h[14][2] = sqrt(5.0 / 7.0) / 4.0;
  //H + 2 = sqrt5 / 2 * (XXZZZ - YYZZZ) - sqrt(35 / 3) / 4 * (XXXXZ - YYYYZ)
  h[11][3] = sqrt(5.0) / 2.0;
  h[2][3] = -sqrt(5.0) / 2.0;
  h[18][3] = -sqrt(35.0 / 3.0) / 4.0;
  h[4][3] = sqrt(35.0 / 3.0) / 4.0;
  //H - 2 = sqrt(5 / 3) * XYZZZ - sqrt(5 / 12) * (XXXYZ + XYYYZ)
  h[7][4] = sqrt(5.0 / 3.0);
  h[16][4] = -sqrt(5.0 / 12.0);
  h[9][4] = -sqrt(5.0 / 12.0);
  //H + 3 = sqrt(5 / 6) * XXXZZ - sqrt(3 / 2) * XYYZZ - sqrt(35 / 2) / 8 * (XXXXX - XYYYY) + sqrt(5 / 6) / 4 * XXXYY
  h[15][5] = sqrt(5.0 / 6.0);
  h[8][5] = -sqrt(1.5);
  h[20][5] = -sqrt(17.5) / 8.0;
  h[10][5] = sqrt(17.5) / 8.0;
  h[17][5] = sqrt(5.0 / 6.0) / 4.0;
  //H - 3 = -sqrt(5 / 6) * YYYZZ + sqrt(3 / 2) * XXYZZ - sqrt(35 / 2) / 8 * (XXXXY - YYYYY) - sqrt(5 / 6) / 4 * XXYYY
  h[3][6] = -sqrt(5.0 / 6.0);
  h[12][6] = sqrt(1.5);
  h[19][6] = -sqrt(17.5) / 8.0;
  h[5][6] = sqrt(17.5) / 8.0;
  h[14][6] = -sqrt(5.0 / 6.0) / 4.0;
  //H + 4 = sqrt35 / 8 * (XXXXZ + YYYYZ) - 3 / 4 * sqrt3 * XXYYZ
  h[18][7] = sqrt(35.0) / 8.0;
  h[4][7] = sqrt(35.0) / 8.0;
  h[13][7] = -0.75 * sqrt(3.0);
  //H - 4 = sqrt5 / 2 * (XXXYZ - XYYYZ)
  h[16][8] = sqrt(5.0) / 2.0;
  h[9][8] = -sqrt(5.0) / 2.0;
  //H + 5 = 3 / 8 * sqrt(7 / 2) * XXXXX + 5 / 8 * sqrt(7 / 2) * XYYYY - 5 / 4 * sqrt(3 / 2) * XXXYY
  h[20][9] = 3.0 / 8.0 * sqrt(3.5);
  h[10][9] = 5.0 / 8.0 * sqrt(3.5);
  h[17][9] = -1.25 * sqrt(1.5);
  //H - 5 = 3 / 8 * sqrt(7 / 2) * YYYYY + 5 / 8 * sqrt(7 / 2) * XXXXY - 5 / 4 * sqrt(3 / 2) * XXYYY
  h[5][10] = 3.0 / 8.0 * sqrt(3.5);
  h[19][10] = 5.0 / 8.0 * sqrt(3.5);
  h[14][10] = -1.25 * sqrt(1.5);
  return true;
}


const int type_vector[168]{
  0, 0, 0,
  1, 0, 0,
  0, 1, 0,
  0, 0, 1,
  2, 0, 0,
  0, 2, 0,
  0, 0, 2,
  1, 1, 0,
  1, 0, 1,
  0, 1, 1,
  3, 0, 0,
  0, 3, 0,
  0, 0, 3,
  2, 1, 0,
  2, 0, 1,
  0, 2, 1,
  1, 2, 0,
  1, 0, 2,
  0, 1, 2,
  1, 1, 1,
  0, 0, 4,
  0, 1, 3,
  0, 2, 2,
  0, 3, 1,
  0, 4, 0,
  1, 0, 3,
  1, 1, 2,
  1, 2, 1,
  1, 3, 0,
  2, 0, 2,
  2, 1, 1,
  2, 2, 0,
  3, 0, 1,
  3, 1, 0,
  4, 0, 0,
  0, 0, 5,
  0, 1, 4,
  0, 2, 3,
  0, 3, 2,
  0, 4, 1,
  0, 5, 0,
  1, 0, 4,
  1, 1, 3,
  1, 2, 2,
  1, 3, 1,
  1, 4, 0,
  2, 0, 3,
  2, 1, 2,
  2, 2, 1,
  2, 3, 0,
  3, 0, 2,
  3, 1, 1,
  3, 2, 0,
  4, 0, 1,
  4, 1, 0,
  5, 0, 0 };

void type2vector(
  const int& index,
  int* vector)
{
  if (index < 1 || index > 56) {
    vector[0] = -1;
    vector[1] = -1;
    vector[2] = -1;
    return;
  }
  const int temp = index - 1;
  vector[0] = type_vector[temp * 3];
  vector[1] = type_vector[temp * 3 + 1];
  vector[2] = type_vector[temp * 3 + 2];
}

bool read_fchk_integer_block(ifstream& in, string heading, vector<int>& result, bool rewind)
{
  if (result.size() != 0) result.clear();
  string line = go_get_string(in, heading, rewind);
  int limit = read_fchk_integer(line);
  int run = 0;
  int temp;
  getline(in, line);
  while (run < limit) {
    if (in.eof()) return false;
    temp = stoi(line.substr(12 * (run % 6), 12 * (run % 6 + 1)));
    result.push_back(temp);
    run++;
    if (run % 6 == 0)
      getline(in, line);
  }
  return true;
};
bool read_fchk_double_block(ifstream& in, string heading, vector<double>& result, bool rewind)
{
  if (result.size() != 0) result.clear();
  string line = go_get_string(in, heading, rewind);
  int limit = read_fchk_integer(line);
  int run = 0;
  double temp;
  getline(in, line);
  while (run < limit) {
    if (in.eof()) return false;
    temp = stod(line.substr(16 * (run % 5), 16 * (run % 5 + 1)));
    result.push_back(temp);
    run++;
    if (run % 5 == 0)
      getline(in, line);
  }
  return true;
};
int read_fchk_integer(string in)
{
  return stoi(in.substr(49, in.length() - 49));
};
double read_fchk_double(string in)
{
  return stod(in.substr(49, in.length() - 49));
};
int read_fchk_integer(std::ifstream& in, std::string search, bool rewind)
{
  string temp = go_get_string(in, search, rewind);
  return stoi(temp.substr(49, temp.length() - 49));
};
double read_fchk_double(std::ifstream& in, std::string search, bool rewind)
{
  string temp = go_get_string(in, search, rewind);
  return stod(temp.substr(49, temp.length() - 49));
};
void swap_sort(std::vector<int> order, std::vector< std::complex<double> >& v)
{
  int i = 0;
  while (i < v.size() - 1) {
    int new_index = 0;
    for (int j = i; j < v.size(); j++)
      if (order[j] < order[i])
        new_index++;
    if (new_index > 0) {
      std::complex<double> temp = v[i];
      v[i] = v[i + new_index];
      v[i + new_index] = temp;
      int temp2 = order[i];
      order[i] = order[i + new_index];
      order[i + new_index] = temp2;
    }
    else
      i++;
  }
}

void swap_sort_multi(std::vector<int> order, std::vector<std::vector<int>>& v)
{
  int i = 0;
  std::vector<int> temp;
  temp.resize(v.size());
  while (i < v.size() - 1) {
    int new_index = 0;
#pragma omp parallel for reduction(+:new_index)
    for (int j = i; j < v.size(); j++)
      if (order[j] < order[i])
        new_index++;
    if (new_index > 0) {
#pragma omp parallel for
      for (int run = 0; run < v.size(); run++) {
        temp[run] = v[run][i];
        v[run][i] = v[run][i + new_index];
        v[run][i + new_index] = temp[run];
      }
      int temp2 = order[i];
      order[i] = order[i + new_index];
      order[i + new_index] = temp2;
    }
    else
      i++;
  }
}

double get_lambda_1(double* a)
{
  vector<double> bw, zw;
  //int run = 0;
  double eig1, eig2, eig3;
  const double p1 = a[1] * a[1] + a[2] * a[2] + a[5] * a[5];
  if (p1 == 0) {
    eig1 = a[0];
    eig2 = a[4];
    eig3 = a[8];
    if ((eig1 < eig2 && eig2 < eig3) || (eig3 < eig2 && eig2 < eig1))
      return eig2;

    else if ((eig2 < eig1 && eig1 < eig3) || (eig3 < eig1 && eig1 < eig2))
      return eig1;

    else
      return eig3;
  }
  else {
    const double q = (a[0] + a[4] + a[8]) / 3;
    const double p2 = pow(a[0] - q, 2) + pow(a[4] - q, 2) + pow(a[8] - q, 2) + 2 * p1;
    const double p = sqrt(p2 / 6);
    const double B[9]{
    a[0] - q,
    a[1],
    a[2],
    a[3],
    a[4] - q,
    a[5],
    a[6],
    a[7],
    a[8] - q };
    const double r = (B[0] * B[4] * B[8]
      + B[1] * B[5] * B[6]
      + B[3] * B[4] * B[7]
      - B[0] * B[5] * B[7]
      - B[1] * B[3] * B[8]
      - B[2] * B[4] * B[6]) / 2;
    double phi;
    if (r <= -1)
      phi = PI / 3;
    else if (r >= 1)
      phi = 0;
    else
      phi = acos(r) / 3;

    eig1 = q + 2 * p * cos(phi);
    eig3 = q + 2 * p * cos(phi + 2 * PI / 3);
    eig2 = 3 * q - eig1 - eig3;
    if ((eig1 < eig2 && eig2 < eig3) || (eig3 < eig2 && eig2 < eig1))
      return eig2;

    else if ((eig2 < eig1 && eig1 < eig3) || (eig3 < eig1 && eig1 < eig2))
      return eig1;

    else
      return eig3;
  }
};

double get_decimal_precision_from_CIF_number(string& given_string) {
  int len = (int) given_string.length();
  int open_bracket = -1;
  int close_bracket = -1;
  int decimal_point = -1;
  //const char* gs = given_string.c_str();
  for (int i = 0; i < len; i++) {
    if (given_string[i] == '(' && open_bracket == -1) {
      open_bracket = i;
    }
    else if (given_string[i] == ')' && close_bracket == -1) {
      close_bracket = i;
    }
    else if (given_string[i] == '.' && decimal_point == -1) {
      decimal_point = i;
    }
  }
  double result = 0;
  int precision = 1;
  int size_of_precision = 1;
  if (open_bracket != -1 && close_bracket != -1) {
    size_of_precision = close_bracket - open_bracket - 1;
    string temp = given_string.substr(open_bracket + 1, size_of_precision);
    precision = stoi(temp);
  }
  int digits = 0;
  if (open_bracket != -1 && close_bracket != -1) {
    if (decimal_point != -1) {
      digits = open_bracket - decimal_point - 1;
    }
    else {
      digits = close_bracket - open_bracket - 1;
    }
    if (digits == 0) {
      return 0.001;
    }
    result = abs(precision * pow(10, -digits));
    return result;
  }
  else return 0.005;
};