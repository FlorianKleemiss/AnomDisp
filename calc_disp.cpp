#include "convenience.h"
#include "cell.h"
#include "spherical_density.h"
#include "matrix_coefficient.h"
#include "intensity_calculation.h"
#include "brennan.h"
#include <gsl/gsl_integration.h>

void plot(const vec &x, 
          const std::vector<vec> &y,
          const std::vector<std::string> &labels,
          const std::string &data_filename = "data.txt",
          const std::string &plot_filename = "",
          const double &xmin = 0.15,
          const double &xmax = 3.0)
{
  std::string command_filename = "commands.txt";
  std::ofstream command;
  std::ofstream data;
  std::string lpfn;
  if(plot_filename == "") lpfn = get_basename_without_ending(data_filename) + ".png";
  else lpfn = plot_filename;

//  Create the data file.
  for(int i=0; i<y.size(); i++){
    std::string ldfn = data_filename + "_" + std::to_string(i);
    data.open ( ldfn.c_str ( ) );
    for (int j = 0; j < x.size(); j++){
      data << "  " << x[j] << "  " << y[i][j] << "\n";
    }
    data.close ( );
  }

  std::cout << "\n";
  std::cout << "  plot: data stored in '" << data_filename << " " << 0 << "-" << y.size() << "'\n";
       
//  Create the command file.

  command.open ( command_filename.c_str ( ) );
  command << "# " << command_filename << "\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < " << command_filename << "\n";
  command << "#\n";
  command << "set term png size 1920,1080 font 'Helvetica,25'\n";
  command << "set output '" << lpfn << "'\n";
  command << "set xlabel 'wavelength'\n";
  command << "set ylabel 'fdp'\n";
  command << "set logscale x\n";
  command << "set xrange [" << xmin << ":" << xmax << "]\n";
  command << "set title 'fdp plot for " << data_filename << "'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  for(int i=0; i<y.size(); i++){
    std::string ldfn = data_filename + "_" + std::to_string(i);
    if (i==0)
      command << "plot '" << ldfn << "' using 1:2 with lines title '" << labels[i] << "'";
    else 
      command << "     '" << ldfn << "' using 1:2 with lines title '" << labels[i] << "'";
    if (i<y.size()-1)
      command << ", \\\n";
    else
      command << "\n";
  }
  command << "quit\n";

  command.close ( );

  std::cout << "  plot: plot commands stored in '"
       << command_filename << "'\n";

  return;
}

int main (int argc, char** argv){
  int steps = 1000;
  vec x_axis(steps);
  std::vector<vec> y_axis;
  y_axis.push_back(vec(steps));
  y_axis.push_back(vec(steps));
  std::vector<std::string> labels({"Kleemiss,Peyerimhoff", "Brennan,Cowan"});
  double min, max, t0 = 0;
  min = 0.15;
  max = 3.0;
  int Z = 52;
  int sel;
  std::cout << "Which atomic number do you want to calculate for? ";
  std::cin >> sel;
  if (sel > 0) Z = sel;
#pragma omp parallel for
  for(int i=0; i<steps; i++){
    double wl = min+(max-min)/steps*i;
    x_axis[i] = wl;
    y_axis[0][i] = kleemiss_peyerimhoff::at_angstrom(wl,t0,Z,5,4);
    y_axis[1][i] = brennan::at_angstrom(wl,Z)[1];
  }
  plot(x_axis,y_axis,labels,atnr2letter(Z)+".fdp","",min,max);

  std::system("gnuplot\\bin\\gnuplot.exe commands.txt");
  std::cout << "Plotting finished!" << std::endl;
  std::remove("commands.txt"); // delete file
  std::cout << "Done!" << std::endl;
  return 0;

}