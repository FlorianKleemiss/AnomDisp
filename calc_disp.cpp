#include "convenience.h"
#include "cell.h"
#include "spherical_density.h"
#include "matrix_coefficient.h"
#include "intensity_calculation.h"
#include "brennan.h"
#include <gsl/gsl_integration.h>

void plot(const vec &x, 
          const vec &y, 
          const std::string &data_filename = "data.txt",
          const std::string &plot_filename = "",
          const double &xmin = 0.15,
          const double &xmax = 3.0)
{
  std::string command_filename = "commands.txt";
  std::ofstream command;
  std::ofstream data;
  int j;
  std::string lpfn;
  if(plot_filename == "") lpfn = get_basename_without_ending(data_filename) + ".png";
  else lpfn = plot_filename;

  std::cout << "\n";
  std::cout << "plot:\n";
  std::cout << "  Write command and data files that can be used\n";
  std::cout << "  by gnuplot for a plot.\n";

//  Create the data file.

  data.open ( data_filename.c_str ( ) );

  for ( j = 0; j < x.size(); j++ )
  {
    data << "  " << x[j]
         << "  " << y[j] << "\n";
  }

  data.close ( );

  std::cout << "\n";
  std::cout << "  plot: data stored in '" << data_filename << "'\n";
       
//  Create the command file.

  command.open ( command_filename.c_str ( ) );

  command << "# " << command_filename << "\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < " << command_filename << "\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output '" << lpfn << "'\n";
  command << "set xlabel 'wavelength'\n";
  command << "set ylabel 'fdp'\n";
  command << "set logscale x\n";
  command << "set xrange [" << xmin << ":" << xmax << "]\n";
  command << "set title 'fdp plot for " << data_filename << "'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:2 with lines\n";
  command << "quit\n";

  command.close ( );

  std::cout << "  plot: plot commands stored in '"
       << command_filename << "'\n";

  return;
}

int main (int argc, char** argv){
  int steps = 1000;
  vec x_axis(steps);
  vec y_axis(steps);
  double min, max, t0 = 0;
  min = 0.15;
  max = 5.0;
  int Z = 52;
  int sel;
  std::cout << "Which atomic number do you want to calculate for? ";
  std::cin >> sel;
  if (sel > 0) Z = sel;
#pragma omp parallel for
  for(int i=0; i<steps; i++){
    double wl = min+(max-min)/steps*i;
    x_axis[i] = wl;
    double nu = speed_of_light / wl * 1E10;
    y_axis[i] = calc_sum(nu,t0,7,1,Z);
  }
  plot(x_axis,y_axis,atnr2letter(Z)+".fdp","",min,max);

  std::system("gnuplot\\bin\\gnuplot.exe commands.txt");
  return 0;

}