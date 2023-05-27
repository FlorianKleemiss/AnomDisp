#include "convenience.h"
#include "cell.h"
#include "spherical_density.h"
#include "matrix_coefficient.h"
#include <gsl/gsl_integration.h>

double f (double x, void * params) {
  double alpha = *(double *) params;
  return log(alpha*x) / sqrt(x);
}
typedef std::vector<double> vec;

void plot(const vec &x, const vec &y, const std::string &data_filename = "data.txt")
{
  std::string command_filename = "commands.txt";
  std::ofstream command;
  std::ofstream data;
  int j;
  std::string plot_filename = get_basename_without_ending(data_filename) + ".png";

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
  command << "set output '" << plot_filename << "'\n";
  command << "set xlabel 'X'\n";
  command << "set ylabel 'Y'\n";
  command << "set title 'Plot using gnuplot'\n";
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
  double result, error;
  double expected = -4.0;
  double alpha = 1.0;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  F.function = &f;
  F.params = &alpha;

  gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
                        w, &result, &error); 

  printf ("result          = % .18f\n", result);
  printf ("exact result    = % .18f\n", expected);
  printf ("estimated error = % .18f\n", error);
  printf ("actual error    = % .18f\n", result - expected);
  printf ("intervals =  %d\n", int(w->size));

  gsl_integration_workspace_free (w);

  vec x {0,1};
  vec y {0,1};
  plot(x,y);

  std::system("gnuplot\\bin\\gnuplot.exe commands.txt");
  std::cout << "Works10!" << std::endl;
  return 0;

}