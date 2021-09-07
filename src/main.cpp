#include "all.h"

extern int cov_in_motions_main(int argc, char** argv);

void PrintHelp(std::string aArg)
{
    std::cerr << "Usage: " << aArg << " "
              << "Options:\n"
              << "\t-h,--help\tShow this help message\n"
              << "\trelm.txt \tmotion file\n"
              << "\tcoords.txt\t2d coordinate file\n";
            //  << "\tShow?\t 1 or 0 for verbose print"
}

int main(int argc, char** argv)
{
  std::cout << "cov_in_motions" << "\n";

  if (argc < 2)
  {
    PrintHelp(argv[0]);
    std::cout << argc << "\n";
    return 1;
  }

  std::string aArg1 = std::string(argv[1]);

  if ((aArg1 == "-h") || (aArg1 == "--help") || (aArg1 == "-help"))
  {
      std::cout << aArg1 << "\n";
      PrintHelp(argv[0]);
      return 1;
  }


  return cov_in_motions_main(argc,argv);


}
