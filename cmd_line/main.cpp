#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <tuple>
#include <vector>
#include <exception>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "../src/family.hpp"
#include "../src/edgenet_gaussian.hpp"
#include "../src/graph_penalized_linear_model_cv_data.hpp"
#include "../src/edgenet_gaussian_model_selection.hpp"

int main(int argc, char *argv[])
{
    std::cout.imbue(std::locale("en_US.UTF-8"));

    std::string design_filename;
    std::string response_filename;
    std::string parameter_filename;
    std::string out_filename;

    double threshold;
    uint32_t maxit;
    uint8_t nfolds;

    boost::program_options::options_description generic("Generic options");
    generic.add_options()("help,h", "Print this help.");

    boost::program_options::options_description files("Files");
    files.add_options()(
",d",
                        boost::program_options::value<decltype(design_filename)>(
            &design_filename)
            ->required(),
        "Filename of the design matrix X.")
        (",r",
        boost::program_options::value<decltype(response_filename)>(
            &response_filename)
            ->required(),
        "Filename of the reponse matrix Y.")(
        ",p",
        boost::program_options::value<decltype(parameter_filename)>(
            &parameter_filename),
        "Filename of an optional file of parameters.")(
        ",o",
        boost::program_options::value<decltype(out_filename)>(&out_filename)
            ->required(),
        "Filename of the output file.");

    boost::program_options::options_description config("Configuration");
    config.add_options()(
        "maxit,m",
        boost::program_options::value<decltype(maxit)>(&maxit)->default_value(
            100000, "100000"),
        "Maximum number of iterations of coordinate descent."
        " You should choose a sufficiently large number.")(
        "nfolds,n",
        boost::program_options::value<decltype(nfolds)>(&nfolds)
          ->default_value(10, "10"),
        "The number of cross-validation folds. "
          "This can be maximal the number of rows of X/Y and minimal 3.")(
          "threshold,t",
        boost::program_options::value<decltype(threshold)>(&threshold)
            ->default_value(0.0000001, "0.0000001"),
        "Convergence threshold for coordinate descent. "
          "Anything below 0.0001 should suffice.");

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(files);
        po::positional_options_description p;
  p.add("input-file", -1);
    po::variables_map global_options;

    try
    {
        po::store(po::command_line_parser(argc, argv)
                  .options(cmdline_options).positional(p)
                  .run(),
            global_options);
    }
    catch (const po::error const& e)
    {
        std::cerr << "asdasdsa" << '\n';
        std::cerr << e.what() << '\n';
    }

    boost::program_options::notify(global_options);
    std::cout << design_filename << "\n"
              << response_filename << "\n"
              << parameter_filename << "\n"
              << out_filename << "\n";

    //    const int n = 100;
    //    const int p = 100;
    //    const int q = 10;
    //
    //    double *x = new double[n * p];
    //    double *y = new double[n * q];
    //    double *gx = new double[p * p];
    //    double *gy = new double[q * q];
    //
    //    for (int j = 0; j < n * p; ++j) x[j] = 1;
    //    for (int j = 0; j < n * q; ++j) y[j] = 2;
    //    for (int k = 0; k < p * p; ++k) gx[k] = 1;
    //    for (int k = 0; k < q * q; ++k) gy[k] = 1;
    //
    //    /*  netreg::graph_penalized_linear_model_data dat
    //          (x, y, gx, gy, n, p, q,
    //           10, 1.0, 1.0, 1.0, 100000, 0.0001,
    //           netreg::family::GAUSSIAN);
    //
    //      netreg::edgenet_gaussian ed;
    //      ed.run(dat);
    //      std::cout << "done1" << std::endl;*/
    //
    //    netreg::graph_penalized_linear_model_cv_data data(
    //        x, y, gx, gy, n, p, q,
    //         -1, 1.0, 0.0, 0.0, 100000, 0.0000000001, 5,
    //         netreg::family::GAUSSIAN);
    //
    //    std::map<std::string, double> m = e.regularization_path(data, 1000,
    //                                                            0.001);
    //    std::cout << "Approx: " << m["lambda"] << " " << m["psigx"] << " "
    //              << m["psigy"] << std::endl;
    //
    //    // std::map<std::string, double> m = e.regularization_path(data,
    //    false, 1000, 0.001);
    //    // std::cout << "True: " <<  m["lambda"] << " " << m["psigx"] << " "
    //    << m["psigy"] << std::endl;
    //
    //    delete[] x;
    //    delete[] y;
    //    delete[] gx;
    //    delete[] gy;

    return 0;
}
