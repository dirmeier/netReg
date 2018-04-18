/**
 * netReg: graph-regularized linear regression models.
 * <p>
 * Copyright (C) 2015 - 2016 Simon Dirmeier
 * <p>
 * This file is part of netReg.
 * <p>
 * netReg is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * netReg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with netReg. If not, see <http://www.gnu.org/licenses/>.
 *
 * @author: Simon Dirmeier
 * @email: simon.dirmeier@gmx.de
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <exception>
#include <algorithm>
#include <map>
#include <iterator>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "armadillo"

#include "../src/family.hpp"
#include "../src/params.hpp"
#include "../src/stat_functions.hpp"

#include "../src/data_factory.hpp"
#include "../src/graph_model_data.hpp"
#include "../src/graph_model_cv_data.hpp"


#include "../src/edgenet.hpp"
#include "../src/edgenet_model_selection.hpp"
#include "../src/cross_validation.hpp"

static const char* netReg =
  "\nnetReg - a network-regularized generalized regression model";

struct data_set
{
    std::vector<double> data;
    unsigned int nrow;
    unsigned int ncol;
};

static data_set read_tsv(const std::string& file);

static void fit(struct data_set& X,
                struct data_set& Y,
                const std::string& gx_filename,
                const std::string& gy_filename,
                const std::string& outfile,
                double lambda,
                double psi,
                double phi,
                double threshold,
                uint32_t maxit,
                uint32_t nfolds);

static std::map<std::string, double> modelselection(
  struct data_set& X,
  struct data_set& Y,
  const std::string& gx_filename,
  const std::string& gy_filename,
  const std::string& outfile,
  double lambda,
  double psi,
  double phi,
  double threshold,
  uint32_t maxit,
  uint32_t nfolds,
  uint32_t bobit,
  double epsilon);

static netreg::graph_model_data get_fit_data(
  struct data_set& X,
  struct data_set& Y,
  const std::string& gx_filename,
  const std::string& gy_filename);

static netreg::graph_model_cv_data get_ms_data(
  struct data_set& X,
  struct data_set& Y,
  const std::string& gx_filename,
  const std::string& gy_filename,
  uint32_t nfolds);

static void check_matrix(const data_set& aff,
                         const data_set& m,
                         const std::string& m_name,
                         unsigned int cl);

static void check_params(double lambda,
                         double psi,
                         double phi,
                         uint32_t maxit,
                         double threshold,
                         uint32_t nfolds,
                         double epsilon,
                         const struct data_set& X,
                         const struct data_set& Y);

int main(int argc, char* argv[])
{
    std::string design_filename = "";
    std::string response_filename = "";
    std::string gx_filename = "";
    std::string gy_filename = "";
    std::string out_filename = "";

    double threshold = .0000001;
    uint32_t maxit = 100000;

    double lambda = .1;
    double psi = .0;
    double phi = .0;

    uint32_t nfolds = 10;
    uint32_t bobit = 1000;
    double epsilon = 0.0001;

    // clang-format off
    boost::program_options::options_description desc("Arguments");
    desc.add_options()
          ("help,h", "Print this help.")
          ("design,d",
           boost::program_options::value<decltype(design_filename)>(
             &design_filename)->required(),
           "Filename of the design matrix X.")
          ("gx,u",
           boost::program_options::value<decltype(gx_filename)>(
             &gx_filename),
           "Filename of the affinity matrix GX for X.")
          ("reponse,r",
           boost::program_options::value<decltype(response_filename)>(
             &response_filename)->required(),
           "Filename of the reponse matrix Y.")
          ("gy,v",
           boost::program_options::value<decltype(gy_filename)>(
             &gy_filename),
           "Filename of the affinity matrix GY for Y.")
          ("lambda,l",
           boost::program_options::value<decltype(lambda)>(
             &lambda)->default_value(1, "1"),
           "LASSO penalization parameter.")
          ("psi,s",
           boost::program_options::value<decltype(psi)>(
             &psi)->default_value(0, "0"),
           "Penalization parameter for affinity matrix GX.")
          ("xi,x",
           boost::program_options::value<decltype(phi)>(
             &phi)->default_value(0, "0"),
           "Penalization parameter for affinity matrix GY.")
          ("maxit,m",
           boost::program_options::value<decltype(maxit)>(
             &maxit)->default_value(100000, "100000"),
           "Maximum number of iterations of coordinate descent."
             " You should choose a sufficiently large number.")
          ("threshold,t",
           boost::program_options::value<decltype(threshold)>(&threshold)
             ->default_value(0.0000001, "0.0000001"),
           "Convergence threshold for coordinate descent. "
             "Anything below 0.0001 should suffice.")
          ("outfile,o",
           boost::program_options::value<decltype(out_filename)>(
             &out_filename)->required(),
           "Filename of the output file.");

    boost::program_options::options_description mods("Model selection");
    mods.add_options()
          ("modelselection",
           "Use modelselection, i.e. estimation of optimal shrinkage "
             "parameters using crossvalition, "
             "before doing the estimation of coefficients.")
          ("nfolds,n",
           boost::program_options::value<decltype(nfolds)>(
             &nfolds)->default_value(
             10, "10"),
           "The number of cross-validation folds. "
             "This can be maximal the number of rows of X/Y and minimal 3.")
          ("epsilon,e",
           boost::program_options::value<decltype(epsilon)>(&epsilon)
             ->default_value(0.001, "0.001"),
           "Convergence threshold for the BOBYQA algorithm, "
             "i.e. the stop criterion for the model selection.")
          ("bobit,b",
           boost::program_options::value<decltype(bobit)>(
             &bobit)->default_value(1000, "1000"),
           "Maximal number of iterations for the BOBYQA algorithm.");
    // clang-format on

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(desc).add(mods);
    boost::program_options::positional_options_description p;
    p.add("input-file", -1);
    boost::program_options::variables_map vm;

    try
    {
        boost::program_options::store(
          boost::program_options::command_line_parser(argc, argv)
            .options(cmdline_options)
            .positional(p)
            .run(),
          vm);
        if (vm.count("help") || argc <= 1)
        {
            std::cout << netReg << "\n";;
            std::cout << cmdline_options << std::endl;
            exit(EXIT_SUCCESS);
        }
        boost::program_options::notify(vm);
    }
    catch (std::exception& e)
    {
        std::cout << netReg << "\n";
        std::cout << cmdline_options << std::endl;
        std::cout << e.what() << std::endl;
        exit(EXIT_SUCCESS);
    }

    struct data_set X = read_tsv(design_filename);
    struct data_set Y = read_tsv(response_filename);
    check_params(lambda, psi, phi, maxit, threshold, nfolds, epsilon, X, Y);

    if (vm.count("modelselection"))
    {
        if (X.nrow > 1000 || X.ncol > 500)
            nfolds = 5;

        std::cout << "Doing model selection.\n";
        std::map<std::string, double> opt = modelselection(
          X, Y, gx_filename, gy_filename,
          out_filename,
          lambda, psi, phi,
          threshold, maxit, nfolds,
          bobit, epsilon);

        lambda = opt["lambda"];
        psi = opt["psigx"];
        phi = opt["psigy"];
    }

    std::cout << "Fitting model" << std::endl;
    fit(X, Y, gx_filename, gy_filename,
        out_filename,
        lambda, psi, phi,
        threshold, maxit, nfolds);

    return EXIT_SUCCESS;
}

data_set read_tsv(const std::string& file)
{
    std::ifstream in(file.c_str());
    if (!in.is_open())
    {
        std::cerr << "IO-error with file: " << file << std::endl;
        exit(EXIT_FAILURE);
    }

    struct data_set data;

    std::string line;
    std::vector<std::string> strs;
    while (getline(in, line))
    {
        if (boost::starts_with(line, "#")) continue;
        boost::split(strs, line, boost::is_any_of("\t"));
        data.ncol = strs.size();
        for (std::vector<std::string>::size_type i = 0; i < strs.size(); ++i)
        {
            data.data.push_back(atof(strs[i].c_str()));
        }
    }
    in.close();
    data.nrow = data.data.size() / data.ncol;

    return data;
}

void fit(struct data_set& X,
         struct data_set& Y,
         const std::string& gx_filename,
         const std::string& gy_filename,
         const std::string& outfile,
         double lambda,
         double psi,
         double phi,
         double threshold,
         uint32_t maxit,
         uint32_t nfolds)
{
    netreg::graph_model_data dat = get_fit_data(
      X, Y, gx_filename, gy_filename);

    netreg::params pars = netreg::params()
      .lambda(lambda)
      .psigx(psi)
      .psigy(phi)
      .thresh(threshold)
      .niter(maxit);

    netreg::edgenet edge(dat, pars);

    arma::Mat<double> coef = edge.run();
    arma::Col<double> intr =
      netreg::intercept(dat.design(), dat.response(), coef);

    std::string coeffile = outfile.substr(0, outfile.find_last_of('.')) +
                           "_coefficients" +
                           outfile.substr(outfile.find_last_of('.'));
    std::string intfile = outfile.substr(0, outfile.find_last_of('.')) +
                          "_intercepts" +
                          outfile.substr(outfile.find_last_of('.'));
    std::cout << "Writing coefficients to:" << coeffile << std::endl;
    std::cout << "Writing intercept to:" << intfile << std::endl;

    coef.save(coeffile, arma::csv_ascii);
    intr.save(intfile, arma::csv_ascii);
}

std::map<std::string, double> modelselection(struct data_set& X,
                                             struct data_set& Y,
                                             const std::string& gx_filename,
                                             const std::string& gy_filename,
                                             const std::string& outfile,
                                             double lambda,
                                             double psi,
                                             double phi,
                                             double threshold,
                                             uint32_t maxit,
                                             uint32_t nfolds,
                                             uint32_t bobit,
                                             double epsilon)
{
    netreg::graph_model_cv_data dat = get_ms_data(
      X, Y, gx_filename, gy_filename, nfolds);

    netreg::params pars = netreg::params()
      .lambda(lambda)
      .psigx(psi)
      .psigy(phi)
      .do_lambda(true)
      .do_psigx(true)
      .do_psigy(true)
      .thresh(threshold)
      .niter(maxit)
      .optim_niter(bobit)
      .optim_epsilon(epsilon);

    std::map<std::string, double> m =\
      netreg::model_selection<netreg::cross_validator>(dat, pars);

    std::string paramfile = outfile.substr(0, outfile.find_last_of('.')) +
                            "_optimal_shrinkage_params" +
                            outfile.substr(outfile.find_last_of('.'));
    std::cout << "Writing shrinkage parameters to:" << paramfile << std::endl;
    std::ofstream stream(paramfile.c_str());
    for (auto& kv : m)
    {
        stream << kv.first << "\t" << kv.second << '\n';
    }
    stream.close();

    return m;
}

netreg::graph_model_data get_fit_data(
  struct data_set& X,
  struct data_set& Y,
  const std::string& gx_filename,
  const std::string& gy_filename)
{
    std::string G = "gaussian";

    if (gx_filename != "" && gy_filename != "")
    {
        struct data_set GX = read_tsv(gx_filename);
        struct data_set GY = read_tsv(gy_filename);
        check_matrix(GX, X, "GX", X.ncol);
        check_matrix(GY, Y, "GY", Y.ncol);

        return netreg::data_factory::build_data(
          X.data.data(), Y.data.data(),
          GX.data.data(), GY.data.data(),
          X.nrow, X.ncol, Y.ncol,
          G);
    }
    else if (gx_filename != "")
    {
        struct data_set GX = read_tsv(gx_filename);
        check_matrix(GX, X, "GX", X.ncol);

        return netreg::data_factory::build_data(
          X.data.data(), Y.data.data(),
          GX.data.data(), std::vector<double>(1).data(),
          X.nrow, X.ncol, Y.ncol,
          G);
    }
    else if (gy_filename != "")
    {
        struct data_set GY = read_tsv(gy_filename);
        check_matrix(GY, Y, "GY", Y.ncol);

        return netreg::data_factory::build_data(
          X.data.data(), Y.data.data(),
          std::vector<double>(1).data(), GY.data.data(),
          X.nrow, X.ncol, Y.ncol,
          G);
    }
    else
    {
        return netreg::data_factory::build_data(
          X.data.data(), Y.data.data(),
          std::vector<double>(1).data(), std::vector<double>(1).data(),
          X.nrow, X.ncol, Y.ncol,
          G);
    }
}

netreg::graph_model_cv_data get_ms_data(
  struct data_set& X,
  struct data_set& Y,
  const std::string& gx_filename,
  const std::string& gy_filename,
  uint32_t nfolds)
{
    int nf = static_cast<int>(nfolds);
    std::string G = "gaussian";

    if (gx_filename != "" && gy_filename != "")
    {
        struct data_set GX = read_tsv(gx_filename);
        struct data_set GY = read_tsv(gy_filename);
        check_matrix(GX, X, "GX", X.ncol);
        check_matrix(GY, Y, "GY", Y.ncol);

        return netreg::data_factory::build_cv_data(
          X.data.data(), Y.data.data(),
          GX.data.data(), GY.data.data(),
          X.nrow, X.ncol, Y.ncol,
          G, nf);
    }
    else if (gx_filename != "")
    {
        struct data_set GX = read_tsv(gx_filename);
        check_matrix(GX, X, "GX", X.ncol);

        return netreg::data_factory::build_cv_data(
          X.data.data(), Y.data.data(),
          GX.data.data(), std::vector<double>(1).data(),
          X.nrow, X.ncol, Y.ncol,
          G, nf);
    }
    else if (gy_filename != "")
    {
        struct data_set GY = read_tsv(gy_filename);
        check_matrix(GY, Y, "GY", Y.ncol);

        return netreg::data_factory::build_cv_data(
          X.data.data(), Y.data.data(),
          std::vector<double>(1).data(), GY.data.data(),
          X.nrow, X.ncol, Y.ncol,
          G, nf);
    }
    else
    {
        return netreg::data_factory::build_cv_data(
          X.data.data(), Y.data.data(),
          std::vector<double>(1).data(),
          std::vector<double>(1).data(),
          X.nrow, X.ncol, Y.ncol,
          G, nf);
    }
}

void check_params(double lambda,
                  double psi,
                  double phi,
                  uint32_t maxit,
                  double threshold,
                  uint32_t nfolds,
                  double epsilon,
                  const struct data_set& X,
                  const struct data_set& Y)
{
    if (X.nrow != Y.nrow)
    {
        std::cerr << "X and Y dont have the same number of observations.\n";
        exit(EXIT_FAILURE);
    }

    if (psi < 0 || phi < 0)
    {
        std::cerr << "Penalization parameters (lambda/psi/phi) need to be >0."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (maxit < 1000)
    {
        std::cerr << "Maxit should be at least 1000." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (threshold < 0)
    {
        std::cerr << "Threshold ahs to be >0." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (nfolds > X.nrow)
    {
        std::cerr << "Nfolds needs to be >0 and <nrow(X)." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (epsilon < 0)
    {
        std::cerr << "Epsilon needs to be >0." << std::endl;
        exit(EXIT_FAILURE);
    }
}

void check_matrix(const data_set& aff,
                  const data_set& m,
                  const std::string& m_name,
                  unsigned int cl)
{
    if (aff.nrow != aff.ncol || aff.nrow != m.ncol)
    {
        std::cerr << m_name << " needs to be a (" << cl << " x " << cl
                  << ")-dimensional matrix.\n";
        exit(EXIT_FAILURE);
    }
}
