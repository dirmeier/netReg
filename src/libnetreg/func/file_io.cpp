
/**
 * Author: Simon Dirmeier
 * Date: 3/17/16
 * Email: uccd_.dirmeier@bsse.ethz.ch
 */

#include "file_io.hpp"

#include <iostream>
#include <fstream>

#include "utility_functions.hpp"
#include "../exceptions/io_exception.hpp"

namespace netreg
{
    matrix<double> read_matrix(const std::string &file)
    {
        std::vector< std::vector<double> > content = read_content(file);
        return to_matrix(content);
    }

    std::vector< std::vector<double> >
        read_content(const std::string &file)
    {
        std::vector < std::vector<double> > content;
        std::ifstream ifs(file);
        if (!ifs.is_open())
            throw io_exception();
        std::string line;
        while (std::getline(ifs, line))
        {
            std::vector<double> vals = to_double(line, '\t');
            content.push_back(vals);
        }
        ifs.close();
        return content;
    }

    void write_results(linear_model_data & data, std::string& output)
    {
        std::ofstream of(output);
        cvector<double>& intercepts = data.intercept();
        matrix<double> & coeffs = data.coefficients();
        of << "#INTERCEPTS\n";
        for (std::vector<double>::size_type i = 0;
            i < intercepts.size();
            i++)
        {
            of << intercepts[i];
            if (i < intercepts.size() - 1)
                of << "\t";
        }
        of << "\n";
        of << "#COEFFICIENTS\n";
        for(int i = 0; i < coeffs.n_rows; i++)
        {
            for(int j = 0; j < coeffs.n_cols; j++)
            {
                of << coeffs(i, j);
                if (j < coeffs.n_cols - 1)
                    of << "\t";
            }
            if (i < coeffs.n_rows - 1)
                of << "\n";
        }
        of.close();
    }

    void write_results(pareto_optimal_point<std::string, double>& data, std::string& output)
    {
        std::ofstream of(output);
        of << "#PARETO OPTIMAL VALUES\n";
        for (int i = 0; i < data.size(); ++i)
        {
            std::pair<std::string, double> entry = data[i];
            of << entry.first << "\t" << entry.second << "\n";
        }
        of.close();
    }
};