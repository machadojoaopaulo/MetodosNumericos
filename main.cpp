#include <iostream>
#include <math.h>
#include <string.h>
#include "exprtk.hpp"

using namespace std;


double *processFunctions(string functions[3], double estimatives[3]) {
    static double totals[3];

    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    symbol_table_t symbol_table;
    symbol_table.add_variable("x", estimatives[0]);
    symbol_table.add_variable("y", estimatives[1]);
    symbol_table.add_variable("z", estimatives[2]);
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;

    int i;

    for (i = 0; i < 3; i++) {
        parser.compile(functions[i], expression);
        totals[i] = expression.value();

    }

    return totals;
}


void processJacobian(double jacobian[3][3], double estimative[3]) {
    jacobian[0][0] = (2 * (estimative[0] + 1));
    jacobian[0][1] = (2 * (estimative[1] + 1));
    jacobian[0][2] = (2 * (estimative[2] + 1));

    jacobian[1][0] = (2 * (estimative[0]));
    jacobian[1][1] = (2 * (estimative[1]));
    jacobian[1][2] = (2 * (estimative[2]));

    jacobian[2][0] = 0;
    jacobian[2][1] = 1;
    jacobian[2][2] = -1;
}

int main() {

    double sum, fMax, xMax, velocity;
    int i, j, k, cont = 0, para = 0;
    double X[3], estimative[3], delta[3];
    string functions[3];
    double *functionResults;
    double jacobianResults[3][3];

    functions[0] = "(-1) * (((x+1)^2) + ((y+1)^2) + ((z+1)^2) - 12)";
    functions[1] = "(-1) * ((x^2) + (y^2) + (z^2) - 3)";
    functions[2] = "(-1) * (y-z)";

    X[0] = estimative[0] = 1.1;
    X[1] = estimative[1] = 0.903;
    X[2] = estimative[2] = 1.1;

    do {
        processJacobian(jacobianResults, estimative);
        functionResults = processFunctions(functions, estimative);

        for (i = 0; i < 3; i++) {
            for (k = i + 1; k < 3; k++) {
                sum = jacobianResults[k][i] / jacobianResults[i][i];
                jacobianResults[k][i] = 0;
                for (j = i + 1; j < 3; j++) {
                    jacobianResults[k][j] = jacobianResults[k][j] - (sum * jacobianResults[i][j]);
                }
                functionResults[k] = functionResults[k] - (sum * functionResults[i]);
            }
        }

        sum = 0;
        for (i = 2; i >= 0; i--) {

            sum = 0;
            for (k = (i + 1); k < 3; k++) {
                sum = sum + (jacobianResults[i][k] * delta[k]);
            }

            delta[i] = (functionResults[i] - sum) / jacobianResults[i][i];

        }

        for (i = 0; i < 3; i++) {
            X[i] = estimative[i];
            estimative[i] = estimative[i] + delta[i];
        }

        fMax = sqrt(pow(functionResults[0], 2) + pow(functionResults[1], 2) + pow(functionResults[2], 2));
        xMax = sqrt(pow(estimative[0] - X[0], 2) + pow(estimative[1] - X[1], 2) + pow(estimative[2] - X[2], 2));

        para++;
    } while (!(fMax < pow(10, -9) && xMax < pow(5, -8)));

    for (i = 0; i < 3; i++) {
        std::cout << "| " << estimative[i] << " |" << "\n";
    }

    std::cout << "Foram realizados " << para << " ciclos de iteração.";
}