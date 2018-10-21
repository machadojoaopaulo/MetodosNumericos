#include <iostream>
#include <math.h>
#include <string.h>


using namespace std;

// classe para determinar a Funçao(x)
double *processFunctions(double estimative[3]) {
    static double totals[3];

    totals[0] = (-1) * (pow((estimative[0]+1),2) + pow((estimative[1]+1),2) + pow((estimative[2]+1), 2) - 12);
    totals[1] = (-1) * (pow(estimative[0], 2) + pow(estimative[1], 2) + pow(estimative[2], 2) - 3);
    totals[2] = (-1) * (estimative[1]-estimative[2]);

    return totals;
}

// classe para determinar a Jacobinada(x)
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

    FILE *fptr;
    double sum, fMax, xMax, velocity;
    int i, j, k, cont = 0, para = 0;
    double X[3], estimative[3], delta[3], deltaCount[3];
    string functions[3];
    double *functionResults;
    double jacobianResults[3][3];

    X[0] = estimative[0] = 1.1;
    X[1] = estimative[1] = 0.903;
    X[2] = estimative[2] = 1.1;


    //Abre o arquivo onde salva os erros
    fptr = fopen("/home/joaopaulo/CLionProjects/MetodosNumericos/erro.txt", "w");
    if (fptr == NULL) {
        printf("Arquivo nao existe! \n");
        return 0;
    }

    do {
        processJacobian(jacobianResults, estimative);
        functionResults = processFunctions(estimative);




        // Utilizando eleminaçao gaussiana
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

        // Determinado o delta
        for (i = 2; i >= 0; i--) {

            sum = 0;
            for (k = (i + 1); k < 3; k++) {
                sum = sum + (jacobianResults[i][k] * delta[k]);
            }

            delta[i] = (functionResults[i] - sum) / jacobianResults[i][i];

        }


        // Determinando a nova estimativa
        for (i = 0; i < 3; i++) {
            X[i] = estimative[i];
            estimative[i] = estimative[i] + delta[i];
        }
        // Determinando o criterio de parada
        fMax = sqrt(pow(functionResults[0], 2) + pow(functionResults[1], 2) + pow(functionResults[2], 2));
        xMax = sqrt(pow(estimative[0] - X[0], 2) + pow(estimative[1] - X[1], 2) + pow(estimative[2] - X[2], 2));



        // Determinar o erro e escreve-lo em um arquivo
        if(para >=1){
            velocity = sqrt(pow(delta[0],2)+pow(delta[1],1)+pow(delta[2],2)) / sqrt(pow(deltaCount[0],2)+pow(deltaCount[1],1)+pow(deltaCount[2],2));
            fprintf(fptr, "%lf \n", velocity);
            }

        // Um vetor recebe o delta antigo, para podermos calcular o erro
        for(i=0; i<3; i++){
            deltaCount[i] = delta[i];
        }
        para++;

        // Aplicando o criterio de parada
    } while (!(fMax < pow(10, -9) && xMax < pow(5, -8)));


    // imprimindo a soluçao
    for (i = 0; i < 3; i++) {
        std::cout << "| " << estimative[i] << " |" << "\n";
    }

    fclose(fptr);
}