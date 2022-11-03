#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <bits/stdc++.h>
#include <chrono>
#include <ctime>
#include <thread>

using namespace std;

vector<vector<float>> pearsonCorrelationPreprocess(vector<vector<float>> expression) {

    int NumberOfSamples = expression[0].size();

    float sum = 0, squareSum = 0;
    vector<vector<float>> res;

    for (int i = 0; i < expression.size(); i++) {
        res.push_back(vector<float>());

        for (int j = 0; j < NumberOfSamples; j++) {

            sum = sum + expression[i][j];
            squareSum = squareSum + expression[i][j] * expression[i][j];

        }

        res[i].push_back(sum);
        res[i].push_back(squareSum);
        sum = 0, squareSum = 0;
    }
    return res;
}

void pearsonCorrelationCalc(vector<vector<float>> expression, vector<vector<float>> expressionPreprocess, int start, int end,
                            vector<string> gene, int threadNum) {
    int NumberOfSamples = expression[0].size();

    string filename(
            "C:\\Users\\Hamed\\CLionProjects\\correlation\\Result\\result-thread-" + to_string(threadNum) + "-split0.csv");
    ofstream file_out;

    file_out.open(filename, std::ios_base::trunc);

    for (int i = start; i < end; i++) {

        ///First arg
        vector<float> firstGene(expression[i].begin(), expression[i].end());
        for (int j = i + 1; j < expression.size(); j++) {

            ///Second arg
            vector<float> secondGene(expression[j].begin(), expression[j].end());

            float sum_X = expressionPreprocess[i][0], sum_Y = expressionPreprocess[j][0], sum_XY = 0;
            float squareSum_X = expressionPreprocess[i][1], squareSum_Y = expressionPreprocess[j][1];

            for (int k = 0; k < NumberOfSamples; k++) {

                /// sum of X[i] * Y[i].
                sum_XY = sum_XY + firstGene[k] * secondGene[k];

            }

            /// use formula for calculating correlation coefficient.
            float corr = (float) (NumberOfSamples * sum_XY - sum_X * sum_Y)
                         / sqrt((NumberOfSamples * squareSum_X - sum_X * sum_X)
                                * (NumberOfSamples * squareSum_Y - sum_Y * sum_Y));

            file_out << gene[i] << "," << gene[j] << "," << corr << endl;
        }
    }
}

/*
void pearsonCorrelationCalc(vector<vector<float>> expression, int start, int end, vector<string> gene, int threadNum) {
    int NumberOfSamples = expression[0].size();

    string filename("C:\\Users\\Hamed\\CLionProjects\\correlation\\Result\\result-thread-" + to_string(threadNum) +
                    "-split0.csv");
    ofstream file_out;

    file_out.open(filename, std::ios_base::trunc);

    for (int i = start; i < end; i++) {

        ///First arg
        vector<float> firstGene(expression[i].begin(), expression[i].end());
        for (int j = i + 1; j < expression.size(); j++) {

            ///Second arg
            vector<float> secondGene(expression[j].begin(), expression[j].end());

            float sum_X = 0, sum_Y = 0, sum_XY = 0;
            float squareSum_X = 0, squareSum_Y = 0;

            for (int i = 0; i < NumberOfSamples; i++) {

                /// sum of elements of array X.
                sum_X = sum_X + firstGene[i];

                /// sum of elements of array Y.
                sum_Y = sum_Y + secondGene[i];

                /// sum of X[i] * Y[i].
                sum_XY = sum_XY + firstGene[i] * secondGene[i];

                /// sum of square of array elements.
                squareSum_X = squareSum_X + firstGene[i] * firstGene[i];
                squareSum_Y = squareSum_Y + secondGene[i] * secondGene[i];
            }

            /// use formula for calculating correlation coefficient.
            float corr = (float) (NumberOfSamples * sum_XY - sum_X * sum_Y) /
                         sqrt((NumberOfSamples * squareSum_X - sum_X * sum_X) *
                              (NumberOfSamples * squareSum_Y - sum_Y * sum_Y));

            file_out << gene[i] << "," << gene[j] << "," << corr << endl;
        }
    }
}
*/

void printExpression(vector<vector<float>> expression) {

    for (int i = 0; i < expression.size(); i++) {
        for (int j = 0; j < expression[i].size(); j++) {
            cout << expression[i][j] << "\n";
        }
    }

}

void printGenes(vector<string> gene) {

    for (int i = 0; i < gene.size(); i++) {
        cout << gene[i] << "\n";
    }

}

void printSamples(vector<string> sample) {

    for (int i = 0; i < sample.size(); i++) {
        cout << sample[i] << "\n";
    }

}

int main() {

    string fname = "C:\\Users\\Hamed\\CLionProjects\\correlation\\Split\\BeatAML_isoforms_0.csv";
    //string fname = "C:\\Users\\Hamed\\CLionProjects\\correlation\\Split\\BeatAML_isoforms.csv";

    vector<vector<float>> expression;

    vector<float> row;
    string line, word;
    bool isHeader = 1;
    bool isGen = 0;
    vector<string> gene;
    vector<string> sample;

    fstream file(fname, ios::in);
    if (file.is_open()) {
        while (getline(file, line)) {
            row.clear();

            stringstream str(line);

            while (getline(str, word, ',')) {
                if (isHeader) {
                    sample.push_back(word);
                } else if (isGen) {
                    gene.push_back(word);
                    isGen = 0;
                } else if (!isGen & !isHeader) {
                    row.push_back(stof(word));
                }
            }

            if (!isGen & !isHeader) {
                expression.push_back(row);
            }

            isGen = 1;
            isHeader = 0;
        }
    } else
        cout << "Could not open the file\n";
    vector<vector<float>> expressionPreprocess = pearsonCorrelationPreprocess(expression);

    auto start = chrono::system_clock::now();
    time_t start_time = chrono::system_clock::to_time_t(start);
    cout << "Computation started at: " << ctime(&start_time);


    int algorithm = 1;
    cout << "Choose your desired algorithm for correlation: \n 1. Pearson correlation \n 2. "
            "Kendall rank correlation \n 3. Spearman correlation \n 4. Point-Biserial correlation \n";
    unsigned int numThreads = std::thread::hardware_concurrency();
    cout << "number of recommended threads: " << numThreads << "\n";

    int threadCount = 0;

    /*int algorithm;
    cin >> algorithm;*/

    ///Defining a vector of pairs: thread and callback value

    switch (algorithm) {
        case 1  : {
            cout << "You chose Pearson correlation \n";

            cout << "number of rows:" << expression.size();

            /*int a;
            cin >> a;*/

            thread myThreads[numThreads];

            //vector<int> fair = {0, 1500, 3000, 4500, 9000, 12000, 18000, 24000, 31266};
            vector<int> fair = {0, 50, 100, 150, 300, 400, 600, 800, 999};

            for (int i = 0; i < numThreads; i++) {
                myThreads[i] = thread(pearsonCorrelationCalc, expression, expressionPreprocess, fair[i], fair[i + 1], gene, i);
            }
            for (int i = 0; i < numThreads; i++) {
                myThreads[i].join();
            }
        }
            break; //optional
        case 2  :
            cout << "You chose Kendall rank correlation";
            break; //optional
        case 3  :
            cout << "You chose Spearman correlation";
            break; //optional
        case 4  :
            cout << "You chose Point-Biserial correlation";
            break; //optional

            // you can have any number of case statements.
        default : //Optional
            cout << "invalid input";
    }

    auto end = chrono::system_clock::now();
    time_t end_time = chrono::system_clock::to_time_t(end);
    cout << "Computation finished at: " << ctime(&end_time);
    return 0;
}

