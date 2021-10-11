#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm> 
#include <numeric>
#include <functional> 
#include <cctype>
#include <locale>
#include <math.h>

using namespace std ;

// read data from input txt. files
vector<string> inputReader(string input_txt) 
{
    vector<string> temp;

    ifstream input;
    std::ifstream file(input_txt);
    std::string str;
    while (std::getline(file, str)) {
        // remove all comma found
        str.erase(std::remove(str.begin(), str.end(), ','), str.end());
        temp.push_back(str);
    } 
    input.close();
    return temp;
}

// lines of query are first ead into a vector of string
// now assign the strings of each line into a vector
vector<vector<string> > inputTrimmer(vector<string> input_txt)
{   
    vector<vector<string> > trimmed_data; 
    for (int i = 0; i < input_txt.size(); i++)
    {
        vector<string> tmp_str;
        stringstream str_s(input_txt[i]);
        int j = 0;

        while(str_s.good())
        {
            string tmp;
            // split by space
            getline(str_s, tmp, ' ');
            tmp_str.push_back(tmp);
        }

        trimmed_data.push_back(tmp_str);
    }

    return trimmed_data;
}

// convert the data of AA matrix from string to integer
vector<vector<int> > aam_to_Int(vector<vector<string> > aamatrix)
{   
    vector<vector<int> > columns;
    int column_num = aamatrix.size();

    for (int i = 0; i < column_num; i++)
    {   
        vector<int> tmp;
        for (int j = 0; j < aamatrix[i].size()-1; j++)
        {   
            int toInt;
            istringstream(aamatrix[i][j]) >> toInt;
            tmp.push_back(toInt);
        }
        columns.push_back(tmp);
    }

    return columns;
} 

// function to initialize the ca matrix by adding the first two columns of AA matrix to CA matrix
void initialize_ca_matrix(vector<vector<int> > aa_Matrix, vector<vector<int> > ca_Matrix)
{
    ca_Matrix[0] = aa_Matrix[0];
    ca_Matrix[1] = aa_Matrix[1];
}

// function to calculate bond for two attributes
int calculateBond(vector<vector<int> > aa_Matrix, vector<int> ca_order, int first_attribute, int second_attribute, int middle_attribute_index)
{
    if (first_attribute < 0 || second_attribute < 0)
    {
        return 0;
    }

    if (first_attribute > middle_attribute_index || second_attribute > middle_attribute_index)
    {
        return 0;
    }

    int cont = 0;
    for (int i = 0; i < aa_Matrix[0].size(); i++)
    {   
        int tmp_res = aa_Matrix[i][ca_order[first_attribute]]*aa_Matrix[i][ca_order[second_attribute]];
        cont += tmp_res;
    }
    
    return cont;
}

// function to calculate cont for three attributes
int BEA_algorithm(vector<vector<int> > aa_Matrix, vector<int> ca_order, int first_attribute, int second_attribute, int third_attribute)
{
    int middle_attribute_index = second_attribute;

    int bond1 = 2*calculateBond(aa_Matrix, ca_order, first_attribute, second_attribute, middle_attribute_index);
    int bond2 = 2*calculateBond(aa_Matrix, ca_order, second_attribute, third_attribute, middle_attribute_index);
    int bond3 = 2*calculateBond(aa_Matrix, ca_order, first_attribute, third_attribute, middle_attribute_index);

    int res = bond1 + bond2 - bond3;
    return res;
}

// function to shuffle rows
vector<vector<int> > shuffle_rows(int att_num, vector<vector<int> > aa_Matrix, vector<vector<int> > ca_Matrix, vector<int> ca_order)
{
    for (int i = 0; i < att_num; i++)
    {
        for (int j = 0; j < att_num; j++)
        {
            ca_Matrix[i][j] = aa_Matrix[ca_order[i]][ca_order[j]];
        }
    }

    return ca_Matrix;
}

// function to print the CA matrix
void print_CA_matrix(int att_num, vector<vector<int> > ca_matrix, vector<int> ca_order)
{
    for (int i = 0; i < att_num; i++)
    {
        printf("%s%i\t", "A", ca_order[i]+1);
    }

    cout << "\n";

    for (int i = 0; i < ca_matrix.size(); i++)
    {
        printf("%s%i\t", "A", ca_order[i]+1);
        for (int j = 0; j < ca_matrix[i].size(); j++)
        {
            printf("%i ", ca_matrix[i][j]);
        }
        cout << "\n";
    }
}

// core function to convert CA matrix 
void convertCA_matrix(vector<vector<int> > aaMatrix)
{   
    // store the CA matrix
    vector<vector<int> > CA_matrix; 

    // initialize the CA matrix with all 0
    for (int i = 0; i < aaMatrix[0].size(); i++)
    {
        vector<int> vect(aaMatrix[0].size(), 0);
        CA_matrix.push_back(vect);
    }

    // initialize the order of attributes of ca matrix
    vector<int> ca_order;
    
    for (int i = 0; i < aaMatrix[0].size(); i++)
    {
        ca_order.push_back(i);
    }

    // first initialize the CA matrix
    initialize_ca_matrix(aaMatrix, CA_matrix);

    // initial index
    int start_column = 2;
    // choose the “best” location for attribute AA index
    while (start_column < aaMatrix[0].size())
    {   
        int BEA_res = 0;
        // placement given by maximum cont value
        int largestBEAindex = -1;
        int largest = -1;

        for (int i = 0; i < start_column; i++)
        {
            BEA_res = BEA_algorithm(aaMatrix, ca_order, i - 1, start_column, i);
            if (BEA_res > largest)
            {
                largest = BEA_res;
                largestBEAindex = i; 
            }
        }
        
        // calculate cont(index - 1; index; index+1)
        BEA_res = BEA_algorithm(aaMatrix, ca_order, start_column - 1, start_column, start_column + 1);

        if (BEA_res > largest)
        {
            largestBEAindex = start_column;
        }

        // for j from index to largestBEAindex by -1 , shuffle the two matrices
        for (int j = start_column; j > largestBEAindex; j--)
        {
            CA_matrix[j] = CA_matrix[j-1];
            ca_order[j] = ca_order[j-1];
        }
        
        CA_matrix[largestBEAindex] = aaMatrix[start_column];
        ca_order[largestBEAindex] = start_column;

        // increment the index number
        start_column++;
    }

    // shuffle rows
    CA_matrix = shuffle_rows(aaMatrix[0].size(), aaMatrix, CA_matrix, ca_order);
    print_CA_matrix(aaMatrix[0].size(), CA_matrix, ca_order);
    
}

// main
int main(int argc, char **argv)
{   
    // store the aa input
    vector<string> AAmatrix_input;
    vector<vector<string> > AAmatrix;
    vector<vector<int> > AAmInInt;

    // get data from files
    string aam_file = argv[1];
    
    // process attribute list input
    AAmatrix_input = inputReader(aam_file);
    AAmatrix = inputTrimmer(AAmatrix_input);
    AAmInInt = aam_to_Int(AAmatrix);

    // get the CA matrix
    convertCA_matrix(AAmInInt);
}


