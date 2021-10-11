#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>


using namespace std;

bool debug = true;

vector<string> AAmatrix_input;
vector<vector<string> > AAmatrix;
vector<vector<int> > AAmInInt;

vector<vector<int> > CAmatrix;

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

void addFirstColumnToCA()
{
    CAmatrix.push_back(AAmInInt[0]);
}

int bond (int attribute1, int attribute2, int index, vector<vector<int> > affinity, vector<int> order) 
{
    if (attribute1 < 0 || attribute2 < 0) 
    {
        return 0;
    }
    if (attribute1 > index || attribute2 > index) 
    {
        return 0;
    }

    int sum = 0;

    for (int z = 0; z < affinity[0].size(); z++) 
    {
        sum += affinity[z][order[attribute1]] * affinity[z][order[attribute2]];
    }

    return sum;
}

int contribution (int attribute1, int attribute2, int attribute3, vector<vector<int> > affinity, vector<int> order) 
{
    int index = attribute2;
    return 2*bond(attribute1, attribute2, index, affinity, order)
    + 2*bond(attribute2, attribute3, index, affinity, order)
    - 2*bond(attribute1, attribute3, index, affinity, order);
}

void initialize_ca_matrix(vector<vector<int> > aa_Matrix, vector<vector<int> > ca_Matrix)
{
    ca_Matrix[0] = aa_Matrix[0];
    ca_Matrix[1] = aa_Matrix[1];
}

vector<vector<int> > calculateClusteredAffinity(vector<vector<int> > affinity) 
{
    int no_attributes = affinity[0].size();
    vector<vector<int> > ca; 
    for (int i = 0; i < no_attributes; i++)
    {
        vector<int> vect(no_attributes, 0);
        ca.push_back(vect);
    }

    vector<int> order; // Order of attributes
    for (int i = 0; i < no_attributes; i++) {
        order.push_back(i);
    }
    
    initialize_ca_matrix(affinity, ca);

    int index = 2;

    while (index < no_attributes) 
    { //Choose the 'best' location for attribute affinity[index]
        int max = -1;
        int loc = -1;
        int cont = 0;

        for (int i = 0; i < index; i++) {
            cont = contribution(i-1, index, i, affinity, order);
            if (cont > max) 
            {
                max = cont;
                loc = i;
            }
        }

        cont = contribution(index - 1, index, index + 1, affinity, order);
        if (cont > max) 
        {
            loc = index;
        }

        // Reshuffle Matrix
        for (int j = index; j > loc; j--) 
        {
            ca[j] = ca[j-1];
            order[j] = order[j-1];
        }
        ca[loc] = affinity[index];
        order[loc] = index;

        if (debug) 
        {
            cout << "Order: ";
            for (int i = 0; i < order.size(); i++) 
            {
                cout << order[i] << ' ';
            }
            cout << endl;
        }

        index++;
    }

    if (debug) 
    {
        cout << "Order: ";
        for (int i = 0; i < order.size(); i++) 
        {
            cout << order[i] << ' ';
        }
        cout << endl;

        for (int i = 0; i < ca.size(); i++) 
        {
            for (int j = 0; j < ca[i].size(); j++) 
            {
                cout << ca[i][j] << ' ';
            }
            cout << endl;
        }
        cout << endl;
    }

    // Shuffle columns
    for (int i = 0; i < no_attributes; i++) 
    {
        for (int j = 0; j < no_attributes; j++) 
        {
            ca[i][j] = affinity[order[i]][order[j]];
        }
    }

    for (int i = 0; i < no_attributes; i++) 
    {
        cout << "A" << order[i] + 1 << "\t";
    }
    cout << endl;

    // Print Clustered Affinity Matrix
    for (int i = 0; i < ca.size(); i++) 
    {
        cout << "A" << order[i] + 1 << "\t";
        for (int j = 0; j < ca[i].size(); j++) {
            cout << ca[i][j] << ' ';
        }
        cout << endl;
    }

    return ca;
}

/**
* Run the Bond Energy Algorithm to calculate the clustered affinity matrix CA of an Affinity Matrix
* @param argc The number of arguments supplied
* @param argv The arguments
* @return 0 if the program is successful.
*/
int main (int argc, char *argv[]) 
{

    // get data from files
    string aam_file = argv[1];
    
    // process attribute list input
    AAmatrix_input = inputReader(aam_file);
    AAmatrix = inputTrimmer(AAmatrix_input);
    AAmInInt = aam_to_Int(AAmatrix);

    // vector<vector<int> > affinity = AAmInInt;

    if (debug) 
    { //print affinity
        for (int i = 0; i < AAmInInt.size(); i++) 
        {
            for (int j = 0; j < AAmInInt[i].size(); j++) 
            {
                cout << AAmInInt[i][j] << ' ';
            }
            cout << endl;
        }
    }

    vector<vector<int> > ca = calculateClusteredAffinity(AAmInInt);
}

