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

// store the input text file of attributes
vector<string> attributes_input;
// store the input text file of queries
vector<string> queries_input;
// store the input text file of access frequency
vector<string> access_freq_input;

// the final list of attributes
vector<string> attributes_list;

// store and process the input of queries
vector<vector<string> > trimmed_queries;
vector<vector<string> > queries_word;

// store and process the input of access frequency
vector<string> firstTrim_acc;
vector<vector<string> > finalTrim_acc;
// the final sumup nuber of access frequency
vector<int> total_acc;

// store the attribute usage matrix
vector<vector<int> > attribue_usage_matrix;

// extract the aam column by column
vector<vector<int> > columns;

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

// return a list of attributes, index -> a(index + 1) e.g 0 = a1, 1 = a2...
vector<string> extract_element(vector<string> original_list)
{   
    vector<string> return_element;
    string space_delimiter = " ";
    int position;
    string attribute;

    // skip first line : "LABEL NAME"
    for (int i = 1; i < original_list.size(); i++)
    {
        position = original_list[i].find(space_delimiter);
        // Below commented code is for reading input file on local terminal, due to index difference
        //attribute = original_list[i].substr(position+1, original_list[i].length()-position-2);
        attribute = original_list[i].substr(position+1, original_list[i].length()-position-1);
        return_element.push_back(attribute);
    }
    
    return return_element;
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

// find SELECT, FROM, WHERE indexes
int findKeywordIndex (vector<string> wordList, string keyword)
{
    for (int i = 0; i < wordList.size(); i++)
    {
        if (wordList[i] == keyword)
        {
            return i;
        }
    }
    return -1;  
}

// find if there is "=" and it's index
int findEqualIndex (string str)
{
    int pos = 0;
    int index;
    while((index = str.find("=", pos)) != string::npos) {
        pos = index + 1; //new position is from next element of index
    }
    return pos - 1;
}

// find if there is "SUM("
int findSumIndex (string str)
{
    int pos = 0;
    int index;
    while((index = str.find("SUM(", pos)) != string::npos) {
        pos = index + 1; //new position is from next element of index
    }
    return pos - 1;
}

// function to trim spaces of a string
string trim_keyword(string keyword) 
{
    const char* typeOfWhitespaces = " \t\n\r\f\v";
    keyword.erase(keyword.find_last_not_of(typeOfWhitespaces) + 1);
    keyword.erase(0,keyword.find_first_not_of(typeOfWhitespaces));
    
    return keyword;
}

// find splitting strings to extract informations needed from queries
vector<string> pick_keyWords(vector<string> pre_query)
{   
    vector<string> post_query;
    int select_pos, from_pos, if_where_pos;
    string trimmed_sum, trimmed_equal;

    select_pos = findKeywordIndex(pre_query, "SELECT");
    from_pos = findKeywordIndex(pre_query, "FROM");
    if_where_pos = findKeywordIndex(pre_query, "WHERE");
    
    for (int i = select_pos + 1; i < from_pos; i++)
    {
        post_query.push_back(pre_query[i]);
    }
    
    if (if_where_pos != -1)
    {   
        for (int i = if_where_pos + 1; i < pre_query.size(); i++)
        {   
            post_query.push_back(pre_query[i]);
        }
    }

    for (int i = 0; i < post_query.size(); i++)
    {   
        string substr1;
        int tmp1 = findEqualIndex(post_query[i]);
        if (tmp1 != -1)
        {
            substr1 = post_query[i].substr(0, tmp1);
            post_query[i] = substr1;
        }
    }

    for (int i = 0; i < post_query.size(); i++)
    {   
        string substr2;
        int tmp2 = findSumIndex(post_query[i]);
        if (tmp2 != -1)
        {
            substr2 = post_query[i].substr(4, post_query[i].size() - 5);
            post_query[i] = substr2;
        }
    }

    for (int i = 0; i < post_query.size(); i++)
    {
        post_query[i] = trim_keyword(post_query[i]);
    }
    

    return post_query;
}

// call the extractor pick_keyWords() on every query
vector<vector<string> > extract_queriesWords(vector<vector<string> >original_list)
{
    vector<vector<string> > tmp_list;
    vector<string> tmp_str;
    for (int i = 0; i < original_list.size(); i++)
    {
       tmp_str = pick_keyWords(original_list[i]);
       tmp_list.push_back(tmp_str);
    }

    return tmp_list;
}

// initialize the attribute usage matrix with all 0
vector<vector<int> > create_attribute_usage(int qu_size, int att_size)
{   
    vector<vector<int> > att_usage;
    for (int i = 0; i < qu_size; i++)
    {   
        vector<int> tmp;
        for (int j = 0; j < att_size; j++)
        {
            tmp.push_back(0);
        }
        att_usage.push_back(tmp);
    }

    return att_usage;
}

// increment the index by 1 when match is found and index has value 0
void increment_usage(vector<string> query, int q_num)
{     
    for (int i = 0; i < query.size(); i++)
    {  
       for (int j = 0; j < attributes_list.size(); j++)
       {   
           if (query[i] == attributes_list[j] && attribue_usage_matrix[q_num][j] == 0)
           {    
               attribue_usage_matrix[q_num][j]++;
           } 
       }
    }
}

// call the above function
void update_attribute_usage()
{   
    for (int i = 0; i < queries_word.size(); i++)
    {   
        increment_usage(queries_word[i], i);
    }
}

// function to get column by column from attribute usage matrix
vector<vector<int> > get_column_of_attr_usage(vector<vector<int> > attr_uasge)
{   
    vector<vector<int> > columns;
    int column_num = attr_uasge[0].size();

    for (int i = 0; i < column_num; i++)
    {   
        vector<int> tmp;
        for (int j = 0; j < attr_uasge.size(); j++)
        {
            tmp.push_back(attr_uasge[j][i]);
        }
        columns.push_back(tmp);
    }

    return columns;
} 

// sum up row by row of ACC file
vector<int> calculate_acc(vector<vector<string> > finalTrim_acc)
{   
    vector<int> acc;
    for (int i = 0; i < finalTrim_acc.size(); i++)
    {   
        int tmp = 0;
        for (int j = 0; j < finalTrim_acc[i].size(); j++)
        {  
           int toInt;
           istringstream(finalTrim_acc[i][j]) >> toInt;
           tmp = tmp + toInt;
        }
        acc.push_back(tmp);
    }
    return acc;
}

// multiply columns with queries frequencies
void indexMultiplyIndex()
{
    for (int i = 0; i < columns.size(); i++)
    {   
        for (int j = 0; j < total_acc.size(); j++)
        {
            int res = columns[i][j]*total_acc[j];
            columns[i][j] = res;
        }
    }
}

// perform the upper part of the formula to attribute affinity
// column1[i] * column2[i]
int upperPartOfFormula(vector<int> list1, vector<int> list2)
{   
    int_least32_t final_value;
    vector<int> multiply_two_list;
    for (int i = 0; i < columns[0].size(); i++)
    {   
        multiply_two_list.push_back(list1[i]*list2[i]);
    }

    final_value = std::accumulate(multiply_two_list.begin(), multiply_two_list.end(), 0);
    
    return final_value;
}

// perform the lower part of the formula
// sum up all column values for every column, then sumup value of column1 times that of column2, sqrt the result
float lowerPartOfFormula(vector<int> list1, vector<int> list2)
{
    int sum_of_list1;
    int sum_of_list2;
    float final_value;

    sum_of_list1 = std::accumulate(list1.begin(), list1.end(), 0);
    sum_of_list2 = std::accumulate(list2.begin(), list2.end(), 0);

    final_value = sqrt(sum_of_list1*sum_of_list2);
    
    return final_value;
}

// calculate and convert to an AA Matrix
vector<vector<string >> calculateAAmatrix()
{   
    vector<vector<string >> AAmatrix;
    // initialize the matrix with all 0
    for (int i = 0; i < columns.size(); i++)
    {
        vector<string> vect(columns.size(), "0");
        AAmatrix.push_back(vect);
    }

    for (int i = 0; i < columns.size(); i++)
    {   
        for (int j = 0; j < columns.size(); j++)
        {   
            double value;
            int up = upperPartOfFormula(columns[i],columns[j]);

            float low = lowerPartOfFormula(columns[i], columns[j]);

            // handle zero divided by zero exception
            if (up == 0 || low == 0)
            {
                value = 0;
    
            } else 
            {
                value = up/low;
                // round up the result
                value = ceil(value);
            }   
            
            // eliminate decimal places
            string toString = to_string(value);
            int pos = toString.find(".");
            toString = toString.substr(0, pos);

            AAmatrix[i][j] = toString;
        }
    }

    return AAmatrix;
}

// function to output the AA matrix
void printAAmatrix(vector<vector<string >> AAmatrix)
{  
    for (int i = 0; i < AAmatrix.size(); i++)
    {   
        for (int j = 0; j < AAmatrix[i].size(); j++)
        {   
            printf("%s ", AAmatrix[i][j].c_str());
        }
        cout << "\n";
    }
}

int main(int argc, char **argv)
{   
    // get data from files
    string att_file = argv[1];
    string queries_file = argv[2];
    string access_freq_file = argv[3];

    // process attribute list input
    attributes_input = inputReader(att_file);
    attributes_list = extract_element(attributes_input);

    // process queries input
    queries_input = inputReader(queries_file);
    trimmed_queries = inputTrimmer(queries_input);
    queries_word = extract_queriesWords(trimmed_queries);

    // get attribute usage matrix
    int queriesNum = queries_word.size();
    int attributeNum = attributes_list.size();
    attribue_usage_matrix = create_attribute_usage(queriesNum, attributeNum);
    update_attribute_usage();
    // extract attribute usage column by column 
    columns = get_column_of_attr_usage(attribue_usage_matrix);
    
    // process access frequencies input
    access_freq_input = inputReader(access_freq_file);
    firstTrim_acc = extract_element(access_freq_input);
    finalTrim_acc = inputTrimmer(firstTrim_acc);
    total_acc = calculate_acc(finalTrim_acc);
    
    // multiply the column with queries frequency
    indexMultiplyIndex();
    // get the AA matrix
    vector<vector<string >> AAmatrix = calculateAAmatrix();
    // print the AA matrix
    printAAmatrix(AAmatrix);

}


