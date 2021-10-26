#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <sstream>
#include <math.h>
#include <cmath>
#include <algorithm>

using namespace std;

// ID of medoids
const string medoids_constID = "medoidID";
// ID of nearest point  
const string nearest_pointID = "nearestPoint";

// a class of a packet
class packet
{
public:
    string flow_tag;
    int arrival_time;
    int packet_size;
};

// a class of the cluster
class cluster
{
public:
    int mediodIndex;
    vector<int> nearest_point;

    void add_point(int point_ID)
    {
        nearest_point.push_back(point_ID);
    }

    void delete_points()
    {
        nearest_point.clear();
    }

};

//calculate distance between 2 nodes
double calculate_nodeDistance(vector<vector<double>> flow_data, int p_1, int p_2)
{
    return abs(flow_data[p_1][0] - flow_data[p_2][0]) + abs(flow_data[p_1][1] - flow_data[p_2][1]);
}

// function to get the nearest mediod
int get_nearestMedoid(int point_ID,vector<vector<double>> distance_matrix, vector<int> medoid_ID)
{
    int nearest_med;
    double minimum_distance = 9999999.0;
    double distance;

    for(int i = 0; i < medoid_ID.size(); i++)
    {
        distance = distance_matrix[point_ID][medoid_ID[i]];
        if (distance < minimum_distance)
        {
            minimum_distance = distance;
            nearest_med = medoid_ID[i];
        }
    }
    return nearest_med;
}

// check if a point is a medoid
bool checkIfMedoid(vector<int> medoid_index, int key)
{
    for (int i = 0; i < medoid_index.size(); i++) {

        if (key == medoid_index[i])
        {
            return true;
        }
    }
    return false;
}

// generate distance matrix
vector<vector<double>> get_distanceMatrix(vector<vector<double>> flow_data, int flow_size)
{
	
    vector<vector<double> > matrix;
    for (int i = 0; i < flow_size; i++) 
    {
        vector<double> tmp_distance;

        for (int j = 0; j < flow_size; j++) 
        {
            tmp_distance.push_back(calculate_nodeDistance(flow_data, i, j));
        }

        matrix.push_back(tmp_distance);
    }

    return matrix;
}

// implementing k-medoid clustering algorithm 
void assign_pointInCluster(vector<vector<double>>& flow_data, map<string, vector<int>>& current_medoids, double& total_distance, int total_medoid, int flow_size, double& sum_square_of_error)
{   
    // store ID of KMedoids
	vector<int> medoid_ID;

    for (int i = 0; i < total_medoid; i++)
    {
        medoid_ID.push_back(current_medoids[medoids_constID][i]); // medoids index
    }

    // arrange medoids into cluster
    vector<vector<double>> distance;
    distance = get_distanceMatrix(flow_data, flow_size); // minimum distance for each row of medoids

    //store the nearest medoids for evert point
    vector<int> nearest_MedoidsID;
    for(int i = 0; i < flow_size; i++)
    {
    	int nearest_MedoidID;
    	nearest_MedoidID = get_nearestMedoid(i, distance, medoid_ID);
    	nearest_MedoidsID.push_back(nearest_MedoidID);
    }

    current_medoids[nearest_pointID] = nearest_MedoidsID;

    // calculate the total distance of all clusters
    vector<double> total_minimum_distance;
    for (int i = 0; i < flow_size; i++)
    {
        double point_dist;
        point_dist = calculate_nodeDistance(flow_data, i, nearest_MedoidsID[i]);
        total_minimum_distance.push_back(point_dist);
    }

    for (int i = 0; i < flow_size; i++)
    {
        total_distance += total_minimum_distance[i];
        sum_square_of_error += pow(total_minimum_distance[i],2);
    }

    sum_square_of_error = sum_square_of_error / flow_size;
    
}

// find in the flow
vector<int> search_flow(string key, vector<packet> packets)
{
    vector<int> target_flow_ID;

    for (int i = 0; i < packets.size(); i++) 
    {
        if (packets[i].flow_tag == key)
        {
            target_flow_ID.push_back(i);
        }
    }
    return target_flow_ID;
}

// remove flow
void remove_flow(vector<int> target_flow_ID, vector<packet>& packets)
{
    vector <packet>::iterator Iterator;

    // remove from end to begin of packet vector
    for (int i = target_flow_ID.size()-1; i >= 0; i--) 
    {
        Iterator = packets.begin()+target_flow_ID[i];
        packets.erase(Iterator);
        Iterator = packets.begin();
    }
}

// transfer vector to set
set<int> transfer_vect_to_set(vector<int> input_vector)
{
    set<int> tmp_set;
    // foreach of vector input
    for (int i : input_vector)
    {
        tmp_set.insert(i);
    }
    return tmp_set;
}

// copy map data from other map data
void copy_map(map<string, vector<int>>& sourceMap, map<string, vector<int>>& newMap, vector<string> labels)
{
    for (string key : labels)
    {
        newMap[key] = sourceMap[key];
    }
}

double k_medoids_cluster(map<string, vector<int>>& current_medoids, vector<vector<double>>& flow_data, int flow_size, int total_medoid, vector<int> medoid_index)
{
	map<string, vector<int>> initial_medoids_map;
    set<int> current_medoids_set;
    set<int> initial_medoids;
    vector<int> initial_medoidsID;
    vector<string> labels;
    double total_distance = 0.0;
    double sum_square_error = 0.0;

    current_medoids[medoids_constID] = medoid_index;
    assign_pointInCluster(flow_data, current_medoids, total_distance, total_medoid, flow_size, sum_square_error);

    initial_medoids_map[medoids_constID] = initial_medoidsID;
    current_medoids_set = transfer_vect_to_set(current_medoids[medoids_constID]);
    initial_medoids = transfer_vect_to_set(initial_medoids_map[medoids_constID]);

    map<string, vector<int>>::iterator iterator;

    for (iterator = begin(current_medoids); iterator != end(current_medoids); iterator++)
    {
        labels.push_back(iterator->first);
    }

    //if medoids updated
    while (current_medoids_set != initial_medoids)
    {
        map<string, vector<int>> closest_medoids;
        // from current_medoids to closest_medoids
        copy_map(current_medoids, closest_medoids, labels);
        copy_map(current_medoids, initial_medoids_map, labels);

        map<string, vector<int>> tmp_medoids_map;
        copy_map(current_medoids, tmp_medoids_map, labels);

        for (int i = 0; i < flow_size; i++)
        {
            for (int j = 0; j < total_medoid; j++)
            {
                //if i is not a current medoid
                if (!checkIfMedoid(closest_medoids[medoids_constID], i))
                {
                    double tmp_total_cost = 0.0;
                    double tmp_sum_square_error = 0.0;
                    
                    // copy to tmp_medoids_map
                    copy_map(closest_medoids, tmp_medoids_map, labels);
                    // set i as medoid
                    tmp_medoids_map[medoids_constID][j] = i;
                    
                    assign_pointInCluster(flow_data, tmp_medoids_map, tmp_total_cost, total_medoid, flow_size, tmp_sum_square_error);//calculate temp_cost

                    // check if swapping improves the total distance of the resulting clustering
                    if (total_distance - tmp_total_cost >0 )
                    {
                        total_distance = tmp_total_cost;
                        sum_square_error = tmp_sum_square_error;
                        copy_map(tmp_medoids_map, closest_medoids, labels);
                        i = 0;
        				j = 0;
                    }                    
                }
            }
        }
        
        // copy the updated medoids to current_medoids
        copy_map(closest_medoids, current_medoids, labels);
        current_medoids_set = transfer_vect_to_set(current_medoids[medoids_constID]);
        initial_medoids = transfer_vect_to_set(initial_medoids_map[medoids_constID]);
    }
    return total_distance;
}

// convert string to integer
int strToInt(string str)
{
  int tmp;

  stringstream ss;  
  ss << str;  
  ss >> tmp;  

  return tmp;
}

int main(int argc, const char * argv[]) 
{
    string file_one = argv[1];
    string file_two = argv[2];
    
    // vector of packets
    vector<packet> total_packets;
    vector<vector<double>> flow_data;
    vector<int> flow_index;
    
    // start reading file1.txt
    ifstream input_file1;
    string data;
    input_file1.open(file_one,ios::in);
    getline(input_file1,data);

    while (getline(input_file1, data)) 
    {
        string tmp;
        // read data line by line
        vector<string> factors;
        stringstream ss(data);

        while (ss>>tmp) 
        {
            factors.push_back(tmp);
        }

        packet packet;
        packet.flow_tag = factors[0]+factors[1]+factors[2]+factors[3]+factors[4];
        packet.arrival_time = strToInt(factors[5]);
        packet.packet_size = strToInt(factors[6]);
        total_packets.push_back(packet);
    }
    input_file1.close();
    
    while (total_packets.size() != 0) 
    {
        flow_index = search_flow(total_packets[0].flow_tag,total_packets);

        // remove flow with only one packet
        if (flow_index.size() <= 1) 
        {
            remove_flow(flow_index, total_packets);
        }
        else
        {
            double time = 0;
            double length = 0;
            // save output 
            vector<double> flow_result;
            for (int i = flow_index.size()-1; i > 0; i--) 
            {
                length = length + total_packets[flow_index[i]].packet_size;
                time = time + total_packets[flow_index[i]].arrival_time-total_packets[flow_index[i-1]].arrival_time;
            }
         
            length = length + total_packets[0].packet_size;
            time = time/(flow_index.size()-1);
            length = length/flow_index.size();

            flow_result.push_back(time);
            flow_result.push_back(length);
            flow_data.push_back(flow_result);
            remove_flow(flow_index, total_packets);
        }
    }

    // output to Flow.txt
    ofstream flow_output;
    flow_output.open ("Flow.txt", ios::out | ios::trunc);
    flow_output.flags(ios::fixed);

    for(int i = 0; i < flow_data.size(); i++)
    {
        flow_output << i << " ";
        flow_output << setprecision(2) << flow_data[i][0];
        flow_output << " ";
        flow_output << setprecision(2) << flow_data[i][1] << endl;
    }
    
    flow_output.close();

    int flow_size = flow_data.size();
    int total_medoid;
    vector<int> medoid_index;
    
    ifstream input_file2;
    
    // read file2.txt
    // string file_two="file2.txt";
    input_file2.open(file_two,ios::in);
    getline(input_file2,data);
    // number of medoids
    total_medoid = strToInt(data);

    while (getline(input_file2, data)) 
    {
        int num;
        stringstream ss(data);
        while (ss >> num) 
        {
        	medoid_index.push_back(num);
        }
    }

    input_file2.close();

    // run K-Medoids algorithm
    map<string, vector<int>> current_medoids;
    double total_distance = k_medoids_cluster(current_medoids, flow_data, flow_size, total_medoid, medoid_index);

    ofstream KMedoids_output;

    KMedoids_output.open ("KMedoids.txt",ios::out | ios::trunc);
    KMedoids_output.flags(ios::fixed);
    
    KMedoids_output << setprecision(2) << total_distance << endl;

    vector<cluster> listOfClusters;
    for(int i = 0; i < total_medoid; i++)
    {
        cluster cluster;
        listOfClusters.push_back(cluster);
    }

    for (int i = 0; i < listOfClusters.size(); i++) 
    {
        listOfClusters[i].delete_points();
        listOfClusters[i].mediodIndex = current_medoids[medoids_constID][i];
        for (int j = 0; j < current_medoids[nearest_pointID].size(); j++)
        {
            if (current_medoids[nearest_pointID][j] == listOfClusters[i].mediodIndex)
            {
                listOfClusters[i].add_point(j);
            }
        }
       
    }
    
    for (int i : current_medoids[medoids_constID])
    {
        KMedoids_output << i << " ";
    }

    KMedoids_output << endl;

    for (int i = 0; i < listOfClusters.size(); i++) 
    {
        for (int j = 0; j < listOfClusters[i].nearest_point.size(); j++) {
            // Sets the decimal precision to be used to format floating-point values
            KMedoids_output << setprecision(2) << listOfClusters[i].nearest_point[j] << " ";
        }
        KMedoids_output << endl;
    }

    KMedoids_output.close();
    
    return 0;
}
