#include <iostream>
#include <algorithm>
#include <cmath>
#include "TreeIndex.h"
// #include "../assignment3/VectorDataset.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <random>

using namespace std;

const int MINSIZE = 2;
DataVector::DataVector(int dimension)
{
    setDimension(dimension);
}

DataVector::DataVector(const std::vector<double> &vec)
{
    v = vec;
}

void DataVector::setDimension(int dimension)
{
    v.clear();
    v.resize(dimension, 0.0);
}

DataVector::~DataVector()
{
}

DataVector::DataVector(const DataVector &other)
{
    v = other.v;
}

DataVector &DataVector::operator=(const DataVector &other)
{
    if (this != &other)
    {
        v = other.v;
    }
    return *this;
}

DataVector DataVector::operator+(const DataVector &other)
{
    if (v.size() != other.v.size())
    {
        std::cerr << "Error: Vector dimensions are not equal for addition.\n";
        exit(EXIT_FAILURE);
    }

    DataVector sum(v.size());
    for (int i = 0; i < v.size(); i++)
    {
        sum.v[i] = v[i] + other.v[i];
    }
    return sum;
}

DataVector DataVector::operator-(const DataVector &other)
{
    if (v.size() != other.v.size())
    {
        std::cerr << "Error: Vector dimensions are not equal for subtraction.\n";
        exit(EXIT_FAILURE);
    }

    DataVector diff(v.size());
    for (int i = 0; i < v.size(); i++)
    {
        diff.v[i] = v[i] - other.v[i];
    }
    return diff;
}

double DataVector::operator*(const DataVector &other) const
{
    if (v.size() != other.v.size())
    {
        std::cerr << "Error: Vector dimensions are not equal for dot product.\n";
        exit(EXIT_FAILURE);
    }

    double dotProduct = 0.0;
    for (int i = 0; i < v.size(); i++)
    {
        dotProduct += v[i] * other.v[i];
    }
    return dotProduct;
}

DataVector &DataVector::operator=(const vector<double> &vec)
{
    v = vec;
    return *this;
}

bool DataVector::operator==(const DataVector &other) const
{
    if (v.size() != other.v.size())
    {
        std::cerr << "Error: Vector dimensions are not equal for dot product.\n";
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] != other.v[i])
        {
            return false;
        }
    }
    return true;

}

double DataVector::norm()
{
    return sqrt((*this) * (*this));
}

double DataVector::dist(const DataVector &other)
{
    DataVector diff = *this - other;
    return diff.norm();
}
void DataVector::print()
{
    for (int i = 0; i < v.size(); i++)
    {
        cout << v[i] << " ";
    }
    cout << endl;
}

vector<DataVector> readData(const string &filepath)
{
    vector<DataVector> final; // Create a VectorDataset object to store the final dataset
    ifstream file(filepath);  // Open the file for reading

    // Check if the file is successfully opened
    if (!file.is_open())
    {
        cerr << "Error: File not found\n";
        return final; // Return an empty dataset if the file is not found
    }

    string line;
    getline(file, line); // Read the first line (assuming it contains header information)

    // Read each line of the file
    while (getline(file, line))
    {
        istringstream iss(line);
        vector<double> values;
        double val;
        char comma;

        // Parse the line and extract values separated by commas
        while (iss >> val)
        {
            values.push_back(val);
            // Check for a comma after each value
            if (iss >> comma && comma != ',')
            {
                cerr << "Error: Unexpected character after value in line: " << line << endl;
                return final; // Return an empty dataset in case of unexpected characters
            }
        }

        // Create a DataVector and populate it with the parsed values
        DataVector temp;
        temp = values;
        DataVector dataVector(temp);

        // Add the DataVector to the final dataset
        final.push_back(dataVector);
    }

    file.close(); // Close the file after reading
    return final; // Return the final dataset
}

struct Node
{
    vector<DataVector> v;
    Node* left;
    Node* right;
    double medianval;
    vector<double> axis;
};

class RPTreeIndex : public TreeIndex
{
    Node* root;

    // generate random dimension
    vector<double> randomUnitDirection(size_t dimensions) {
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<double> dis(-1.0, 1.0);

        vector<double> direction(dimensions, 0.0);
        double sum = 0.0;
        for (size_t i = 0; i < dimensions; ++i) {
            direction[i] = dis(gen);
            sum += direction[i] * direction[i];
        }
        double norm = sqrt(sum);
        for (size_t i = 0; i < dimensions; ++i) {
            direction[i] /= norm;
        }
        return direction;
    }

    int randomX(int k){
        return rand()%k;
    }

    // Choose Rule
    // int chooseRule(vector<DataVector>::iterator begin, vector<DataVector>::iterator end)
    // {
    //     int axis = 0;
    //     // calculate variance of each feature and return the feature with maximum variance
    //     for (int i = 0; i < begin->v.size(); i++)
    //     {
    //         //find mean
    //         double mean = 0;
    //         for (auto it = begin; it != end; it++)
    //         {
    //             mean += (*it).v[i];
    //         }
    //         mean /= end - begin;
    //         double variance = 0, newvariance = 0;
    //         for (auto it = begin; it != end; it++)
    //         {
    //             newvariance += pow((*it).v[i] - mean, 2);
    //         }
    //         if (newvariance > variance) {
    //             variance = newvariance;
    //             axis = i;
    //         }
    //     }
    //     return axis;
    // }
    // initial call to build tree
    Node* buildTree(vector<DataVector>& points)
    {
        if (points.empty())
        {
            return NULL;
        }

        return buildTree(points.begin(), points.end());
    }

    // Overloaded buildTree function for iterators (actual implementation)
    Node *buildTree(vector<DataVector>::iterator begin, vector<DataVector>::iterator end)
    {
        int k = begin->v.size();
        //Here k is dimension
        // int axis = chooseRule(begin, end);
        vector<double> axis = randomUnitDirection(k);

        if (end - begin <= MINSIZE)
        {
            Node *newnode = new Node;
            newnode->axis = axis;
            newnode->left = NULL;
            newnode->right = NULL;
            newnode->medianval = -1;
            newnode->v = vector<DataVector>(begin, end); // Store the points in leaf nodes
            return newnode;
        }

        DataVector x = *( begin + randomX(end - begin) );

        DataVector y;
        double maxdist = 0;

        for (auto it = begin; it != end; it++)
        {
            double dist = 0;
            dist = (*it)*x;
            if(dist > maxdist){
                maxdist = dist;
                y = *it;
            }
        }

        
        Node *newnode = new Node;
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<double> dis(-1.0, 1.0);
        double delta = dis(gen)*6*(x.dist(y))/sqrt(k);

        // sort(begin, end, [axis](const DataVector &a, const DataVector &b)
        //      { return a.v[axis] < b.v[axis]; });

        newnode->axis = axis;
        // newnode->medianval = medianiter->v[axis];
        // newnode->v contains all points from begin to end

        sort(begin, end, [axis](const DataVector &a, const DataVector &b)
             { return a*axis < b*axis; });
             
        auto medianiter = begin + (end - begin) / 2;

        newnode->medianval = *medianiter*axis + delta;
        for(auto it = begin; it != end; it++){
            newnode->v.push_back(*it);
        }
        // No need to copy points to left and right vectors, just pass the iterators
        newnode->left = buildTree(begin, medianiter + 1);
        newnode->right = buildTree(medianiter + 1, end);

        return newnode;
    }
    vector<DataVector> search(DataVector &point, Node* node, int k=1)
    {
        if (node == NULL)
        {
            return vector<DataVector>();
        }
        if (node->medianval == -1)
        {
            return node->v;
        }

        double compareval = point*(DataVector(node->axis));
        vector<DataVector> nneighbours, sibling;
        Node *temp;
        // if(compareval==node->medianval){
        //     return node->v;
        // }
        if (compareval <= node->medianval)
        {
            temp = node->left;
            nneighbours = search(point, node->left,k);
            sibling = node->right->v;
        }
        else
        {
            temp = node->right;
            nneighbours = search(point, node->right,k);
            sibling = node->left->v;
        }

        //if the distance of the given point from the farthest point in the current subtree is less than the perpendicular distance of the given point from the median, then return the left subtree else return the current node
        double maxdist = -1;
        for(auto it = nneighbours.begin(); it != nneighbours.end(); it++)
        {
            if (it->dist(point) > maxdist)
            {
                maxdist = it->dist(point);
            }
        }
        double mediandist = abs(node->medianval - point*(DataVector(temp->axis)));

        if (maxdist > mediandist || nneighbours.size() < k) {
            for (auto it = sibling.begin(); it != sibling.end(); it++)
            {
                nneighbours.push_back(*it);

                
            }
            // for(auto it = nneighbours.begin(); it != nneighbours.end(); it++){
            //     (*it).print();
            // }
        }
        return nneighbours;

    }
    RPTreeIndex() : root(nullptr) {}
    static RPTreeIndex *instance;

public:
    static RPTreeIndex *GetInstance()
    {
        if (!instance)
        {
            instance = new RPTreeIndex();
        }
        return instance;
    }

    vector<DataVector> query_search(DataVector &point, int k)
    {
         return search(point, root, k);
    }

    void maketree(vector<DataVector> &points)
    {
        root = buildTree(points);
    }

    //AddData
    void AddData(DataVector &newpoint, vector<DataVector> &points)
    {
        points.push_back(newpoint);
        maketree(points);
    }

    //DeleteData
    void DeleteData(DataVector &newpoint, vector<DataVector> &points)
    {
        points.erase(remove(points.begin(), points.end(), newpoint), points.end());

        clearTree(root);

        maketree(points);
    }

    void clearTree(Node *node)
    {
        if (node == nullptr)
            return;

        clearTree(node->left);
        clearTree(node->right);
        delete node;
    }
};
RPTreeIndex *RPTreeIndex::instance = nullptr;

bool cmp(const pair<DataVector, double> &a, const pair<DataVector, double> &b)
{
    return a.second < b.second;
}
void knearestneighbor(vector<DataVector> &traindata, DataVector &testvector, int k)
{
    // Calculate distances between the test vector and all vectors in the training dataset
    vector<pair<DataVector, double>> distances;
    int n = traindata.size();
    for (int i = 0; i < n; i++)
    {
        double distance = testvector.dist(traindata[i]);
        distances.push_back(make_pair(traindata[i], distance));
    }

    // Sort the distances in ascending order
    sort(distances.begin(), distances.end(), cmp);

    // Create a dataset for the k-nearest neighbors
    for (int i = 0; i < min(k, n); i++)
    {
        cout << "Vector " << i + 1 << ": " << "Distance: " << distances[i].second << endl;
        distances[i].first.print();
    }
}

int main()
{
    // cout << "Reading Data initiated" << endl;
    vector<DataVector> myData = readData("fmnist-test.csv");
    // cout << "Reading data finished" << endl;
    // myData[0].print();
    // cout << myData.size() << endl;

    RPTreeIndex::GetInstance()->maketree(myData);
    //cout << "hello" << endl;
    DataVector testvector;
    int k = 9;
    // vector<double> v={0.0,0.0,0.0,2.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,57.0,67.0,73.0,76.0,76.0,83.0,62.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0,13.0,79.0,128.0,201.0,162.0,161.0,173.0,192.0,172.0,181.0,184.0,108.0,30.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,102.0,139.0,142.0,103.0,115.0,162.0,154.0,165.0,153.0,139.0,129.0,150.0,138.0,171.0,161.0,26.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,106.0,124.0,108.0,103.0,106.0,93.0,100.0,180.0,156.0,147.0,138.0,85.0,157.0,114.0,124.0,154.0,157.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.0,0.0,29.0,134.0,101.0,113.0,108.0,97.0,116.0,81.0,146.0,183.0,164.0,111.0,146.0,131.0,122.0,132.0,145.0,169.0,93.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,65.0,131.0,105.0,101.0,108.0,100.0,104.0,97.0,74.0,206.0,174.0,115.0,150.0,108.0,119.0,146.0,152.0,162.0,141.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,108.0,123.0,112.0,106.0,98.0,100.0,101.0,111.0,72.0,136.0,132.0,112.0,115.0,109.0,142.0,150.0,160.0,162.0,195.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,135.0,122.0,119.0,112.0,83.0,93.0,97.0,106.0,115.0,91.0,109.0,83.0,109.0,125.0,165.0,140.0,193.0,160.0,176.0,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,23.0,135.0,132.0,149.0,104.0,78.0,96.0,97.0,103.0,108.0,108.0,115.0,84.0,114.0,149.0,158.0,147.0,209.0,160.0,178.0,52.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,57.0,124.0,131.0,160.0,103.0,96.0,101.0,109.0,109.0,111.0,108.0,114.0,106.0,113.0,156.0,163.0,156.0,196.0,167.0,174.0,85.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,86.0,116.0,126.0,160.0,132.0,102.0,96.0,104.0,107.0,98.0,93.0,123.0,116.0,112.0,149.0,160.0,181.0,186.0,162.0,162.0,123.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,112.0,114.0,112.0,167.0,146.0,100.0,100.0,101.0,107.0,96.0,100.0,126.0,103.0,120.0,141.0,158.0,167.0,187.0,147.0,148.0,170.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,124.0,118.0,101.0,197.0,89.0,101.0,100.0,96.0,108.0,103.0,108.0,122.0,107.0,127.0,139.0,150.0,119.0,196.0,145.0,142.0,179.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,131.0,115.0,97.0,198.0,47.0,108.0,96.0,96.0,106.0,119.0,102.0,120.0,103.0,123.0,140.0,150.0,81.0,206.0,146.0,135.0,194.0,0.0,0.0,0.0,0.0,0.0,0.0,4.0,136.0,118.0,98.0,193.0,21.0,109.0,96.0,106.0,102.0,129.0,100.0,119.0,102.0,118.0,141.0,153.0,45.0,216.0,148.0,139.0,205.0,12.0,0.0,0.0,0.0,0.0,0.0,24.0,137.0,117.0,101.0,187.0,25.0,113.0,92.0,112.0,94.0,120.0,105.0,127.0,97.0,115.0,142.0,150.0,19.0,213.0,148.0,131.0,204.0,36.0,0.0,0.0,0.0,0.0,0.0,40.0,141.0,109.0,117.0,158.0,12.0,124.0,90.0,116.0,96.0,117.0,113.0,136.0,94.0,111.0,147.0,163.0,31.0,191.0,146.0,134.0,206.0,58.0,0.0,0.0,0.0,0.0,0.0,39.0,129.0,102.0,147.0,119.0,0.0,120.0,90.0,119.0,101.0,116.0,106.0,127.0,102.0,114.0,146.0,180.0,20.0,146.0,159.0,129.0,204.0,70.0,0.0,0.0,0.0,0.0,0.0,58.0,126.0,111.0,162.0,101.0,21.0,132.0,89.0,123.0,108.0,106.0,102.0,126.0,111.0,116.0,141.0,178.0,50.0,112.0,169.0,138.0,164.0,84.0,0.0,0.0,0.0,0.0,0.0,108.0,134.0,134.0,185.0,84.0,73.0,123.0,87.0,122.0,111.0,102.0,106.0,128.0,114.0,109.0,134.0,174.0,97.0,98.0,160.0,147.0,175.0,135.0,0.0,0.0,0.0,0.0,0.0,97.0,142.0,169.0,198.0,37.0,84.0,107.0,97.0,122.0,112.0,109.0,112.0,126.0,102.0,96.0,131.0,172.0,128.0,51.0,254.0,168.0,164.0,111.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,12.0,46.0,2.0,119.0,102.0,109.0,123.0,96.0,116.0,122.0,132.0,106.0,106.0,129.0,151.0,184.0,26.0,20.0,8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,13.0,137.0,101.0,114.0,100.0,106.0,134.0,135.0,129.0,104.0,111.0,123.0,151.0,194.0,74.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,3.0,0.0,73.0,132.0,101.0,104.0,102.0,140.0,108.0,108.0,131.0,106.0,119.0,113.0,142.0,167.0,135.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,112.0,118.0,107.0,104.0,139.0,109.0,97.0,123.0,138.0,107.0,127.0,120.0,136.0,161.0,159.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.0,0.0,142.0,148.0,112.0,105.0,101.0,83.0,125.0,123.0,143.0,104.0,115.0,100.0,126.0,168.0,178.0,7.0,0.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.0,0.0,9.0,126.0,184.0,200.0,167.0,163.0,171.0,150.0,167.0,156.0,174.0,197.0,182.0,162.0,61.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,30.0,72.0,95.0,109.0,111.0,111.0,106.0,101.0,71.0,12.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    testvector = myData[1];
    vector<DataVector> smallset = RPTreeIndex::GetInstance()->query_search(testvector, k);
    //cout << "bye" << endl;
    //smallset[1].print();
    // RPTreeIndex::GetInstance()->print();

    // for (int i = 0; i < smallset.size(); i++)
    // {
    //     cout << "Vector " << i + 1 << ": ";
    //     smallset[i].print();
    // }
    knearestneighbor(smallset, testvector, k);

    cout<<"Deleted here\n";

    //delete one vector
    RPTreeIndex::GetInstance()->DeleteData(myData[1], myData);
    testvector = myData[2];
    
    smallset = RPTreeIndex::GetInstance()->query_search(testvector, k);
    //cout << "bye" << endl;
    //smallset[1].print();
    // RPTreeIndex::GetInstance()->print();

    // for (int i = 0; i < smallset.size(); i++)
    // {
    //     cout << "Vector " << i + 1 << ": ";
    //     smallset[i].print();
    // }
    knearestneighbor(smallset, testvector, k);

    return 0;
}