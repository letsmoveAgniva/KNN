#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
using namespace std;

class DataVector {
    public:
    vector<double> v;
    // Constructor: initializes a DataVector with the given dimension (default is 0)
    DataVector(int dimension = 0);

    // Destructor: cleans up resources, called when an object goes out of scope
    ~DataVector();

    // Copy constructor: creates a new DataVector as a copy of another DataVector
    DataVector(const DataVector& other);
    // Constructor for vector<double>: creates a new DataVector from a vector<double>
    DataVector(const std::vector<double>&);

    // Copy assignment operator: assigns the values of another DataVector to this DataVector
    DataVector& operator=(const DataVector& other);

    // Assignment operator for vector<double>: assigns the values of a vector<double> to this DataVector
    DataVector& operator=(const vector<double>& vec);
    

    // Setter method to set the dimension of the DataVector
    void setDimension(int dimension = 0);

    // Addition operator overloading: performs element-wise addition of two DataVectors
    DataVector operator+(const DataVector& other);

    // Subtraction operator overloading: performs element-wise subtraction of two DataVectors
    DataVector operator-(const DataVector& other);

    // Dot product operator overloading: calculates the dot product of two DataVectors
    double operator*(const DataVector& other) const;

    // comparison operator overloading
    bool operator==(const DataVector& other) const;

    // Method to calculate the Euclidean norm of the DataVector
    double norm();

    // Method to calculate the Euclidean distance between two DataVectors
    double dist(const DataVector& other);

    // Method to print the elements of the DataVector
    void print();
};


class TreeIndex
{
protected:
    TreeIndex() {}
    ~TreeIndex() {}

public:
    static TreeIndex &GetInstance();
};


// class KDTreeIndex : public TreeIndex
// {
// public:
//     static KDTreeIndex &GetInstance();

// private:
//     KDTreeIndex() {}
//     ~KDTreeIndex() {}
// };
// class RPTreeIndex : public TreeIndex
// {
// public:
//     static RPTreeIndex &GetInstance();

// private:
//     RPTreeIndex() {}
//     ~RPTreeIndex() {}
// };