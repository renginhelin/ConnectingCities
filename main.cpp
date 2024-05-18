/**
 * BLG 336E - Analysis of Algorithms II
 * Assignment 2
 * Rengin Helin Yalçın
 * 150220767
*/

#include <iostream>   // For input/output
#include <cmath>      // For distance calculations
#include <vector>     // For vector operations
#include <algorithm>  // For sorting operations
#include <limits>     // For min(), max() functions
#include <fstream>    // To read/write files
#include <chrono>     // For time measurement

//(These libraries are sufficient for this assignment.)

using namespace std;
using namespace std::chrono;

// Structure to represent a point in 2D space
struct Point {
    double x, y;
    // Define the == operator for comparing Point objects
    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }
};

// Helper function to calculate distance between two points
double distance(Point p1, Point p2) {
    
     return sqrt(pow((p2.x - p1.x), 2) + pow((p2.y - p1.y), 2));

}

// Helper function to compare points based on x-coordinate, can use this for sorting.
bool compareX(Point p1, Point p2) {

     return p1.x < p2.x;

}

// Helper function to compare points based on y-coordinate, can use this for sorting.
bool compareY(Point p1, Point p2) {

     return p1.y < p2.y;

}

// Function to find the closest pair of points using brute force.
pair<Point, Point> bruteForceClosestPair(vector<Point>& points, int start, int end) {

    pair<Point, Point> closest;
    double minDist = numeric_limits<double>::max();
    for (int i = start; i < end; ++i) {
        for (int j = i + 1; j < end; ++j) {
            double dist = distance(points[i], points[j]);
            if (dist < minDist) {
                minDist = dist;
                closest = {points[i], points[j]};
            }
        }
    }
    return closest;
}

// The main objective of the assignment. This is the function to find the closest pair of points recursively using divide and conquer.
pair<Point, Point> closestPair(vector<Point>& points, int start, int end) {
    if (end - start <= 3) {
        return bruteForceClosestPair(points, start, end);
    }
        // Divide the points into two halves based on the median x-coordinate
        int mid = (start + end) / 2;
        pair<Point, Point> leftPair = closestPair(points, start, mid);
        pair<Point, Point> rightPair = closestPair(points, mid, end);

        // Determine the minimum distance pair among the left and right pairs
        double leftDist = distance(leftPair.first, leftPair.second);
        double rightDist = distance(rightPair.first, rightPair.second);
        pair<Point, Point> minPair;
        double minDist;
        if (leftDist < rightDist) {
            minPair = leftPair;
            minDist = leftDist;
        } else {
            minPair = rightPair;
            minDist = rightDist;
        }

        // Merge the two sorted halves based on y-coordinate
        vector<Point> strip;
        for (int i = start; i < end; ++i) {
            if (abs(points[i].x - points[mid].x) < minDist) {
                strip.push_back(points[i]);
            }
        }
        sort(strip.begin(), strip.end(), compareY);

        // Check for closer pairs within the strip
        for (size_t i = 0; i < strip.size(); ++i) {
            for (size_t j = i + 1; j < strip.size() && (strip[j].y - strip[i].y) < minDist; ++j) {
                double dist = distance(strip[i], strip[j]);
                if (dist < minDist) {
                    minDist = dist;
                    minPair = {strip[i], strip[j]};
                }
            }
        }
        return minPair;
}


// Function to remove a pair of points from the point vector. Returns the newly formed vector.
vector<Point> removePairFromVector(vector<Point>& point_vector, pair<Point,Point> point_pair) {

    // Use the erase-remove idiom to remove elements matching point_pair from point_vector
    point_vector.erase(
        remove_if(point_vector.begin(), point_vector.end(), [&](const Point& point) {
            return (point == point_pair.first || point == point_pair.second);
        }),
        point_vector.end()
    );
    
    return point_vector;

}

// This function should wrap up the entire algorithm. It should:
//    1- Find the closest pair
//    2- Add the pair to the "pairs" vector
//    3- Remove the paired result from the map
// at each iteration.
// Then it should set "unconnected" variable as the final, unconnected city (if it exists).
void findClosestPairOrder(vector<Point> points) {
    vector<pair<Point,Point>> pairs; //add your pairs to this vector
    Point unconnected = {-1,-1}; //set this as the unconnected point if it exists, leave as is otherwise.

    auto start = high_resolution_clock::now();      // Start measuring time

    if (!is_sorted(points.begin(), points.end(), compareX)) {
        sort(points.begin(), points.end(), compareX);
    }
    // Preallocate memory for the pairs vector
    pairs.reserve(points.size() / 2);

    while (points.size() > 1) {
        // Find the closest pair of points
        pair<Point, Point> closest = closestPair(points, 0, points.size());

        // Ensure the pair is in the correct order before adding it to pairs vector
         if (closest.first.y > closest.second.y || 
            (closest.first.y == closest.second.y && closest.first.x > closest.second.x)) {
            swap(closest.first, closest.second);
        }

        // Add the closest pair to the pairs vector
        pairs.emplace_back(closest);

        // Remove the paired points from the points vector
        points = removePairFromVector(points, closest);
    }

    // If there's one unconnected point left
    if (!points.empty()) {
        unconnected = points[0];
    }

    //before printing, please make sure that within each pair, the city with the smaller y coordinate is printed first...
    //...if the y coordinates are the same, print the city with the smaller x coordinate first.

    //This part prints the output, don't modify.    
    for(size_t i = 0; i < pairs.size(); i++){
        cout << "Pair " << i+1 << ": " << pairs[i].first.x << ", " << pairs[i].first.y << " - " << pairs[i].second.x << ", " << pairs[i].second.y << endl;
    }
    if (unconnected.x != -1){
        cout << "Unconnected " << unconnected.x << ", " << unconnected.y;
    }

    auto stop = high_resolution_clock::now();                       // Stop measuring time
    auto duration = duration_cast<nanoseconds>(stop - start);       // Calculate the duration

}

//Read the coordinates from the file and convert them to a vector. Return the vector as a "vector<Point>"
vector<Point> readCoordinatesFromFile(const string& filename) {
    vector<Point> points;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return points;
    }

    Point p;
    while (file >> p.x >> p.y) {
        points.push_back(p);
    }

    file.close();

    return points;
}

// Main function. Ideally, you don't need to modify this.
int main(int argc, char* argv[]) {
    vector<Point> points = readCoordinatesFromFile(argv[1]);
    findClosestPairOrder(points);
    return 0;
}
