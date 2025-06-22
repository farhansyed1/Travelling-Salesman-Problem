#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <random>

using namespace std;
using namespace chrono;

// Calculates the Euclidean distance between two points a and b
int dist(const pair<double, double> &a, const pair<double, double> &b) {
    double dx = a.first - b.first;
    double dy = a.second - b.second;
    return static_cast<int>(round(sqrt(dx * dx + dy * dy)));
}

// Calculates the distance matrix for all pairs of points
vector<vector<int>> computeDistances(const vector<pair<double, double>> &points) {
    int numNodes = points.size();
    vector<vector<int>> distanceMatrix(numNodes, vector<int>(numNodes)); // Empty matrix

    // Fill the matrix
    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < numNodes; j++) {
            distanceMatrix[i][j] = dist(points[i], points[j]);
        }
    }
    return distanceMatrix;
}

// Generates a random tour using all the nodes
vector<int> generateRandomTour(int numNodes) {
    // Create list with nodes
    vector<int> tour(numNodes); 
    for (int i = 0; i < numNodes; i++) {
        tour[i] = i;
    }
    // Shuffle the list
    random_shuffle(tour.begin(), tour.end());
    return tour;
}

// Greedy nearest neighbor heuristic - O(n^2)
vector<int> greedyAlg(int numNodes, const vector<vector<int>> &distanceMatrix){
    vector<int> tour(numNodes); 
    vector<bool> visitedNodes(numNodes, false);
    
    // Visit first node
    tour[0] = 0;
    visitedNodes[0] = true; 
    
    // Move to nearest unvisited node
    for(int i = 1; i < numNodes; i++){
        int bestNode = -1;
        // Find unvisited node (j) that is closest to current node (i) 
        for(int j = 0; j < numNodes; j++){
            // Node j unvisted + No best node found OR distance from current node to j is less than distance from current node to best node
            if(!visitedNodes[j] && (bestNode == -1 || distanceMatrix[tour[i - 1]][j] < distanceMatrix[tour[i - 1]][bestNode])) {
                bestNode = j; // if so, j becomes new best node
            }
        }
        tour[i] = bestNode; // nearest unvisited node is bestNode
        visitedNodes[bestNode] = true;
    }

    return tour;
}

// Calculate the total distance of a tour
int totalDistance(vector<int> &tour, const vector<vector<int>> &distanceMatrix) {
    int totalDist = 0;
    // Sum the distance between each pair of consecutive nodes
    for (int i = 0; i < tour.size() - 1; i++) {
        totalDist += distanceMatrix[tour[i]][tour[i + 1]];
    }
    // Add the distance from last node to first node
    totalDist += distanceMatrix[tour.back()][tour[0]];
    return totalDist;
}

// Check if Kattis time limit has been reached
bool isWithinTimeLimit(steady_clock::time_point startTime, double maxDuration){
    auto currentTime = steady_clock::now();
    double elapsedTime = duration_cast<milliseconds>(currentTime - startTime).count() / 1000.0;
    
    return elapsedTime < maxDuration;
}

// Swaps nodes i and k in the tour
void twoOptSwap(vector<int> &tour, int i, int k) {
    reverse(tour.begin() + i, tour.begin() + k + 1); // Order of nodes between i and k is reversed
}

// Two Opt optimisation that reduces total distance of a tour by swapping nodes 
vector<int> twoOpt(const vector<int> &initialTour, const vector<vector<int>> &distanceMatrix, steady_clock::time_point startTime, double maxDuration) {
    vector<int> bestTour = initialTour;
    int bestDistance = totalDistance(bestTour, distanceMatrix);
    
    bool improvement = true;
    int maxIterations = 300;
    int iteration = 0;

    // Continue swapping nodes as long as at least one improvement has been made in the nested loops
    while (improvement && iteration < maxIterations) {
        improvement = false;
        // Pick first node
        for (int i = 1; i < bestTour.size() - 1; i++) {
            // Pick second node
            for (int k = i + 1; k < bestTour.size(); k++) {
                // Stop if exceeding time
                if (!isWithinTimeLimit(startTime, maxDuration)) {
                    return bestTour;
                }
                // Perform the swap between the two nodes
                twoOptSwap(bestTour, i, k);
                // If new distance is better, save it. Otherwise, swap back. 
                int newDistance = totalDistance(bestTour, distanceMatrix);
                if (newDistance < bestDistance) {
                    bestDistance = newDistance;
                    improvement = true;
                } else {
                    twoOptSwap(bestTour, i, k);
                }
            }
        }
        iteration++;
    }

    return bestTour;
}

int main() {
    // Read number of nodes
    int numNodes;
    cin >> numNodes;

    // Read coordinates of each node
    vector<pair<double, double>> points(numNodes);
    for (int i = 0; i < numNodes; i++) {
        cin >> points[i].first >> points[i].second;
    }

    // Calculate distance matrix
    vector<vector<int>> distanceMatrix = computeDistances(points);

    // Create an initial greedy tour
    vector<int> initialTour = greedyAlg(numNodes, distanceMatrix);
    
    // Start the clock
    srand(time(0));
    auto startTime = steady_clock::now();
    double maxDuration = 1.96;

    // Optimise the initial tour with Two Opt 
    vector<int> twoOptInitialTour = twoOpt(initialTour, distanceMatrix, startTime, maxDuration);

    vector<int> bestTour = twoOptInitialTour;
    int bestDist = totalDistance(twoOptInitialTour, distanceMatrix);

    // While there is still time, generate a random tour, optimise it, and see if it better than the best found tour. 
    while (isWithinTimeLimit(startTime, maxDuration)) {
        vector<int> randomTour = generateRandomTour(numNodes);
        vector<int> twoOptRandomTour = twoOpt(randomTour, distanceMatrix, startTime, maxDuration);

        int dist = totalDistance(twoOptRandomTour, distanceMatrix);

        if (dist < bestDist) {
            bestTour = twoOptRandomTour;
            bestDist = dist;
        }
    }

    // Output the results
    for (int node : bestTour) {
        cout << node << '\n';
    }

    return 0;
}