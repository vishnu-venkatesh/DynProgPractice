#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <climits>
#include <iomanip>
#include <queue>
#include "LookUpMatrix.h"


/*
In the longest increasing subsequence problem, the input is a 
sequence of numbers a 1 , . . . , a n . A subsequence is any subset 
of these numbers taken in order, of the form a i 1 , a i 2 , . . . , a i k where
1 ≤ i 1 < i 2 < · · · < i k ≤ n, and an increasing subsequence is one 
in which the numbers are getting strictly larger. The task is to find 
the increasing subsequence of greatest length. For instance, 
the longest increasing subsequence of 5, 2, 8, 6, 3, 6, 9, 7 is 2, 3, 6, 9

L[j] is the length of the increasing subsequence ending at index j
L[j] = 1 + L[i] where nums[j] > nums[i], i is the end of the last inc subseq



*/
class LongestIncSubSeq {
public:
    int iter(vector<int>& nums) {
        vector<int> ssLens(nums.size(), 0);

        int i;
        int j;
        int maxLen = 0;

        if(nums.size() <= 1)
            maxLen = (int)nums.size();
        else {

            ssLens[0] = 1;
            i = 0;
            for(j = 1; j < (int)nums.size(); ++j) {
                // Find the max sub seq length thus far
                int maxLenSoFar = 0;
                for(i = 0; i < j; ++i) {
                    if(nums[j] >= nums[i] && maxLenSoFar < ssLens[i])
                        maxLenSoFar = ssLens[i];
                }
                ssLens[j] = 1 + maxLenSoFar;

                maxLen = max(ssLens[j], maxLen);
            }

            for(int ssl: ssLens)
                cout << ssl << " ";
            cout << endl;

        }
        return maxLen;
    }

    void test() {
        vector<int> nums {5, 2, 8, 6, 3, 6, 9, 7};
        cout << "LongestIncSubSeq" << endl;
        cout << "nums = ";
        for(int n : nums)
            cout << n << " ";
        cout << endl;
        cout << "length of longest inc subseq = " << iter(nums) << endl;
    }
};


/*
A natural measure of the distance between two strings is the extent to which they can be
aligned, or matched up. Technically, an alignment is simply a way of writing the strings one
above the other. For instance, here are two possible alignments of SNOWY and SUNNY:

S - N O W Y 
S U N N - Y
Cost: 3

- S N O W - Y
S U N - - N Y
Cost: 5

The “−” indicates a “gap”; any number of these can be placed in either string. The cost of an
alignment is the number of columns in which the letters differ. And the edit distance between
two strings is the cost of their best possible alignment. Do you see that there is no better
alignment of SNOWY and SUNNY than the one shown here with a cost of 3?
Edit distance is so named because it can also be thought of as the minimum number of
edits—insertions, deletions, and substitutions of characters—needed to transform the first
string into the second. For instance, the alignment shown on the left corresponds to three
edits: insert U, substitute O → N, and delete W.
*/
class EditDist {
public:
    // Sub problem is E(i,j): Edit dist of aligning x[0 .. i] and y[0 .. j]
    // Looking at the right most column when trying to align the substrings
    // results in three cases:
    //     x[i]    -      x[i]
    //      -     y[j]    y[j]
    //
    //  In case (3) If x[i] == y[j], nothing needs to be done. (cost 0)
    //  If x[i] != y[j], substitute (y[j] = x[i]) (cost 1)
    //  Case 1:
    //         cost = 1
    //         Remaining is matching x[1..i-1] with y[1..j] or E(i-1, j)
    //  Case 2: 
    //         cost = 1
    //         Remaining is matching x[1 .. i] with y[0..j-1] or E(i, j-1)
    //  Case 3:
    //         cost = 0 if x[i] = y[j], cost = 1 if x[i] != y[j]
    //         Remaining is matching x[1..i-1] with y[1..j-1] or E(i-1, j-1)
    // Min of 3 cases gives E[i,j]
    // E(i,j) = min(1+E(i-1, j), 1+E(i,j-1), diff(x[i], y[j]) + E(i-1, j-1)
    
    int calc(const string& x, const string& y) {
        vector<vector<int>> E(x.length()+1, vector<int>(y.length()+1));
        
        if(x.length() == 0)
            return y.length();
        else if(y.length() == 0)
            return x.length();   

        E[0][0] = 0;
         
        // x's substrs and the empty string
        for(unsigned i = 1; i < x.length()+1; ++i)
            E[i][0] = i;
        for(unsigned j = 1; j < y.length()+1; ++j)
            E[0][j] = j;

        for(unsigned i = 1; i <= x.length(); ++i) {
            for(unsigned j = 1; j <= y.length(); ++j) {
                int case1 = 1 + E[i-1][j];
                int case2 = 1 + E[i][j-1];
                int diff = x[i] != y[j] ? 1 : 0;
                int case3 = diff + E[i-1][j-1];

                E[i][j] = min(min(case1, case2), case3);
            }
        }

        for(vector<int>& row : E) {
            for(int elem : row)
                cout << setw(4) << elem;
            cout << endl;
        }

        return E[x.length()][y.length()];
    }

    void test() {
        string x = "polynomial";
        string y = "exponential";
        int dist = calc(x, y);

        cout << "Edit dist of " << x << " and " << y << " is " << dist << endl;
    }
};


/*
    His bag (or “knapsack”) will hold a total weight of at most W pounds. 
    There are n items to pick from, of weight w1 ,... , wn and 
    dollar value v1, ... , vn . What’s the most valuable combination of items 
    he can fit into his bag? 
*/
class Knapsack {
public:
    struct Item {
        int weight;
        int value;
    };

    // Two possible subproblems: smaller cap, fewer items.
    // First constraint suggests trying smaller cap
    // K(w) = max value acievable with cap of w
    // K(w) = max(v[i] + K(w-w[i])) for i = 1 to n

    int calcWithReps(int cap, vector<Item>& items) {
        vector<int> K(cap+1, 0);

        K[0] = 0;
        for(int w = 1; w <= cap; ++w) {
            int maxValue = 0;
            for(unsigned i = 0; i < items.size(); ++i) {
                if(w - items[i].weight >= 0)
                    maxValue 
                        = max(items[i].value + K[w- items[i].weight], maxValue);
            }
            K[w] = maxValue;
        }
        cout << " Calc table: ";
        for(int v : K)
            cout << setw(4) << v;
        cout << endl;

        return K[cap];
    }

    void testWithReps() {
        vector<Item> items {{6, 30}, {3, 14}, {4, 16}, {2, 9}};
        int cap = 10;
        cout << "Knapsack testWithReps()" << endl;
        cout << "For cap = " << cap << " and for items: " << endl;
        for(unsigned i = 0; i < items.size(); ++i)
            cout << setw(4) << i+1 << ". " 
                 << "weight = " << setw(4) << items[i].weight 
                 << " value = " << setw(4) << items[i].value << endl;
        
        int value = calcWithReps(cap, items);
        cout << "Max value = " << value << endl;

    }

    // Modify sub problems to include extra info.
    // Items 1 .. j
    // K(w, j) = max value produced by items 1 .. j for given cap w
    // K(w, j) = max(K(w-item[j].weight, j-1) + item[j].value,
    //               K(w, j-1))
    // K(0, *) = 0
    // K(*, 0) = 0; // no items
    int calcNoReps(int cap, vector<Item>& items) {
        vector<vector<int>> K(cap+1, vector<int>(items.size()+1, 0));
        for(unsigned j = 0; j <= items.size(); ++j)
            K[0][j] = 0;
        for(int w = 0; w <= cap; ++w)
            K[w][0] = 0;

        for(unsigned itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            unsigned j = itemIdx+1;

            for(int w = 1; w <= cap; ++w) {
                int pickItemVal = 0;
                int notPickItemVal  = K[w][j-1]; // not picking item.
                if(w - items[itemIdx].weight >= 0)
                    pickItemVal = 
                        K[w - items[itemIdx].weight][j-1] + items[itemIdx].value;
                K[w][j] = max(pickItemVal, notPickItemVal); 
            }

        }
        
        cout << "Knapsack no rep table: " << endl;
        for(vector<int>& row : K) {
            for(int elem : row)
                cout << setw(4) << elem;
            cout << endl;
        }
        return K[cap][items.size()];
    }

    void testCalcNoReps() {
        vector<Item> items {{6, 30}, {3, 14}, {4, 16}, {2, 9}};
        int cap = 10;
        cout << "Knapsack testCalcNoReps()" << endl;
        cout << "For cap = " << cap << " and for items: " << endl;
        for(unsigned i = 0; i < items.size(); ++i)
            cout << setw(4) << i+1 << ". " 
                 << "weight = " << setw(4) << items[i].weight 
                 << " value = " << setw(4) << items[i].value << endl;
        
        int value = calcNoReps(cap, items);
        cout << "Max value = " << value << endl;  
    }
};


/*
How do we determine the optimal order, if we want to compute 
A 1 × A 2 × · · · × A n , where the A i ’s are matrices with dimensions 
m 0 × m 1 , m 1 × m 2 , . . . , m n−1 × m n , respectively?
*/
class CostOfMatMult {
public:
    // Each ordering can be represented as a binary tree.
    // A sub problem would be a particular segment of the mutiplication
    // C(i,j) = min cost of mulitplying A[i] .. A[j]
    // We can divide the multpilication into two parts A[i] .. A[k] and
    // A[k+1] .. A[j]. Dimension of A[i] = m[i-1] x m[i] 
    // Multiplying A[i] .. A[k] gives a matrix of dimensions
    // m[i-1] x m[k]. Multiplying A[k+1] .. A[j] gives a matrix of dimensions
    // m[k] x m[j]. So combining the subproblems costs m[i-1] * m[k] * m[j]
    // C(i, j) = min(C(i, k) + C(k+1, j) + m[i-1] * m[k] *m[j]) 
    // C(i, i) = 0 // matrix itself
    int calc(vector<int>& m) {
        vector<vector<int>> C(m.size() + 1, vector<int>(m.size() + 1, 0));
        unsigned i;
        unsigned s; // subproblem size
        unsigned n = m.size();
        unsigned j;

        for(i = 1; i < m.size(); ++i)
            C[i][i] = 0;

        for(s = 1; s <= n-1; ++s) {
            for(i = 1; i <= n-s; ++i) {
                j = i + s;

                int minCost = INT_MAX;
                for(unsigned k = i; k < j; ++k) {
                    int curCost = C[i][k] + C[k+1][j] + (m[i-1] * m[k] * m[j]);
                    minCost = min(minCost, curCost);
                }
                C[i][j] = minCost;
            }
        }

        cout << "Cost table: " << endl;
        for(vector<int>& row : C) {
            for(int elem : row) 
                cout << setw(10) << elem;    
            cout << endl;
        }

        return C[1][n-1]; // THIS IS LISTED AS C[1][n] IN THE BOOK !!!
    }
    
    void test() {
        vector<int> m {50, 20, 1, 10, 100};
        
        cout << "Min cost of matrix multiplies is " << calc(m) << endl;

    }

};


struct EdgeInfo {
    string u; // Origin node
    string v; // dest node
    int dist; // Dist to that node.
    
    EdgeInfo() {}
    EdgeInfo(string orig, string dest, int d) : u(orig), v(dest), dist(d) {}
    
};

struct Node {
    string name;
    vector<EdgeInfo> edges;

    void addEdge(const EdgeInfo& ei) { edges.push_back(ei); }
};

class Graph {

    unordered_map<string, Node> namesToNodes;
public:
    void populate(vector<EdgeInfo>& edges) {
        Node tempNode;
        for(EdgeInfo& e : edges) {
            // Add the node info to the graph.
            if(namesToNodes.find(e.u) == namesToNodes.end()) {
                tempNode.name = e.u;
                namesToNodes[tempNode.name] = tempNode;
            }
            if(namesToNodes.find(e.v) == namesToNodes.end()) {
                tempNode.name = e.v;
                namesToNodes[tempNode.name] = tempNode;
            }
        }

        for(EdgeInfo& e : edges)             
            namesToNodes[e.u].addEdge(e);        
    }

    bool getNode(string name, Node& n) {
        if(namesToNodes.find(name) != namesToNodes.end()) {
            n = namesToNodes[name];
            return true;
        }
        return false;
    }

    vector<string> getAllNodeNames() {
        vector<string> retval(namesToNodes.size());
        int i = 0;
        for(unordered_map<string, Node>::iterator it = namesToNodes.begin();
            it != namesToNodes.end();
            ++it) {
            retval[i] = it->first;
            ++i;
        }

        return retval;
    }

    void print() {
        cout << setw(8) << "Node" << setw(12) << "Neighbors" << endl;
        for(unordered_map<string, Node>::iterator it = namesToNodes.begin();
            it != namesToNodes.end(); 
            ++it) {
            cout << setw(8) << it->first;
            vector<EdgeInfo>& edges = it->second.edges;
            for(EdgeInfo& e : edges) 
                cout << setw(8) << e.v << ":"<< e.dist << " ";
            
            cout << endl;
        }
    }
};



class ShortestPaths {
public:
    struct DistInfo {
        string dest;
        int dist;
        int hops;

        DistInfo() {}

        DistInfo(string theDest, int theDist, int theHops) :
            dest(theDest), dist(theDist), hops(theHops) {}

    };

    struct DistCmp {
        bool operator() (const DistInfo& a, const DistInfo& b) {
            return a.dist < b.dist;
        }
    };

    int distWithMaxHops(Graph& g, string start, string end, int maxHops) {
        int retval;
        unordered_map<string, DistInfo> distMap; // distance to each node
        priority_queue<DistInfo, std::vector<DistInfo>, DistCmp> priQ;

        Node startNode;
        Node endNode;

        if(!g.getNode(start, startNode) || !g.getNode(end, endNode))
            return INT_MAX;

        DistInfo curDist(start, 0, 0);
        distMap[start] = curDist;
        priQ.push(curDist);

        while(!priQ.empty()) {
            curDist = priQ.top();
            priQ.pop();

            if(distMap.find(curDist.dest) == distMap.end() || 
               curDist.dist < distMap[curDist.dest].dist)
                distMap[curDist.dest] = curDist;
            Node curNode; 
            g.getNode(curDist.dest, curNode);
            for(EdgeInfo& e : curNode.edges) {
                if(curDist.hops + 1 <= maxHops) {
                    priQ.push(DistInfo(e.v, 
                                       distMap[curNode.name].dist + e.dist,
                                       curDist.hops + 1));
                }                    
            }
        }

        if(distMap.find(end) != distMap.end())
            retval = distMap[end].dist;
        else
            retval = INT_MAX;

        // ====================== Debug Output ================

        priority_queue<DistInfo, std::vector<DistInfo>, DistCmp>  dbgQ;

        for(unordered_map<string, DistInfo>::iterator it = distMap.begin(); 
            it != distMap.end();
            ++it) {
            dbgQ.push(it->second);
        }

        cout << "Max hops = " << maxHops 
             << " Distances with hops to each node: " << endl;
        while(!dbgQ.empty()) {
            DistInfo d = dbgQ.top();
            dbgQ.pop();
                
            cout << d.dest << ":" << setw(4) << d.dist << " in " << setw(4)
                 << d.hops << " hops" << endl;
        }

        return retval;

    }

    typedef vector<vector<vector<int>>> DistMatType;

    DistMatType makeDistMatType(int iDim, int jDim, int kDim, int initVal) {
        return DistMatType (iDim, 
                             vector<vector<int>>(jDim, 
                                                 vector<int>(kDim,
                                                             initVal)));
    }

    // dist (i, j, 0) = direct dist from i to j (i.e. edge)
    // Let d(i, j, k) = dist from i to j using vertices 1 .. k
    // Either we can go through vertex k and get a shorter path
    // or we don't use vertex k and use only vertices 1 .. k-1
    // to go from i to j
    // d(i, j, k) = min(d(i, j, k-1), 
    //                  d(i, k, k-1) + d(k, j, k-1))                    
    void allPointsShortestPaths(Graph& g) {
        vector<string> nodeNames = g.getAllNodeNames();
        unsigned n = nodeNames.size();
        DistMatType d = makeDistMatType(n+1, n+1, n+1, 1000);
        unsigned i;
        unsigned j;
        unsigned k;
        Node nodeI;
        Node nodeJ;
        Node nodeK;
        for(i = 1; i <= n; ++i) {
            g.getNode(nodeNames[i-1], nodeI);
            for(j = 1; j <= n; ++j) {
                int edgeLen = 2000;
                for(EdgeInfo& e : nodeI.edges)
                    if(e.v == nodeNames[j-1]) {
                        edgeLen = e.dist;
                        break;
                    }
                d[i][j][0] = edgeLen;
            }
        }

        for(k = 1; k <= n; ++k)
            for(i = 1; i <= n; ++i)
                for(j = 1; j <= n; ++j)
                    d[i][j][k] = min(d[i][j][k-1],
                                     d[i][k][k-1] + d[k][j][k-1]);

        // ==================== DEBUG ONLY =====================
        cout << "Final set of distances" << endl;
        cout << setw(10) << " ";
        for(string s : nodeNames)
            cout << setw(10) << s;
        cout << endl;
        for(i = 1; i <= n; ++i) {
            cout << setw(10) << nodeNames[i-1];
            for(j = 1; j <= n; ++j)
                cout << setw(10) << d[i][j][n];
            cout << endl;
        }
                
    }


    /*
    Given the pairwise distances between cities, what is the best order
    in which to visit them, so as to minimize the overall distance traveled?
    Denote the cities by 1, . . . , n, the salesman’s hometown being 1, 
    and let D = (d ij ) be the matrix of intercity distances. 
    The goal is to design a tour that starts and ends at 1, includes
    all other cities exactly once, and has minimum total length.
    */

    vector<vector<int>> getNextSetOfSubsets(vector<vector<int>>& subsets,
                                            int n) {
        vector<vector<int>> retval; 

        for(vector<int>& cur : subsets) {
            vector<int> newss(cur.size()+1);
            unsigned k;
            for(k = 0; k < cur.size(); ++k)
                newss[k] = cur[k];

            for(k = cur[cur.size()-1]+1; k < (unsigned)n; ++k) {
                newss[cur.size()] = k;
                retval.push_back(newss);
            }
        }
        cout << "Returning subsets = ";
        for(vector<int>& subset : retval) {
            cout << "{";
            for(unsigned z = 0; z < subset.size(); ++z) {
                if(z > 0)
                    cout << ",";
                cout << subset[z];
            }
            cout << "} ";               
        }
        cout << endl;

        return retval;
    }

    unsigned long long getSubsetHash(vector<int>& subset) {
        unsigned long long retval = 0;
        for(int e : subset) 
            retval |= 1 << e;
        return retval;            
    }

    inline unsigned long long getSWithoutJ(unsigned long long S, int j) {
        return S & ~(1 << j);
    }

    enum {MOTTU_INF = 65536};

    int tspLength(vector<vector<int>>& D) {
        
        vector<vector<int>> subsets;
        int n = D.size();
        // C[S, j] = 
        LookUpMatrix<int, unsigned long long, int> C;
        
        for(int k = 1; k < n; ++k) {
            vector<int> curSS;
            curSS.push_back(0);
            curSS.push_back(k);
            subsets.push_back(curSS);
        }

        vector<vector<int>> nextSubsets;
        C.put(0, 0x1, 0);
        for(int s = 2; s <= n; ++s) {
            for(vector<int>& subset : subsets) {
                unsigned long long S = getSubsetHash(subset);
                for(int j : subset) {
                    if(j == 0) continue;

                    int minDist = MOTTU_INF;
                    unsigned long long SWoJ = getSWithoutJ(S, j); 
                    for(int i : subset) {
                        if(i == j) continue;                        
                        if(D[i][j] != MOTTU_INF && C.exists(SWoJ, i) &&
                           minDist > C.get(SWoJ,i) + D[i][j])
                            minDist = C.get(SWoJ,i) + D[i][j];
                    }
                    C.put(minDist, S, j);

                }
            }
            subsets = getNextSetOfSubsets(subsets, n);
        }
        
        // subset with 0 .. n-1
        unsigned long long S = 0;
        for(int k = 0; k < n; ++k)
            S |= 1 << k;
        // Close the loop (go from 0 ... j and then from j, back to 0)
        int retval = INT_MAX;
        for(int j = 0; j < n; ++j){
            if(C.exists(S, j) && C.get(S,j) + D[j][0] < retval)
                retval = C.get(S,j) + D[j][0];
        }
        return retval;
    }

    

    void test() {
        Graph g;

        vector<EdgeInfo> edges {{"S", "A", 1},
                                {"S", "D", 5},
                                {"S", "C", 2},
                                {"A", "S", 1},
                                {"A", "D", 5},
                                {"A", "B", 2},
                                {"C", "S", 2},
                                {"C", "D", 3},
                                {"B", "A", 2},
                                {"B", "D", 1},
                                {"B", "T", 4},
                                {"D", "B", 1},
                                {"D", "T", 1},
                                {"T", "B", 4},
                                {"T", "D", 1}};
             
        g.populate(edges);
        g.print();
        //int dist = distWithMaxHops(g, "S", "T", 1);
        //cout << "Dist = " << dist << endl;
        allPointsShortestPaths(g);
    }

    void testTsp() {
        vector<vector<int>>D {{MOTTU_INF, 2, 2, 1, 4},
                              {2, MOTTU_INF, 3, 2, 3},
                              {2, 3, MOTTU_INF, 2, 2},
                              {1, 2, 2, MOTTU_INF, 4},
                              {4, 3, 2, 4, MOTTU_INF}};
        for(vector<int>& row : D) {
            for(int elem : row) {
                cout << setw(4) << elem;
            }
            cout << endl;
        }

        cout << "Min TSP path is of length = " << tspLength(D) << endl;

    }
};



int main() {
    //LongestIncSubSeq lis;
    //lis.test();
    //EditDist ed;
    //ed.test();
    //Knapsack k;
    //k.testWithReps();
    //k.testCalcNoReps();
    //CostOfMatMult comm;
    //comm.test();
    ShortestPaths sp;
    //sp.test();
    sp.testTsp();
    return 0;
}
