#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <climits>
#include <iomanip>
#include <queue>
#include "LookUpMatrix.h"

class ContigSubseq {
public:
    /*
        Input: A list of numbers, a 1 , a 2 , . . . , a n .
        Output: The contiguous subsequence of maximum sum (a subsequence of length zero
                has sum zero). start, end, and maxSum

    */
    bool find(vector<int>& nums, unsigned& start, unsigned& end, int& maxSum) {
        if(nums.size() == 0)
            return false;
        vector<int> sums(nums.size(), 0);
        int acc = 0;
        for(unsigned z = 0; z < nums.size(); ++z) {
            acc += nums[z];
            sums[z] = acc;
        }

        maxSum = INT_MIN;
        unsigned i = 0;
        unsigned j = sums.size() - 1;
        while(i <= j) {
            int curSum = sums[j] - sums[i] + nums[i];
            if(maxSum < curSum) {          
                maxSum = curSum;
                start = i;
                end = j;
            }
            // Here we optimize over results.
            // Greedy approach using nums, results in wrong answer!!!
            if(sums[j] < sums[i])
                --j;
            else
                ++i;
        }
        return true;
    }

    void test() {
        vector<int> nums {5, 15, -30, 10, -5, 40, 10};
        unsigned start = 0, end = 0;
        int maxSum = 0;
        find(nums, start, end, maxSum);

        cout << "For nums = ";
        for(unsigned z = 0; z < nums.size(); ++z) {
            if(z > 0) cout << ", ";
            cout << nums[z];
        }
        cout << endl;
        cout << "max sum = " << maxSum << " start of seq = " << start
             << " end of seq = " << end << endl;
    }
};

// 6.2 ... This can be implemented as a greedy alg ... Dyn Prog is overkill
class LongTrip {
    
};


// 6.3
/*
Yuckdonald’s is considering opening a series of restaurants along 
Quaint Valley Highway (QVH). The n possible locations are along a straight line,
and the distances of these locations from the start of QVH are, in miles and
in increasing order, m 1 , m 2 , . . . , m n . The constraints are as follows:
- At each location, Yuckdonald’s may open at most one restaurant. 
  The expected profit from opening a restaurant at location i is p i , 
  where p i > 0 and i = 1, 2, . . . , n.
- Any two restaurants should be at least k miles apart, 
  where k is a positive integer.

MP[j] = max(MP[j-1], MP[i] + p[j] such that i is the largest index where m[j] - m[i] > k)
We can either pick this location, in which case the profit is MP[i] + p[j] 
where i meets the above condition, or we don't pick this location. In that 
case MP is the same as that of the last location j-1.
*/
class YuckDonald {
public:
    int calcMaxProfit(vector<int>& m, vector<int>& p, int k) {
        vector<int> MP(m.size(), 0);
        
        MP[0] = p[0];
        for(int j = 1; j < (int)m.size(); ++j) {
            int i;
            // We could also search for each location that comes 
            // k miles before using binary search and store these in an
            // array before hand, but this is more complex. 
            for(i = j-1; i >= 0 && (m[j] - m[i] < k); --i)
                continue;
            MP[j] = max(MP[j-1], MP[i] + p[j]);  
        }

        return MP[m.size()-1];

    }


};


/*
You are given a string of n characters s[1 . . . n], which you believe to be 
a corrupted text document in which all punctuation has vanished 
(so that it looks something like “itwasthebestoftimes...”).
You wish to reconstruct the document using a dictionary, which is available 
in the form of a Boolean function dict(·): for any string w,
    true if w is a valid word
    dict(w) = false otherwise .
(a) Give a dynamic programming algorithm that determines whether 
the string s[·] can be reconstituted as a sequence of valid words. 
The running time should be at most O(n 2 ), assuming calls to dict 
take unit time.
(b) In the event that the string is valid, make your algorithm output the corresponding se-
quence of words.
*/
class ReconstructDoc {
    unordered_set<string> dict;
public:
    void addWords(vector<string>& words) {
        for(string w: words)
            dict.insert(w);
    }

    // Assume we're looking at s[i,j]
    // cbr(0, j) = cbr(0, i) and s[i,j] is a word in dict
    // To simplify:
    // cbr[j] = cbr[i-1] and s[i,j] in dict
    // where j goes from 0 to n-1 and i goes from 0 to j
    // indexOfBegOfWord recordes the beginning of the last word that ends at
    // that index.
    bool canBeReconstructed(const string& s, vector<int>& indexOfBegOfWord) {
        vector<bool> cbr(s.length(), false);
        unsigned i;
        unsigned j;

        if(s.length() == 0)
            return true;

        cbr[0] = dict.find(s.substr(0, 1)) != dict.end();
        //cout << "Before first for loop" << endl;    
        for(j = 1; j < s.length(); ++j) {
            cbr[j] = false;
            //cout << "After initing cbr[" << j << "] to false" << endl;
            for(i = j; !cbr[j] && i > 0; --i) {
                if(cbr[i-1]) {
                    //cout << "checking substr (" << i << ", " << j << ")" << endl;
                    cbr[j] = dict.find(s.substr(i, j-i+1)) != dict.end();
                }
            }
            if(i == 0 && !cbr[j]) {
                    //cout << "check substr from beg (" << i << ", " << j << ")" << endl;
                    cbr[j] = dict.find(s.substr(i, j-i+1)) != dict.end();
                    i = -1; // the next if will set the indexOfBegWord correctly
            } 
            if(cbr[j])
                indexOfBegOfWord[j] = i+1;
        }
        return cbr[s.length()-1];
    }

    // Reconstructs the corrupted string in reverse order.
    vector<string> reconstruct(const string& s, vector<int>& indexOfBegOfWord) {
        vector<string> retval;
        for(int j = s.length()-1; j >= 0; j = indexOfBegOfWord[j]-1)
            retval.push_back(s.substr(indexOfBegOfWord[j], j-indexOfBegOfWord[j]+1));
        return retval;    
    }

    void test() {
        vector<string> words {"the", "it", "was", "best", "of", "times", 
                               "quick", "brown", "fox", "jumps", "jumped", 
                               "over",  "lazy", "dogs", "foxes", "dog",
                               "boogie", "is", "beauti"};
        vector<string> testStrs { "",
                                  "dogs",
                                  "itwasthebestoftimes",
                                  "boogieisbeautiful",
                                  "thequickbrownfoxesjumpedoverthelazydogs"
                                };
        addWords(words);
        for(string s : testStrs) {
            vector<int> indexOfBegOfWord(s.length(), -1);
            cout << "test str =  " << s << endl;
            if(canBeReconstructed(s, indexOfBegOfWord)){
                vector<string> reconInRev = reconstruct(s, indexOfBegOfWord);
                cout << "Reconstruction is: ";
                for(int i = reconInRev.size()-1; i >= 0; --i)
                    cout << reconInRev[i] << " ";
            } else {
                cout << "No Reconstruction possible";
            }
            cout << endl;            
        }                                       
    }
};

/*
    Pebbling a checkerboard. We are given a checkerboard which has 4 rows and n columns, and
has an integer written in each square. We are also given a set of 2n pebbles, and we want to
place some or all of these on the checkerboard (each pebble can be placed on exactly one square)
so as to maximize the sum of the integers in the squares that are covered by pebbles. There is
one constraint: for a placement of pebbles to be legal, no two of them can be on horizontally or
vertically adjacent squares (diagonal adjacency is fine).
(a) Determine the number of legal patterns that can occur in any column (in isolation, ignoring
the pebbles in adjacent columns) and describe these patterns.
Call two patterns compatible if they can be placed on adjacent columns to form a legal placement.
Let us consider subproblems consisting of the first k columns 1 ≤ k ≤ n. Each subproblem can
be assigned a type, which is the pattern occurring in the last column.
(b) Using the notions of compatibility and type, give an O(n)-time dynamic programming algo-
rithm for computing an optimal placement. 
*/
class PebbleCheckerBoard {
// Legal arrangements in once col (pebble denoted by x)
    // Row 0:  0    x    0    0    0    x    0    x
    // Row 1:  0    0    x    0    0    0    x    0
    // Row 2:  0    0    0    x    0    x    0    0
    // Row 3:  0    0    0    0    x    0    x    x

    
    public:

    struct Col {
        int sum;
        int patternIn;
        Col() : sum(INT_MIN), patternIn(0) {}
    };

    // Strategy is to compute the sums for all possible patterns for
    // each column.
    // Start at the final column in the board.
    // Go through each pattern.
    // For each subsequent column, pick the pattern with the max sum
    // 
    // We're starting at the end because that is the "freshest" info
    // that should determine all the other choices.
    // We're customizing the optimum for each column independent of 
    // what the optimums in the previous columns were.
    int optimalPlacement(vector<vector<int>>& board, 
                          vector<unsigned>& placements)
    {
        int n = board[0].size();
        // patterns vector: one entry for each type of pattern.
        vector<unsigned> patterns {0x0, 0x1, 0x2, 0x4, 0x8, 0x5, 0xa, 0x9};
        // check compatiblity by anding the patterns and seeing if result is 0
        
        // sums for each column. Each column contains an entry for a pattern.
        vector<vector<int>> sums(n, vector<int>(patterns.size(), 0));
        for(int j = 0; j < n; ++j) {
            for(unsigned p = 0; p < patterns.size(); ++p) {
                int curSum = 0;
                for(int i = 0; i < 4; ++i) {
                    unsigned mask = 1 << i;
                    if(mask & patterns[p])
                        curSum += board[i][j];
                }
                sums[j][p] = curSum;
            }
        }

        vector<Col> choice(n);
        
        //unsigned maxSumP = 0; // the max sum pattern in the last col.
        int maxSum = 0;
        for(unsigned p = 0; p < patterns.size(); ++p) {
            
            choice[n-1].sum = sums[n-1][p];
            choice[n-1].patternIn = p;
            int curSum = choice[n-1].sum;
            for(int j = n-2; j >= 0; --j) {
                for(unsigned q = 0; q < patterns.size(); ++q) {
                    
                    if(choice[j].sum < sums[j][q]) {
                        if(!(patterns[q] & patterns[choice[j-1].patternIn])) {
                            choice[j].sum = sums[j][q];
                            choice[j].patternIn = q; 
                        }                        
                    }
                }
                curSum += choice[j].sum;
            }
            if(maxSum < curSum) {
                maxSum = curSum;
                //maxSumP = p;
                for(int k = n-1; k >= 0; --k) 
                    placements[k] = patterns[choice[k].patternIn];
            }
        }
        return maxSum;
    }
};

// 6.6

/*
    Let us define a multiplication operation on three symbols a, b, c according to the following table;
thus ab = b, ba = c, and so on. Notice that the multiplication operation defined by the table is
neither associative nor commutative.

     a   b   c
-----------------
a |  b   b   a 
  |
b |  c   b   a  
  |
c |  a   c   c


Find an efficient algorithm that examines a string of these symbols, say bbbbac, and decides
whether or not it is possible to parenthesize the string in such a way that the value of the
resulting expression is a. For example, on input bbbbac your algorithm should return yes because
((b(bb))(ba))c = a. 
*/
//
//   R(sym, i, j) = true if s[i..j] can be parenthtisized to = sym, false otherwise.   
//   R(sym, i, j) = R(symp, i, k) * R(symq, k+1, j), where symp * symq = sym for some k.
//
//   R(sym, x, x) = x if s[x] == sym, false otherwise.
struct IndexOrBool {
    int idx;
    IndexOrBool() :  idx(-1) {}
    IndexOrBool(int x) : idx(x) {}

    operator bool() const {
        return idx > -1;
    }
    operator int() const {
        return idx;    
    }
        
};


class SpecialMult {

    
    typedef vector<IndexOrBool> BoolVec;
    typedef vector<BoolVec> BoolMat;
    //typedef unordered_map<char, char> TableColType;
    typedef unordered_map<char, unordered_map<char, char>> TableType;

    void printR(unordered_map<char, BoolMat>& R, int len) {
        for(const auto& pair : R) {
            char sym = pair.first;
            cout << "For symbol = " << sym << endl;
            for(int i = 0; i < len; ++i) {
                for(int j = 0; j < len; ++j) {
                    if(R[sym][i][j])
                        cout <<  int(R[sym][i][j]) << " ";
                    else
                        cout << "- ";
                }
                cout << endl;
            }
            cout << "----------------------" << endl;
        }
 
    }

public:
    bool canEvalTo(char c, string& s, TableType& table) {
        unordered_map<char, BoolMat> R;
        int len = s.length();
        cout << "Initializing R " << endl;
        // Initialize R.
        for(auto& pair : table) {
            R[pair.first] = BoolMat(len, BoolVec(len));
            for(int x = 0; x < len; ++x)
                R[pair.first][x][x] = IndexOrBool(s[x] == pair.first ? x : -1);
        }
        cout << " Done initilizing R " << endl;
        for(int s = 1; s < len; ++s) {
            for(const auto& pair : R) {
            char sym = pair.first;
                for(int i = 0; i < len - s; ++i) {
                    int j = i + s;
                    for(int k = i; k < j && !R[sym][i][j]; ++k) {
                        for(const auto& pairP : R) {
                            char symp = pairP.first;
                            IndexOrBool RForSymP = R[symp][i][k];
                            for(const auto& pairQ : R) {
                                char symq = pairQ.first;
                                IndexOrBool RForSymQ;
                                RForSymQ.idx = -1;
                                if(k+1 <= j)
                                    RForSymQ = R[symq][k+1][j];
                                int val = (RForSymP && RForSymQ) ? k : -1;
                                if(table[symp][symq] != sym)
                                    val = -1;
                                R[sym][i][j] = IndexOrBool(val); 
                                if(R[sym][i][j])
                                    break;
                            }
                            if(R[sym][i][j])
                                break;
                        }
                    }
                }
            }
        }
        
        printR(R, len);
        return R[c][0][len-1];
    }

    void test() {

        TableType table; 
        string syms = "abc";
        string prods = "bbacbaacc";
        cout << " Stuff " << endl;
        for(int i = 0; i < (int)syms.length(); ++i) {
            cout << " Looping over i" << endl;
            table[syms[i]] = unordered_map<char, char>();  
            for(int j = 0; j < (int)syms.length(); ++j) {
                
                cout << "Looping over j " << endl;                               
                table[syms[i]][syms[j]] = prods[i * syms.length() + j];
                cout << " prods index = " << i * syms.length() + j << endl;
            }
            cout << " end of outer loop i = " << i << endl;
        }
        cout << "Exited loop " << endl; 
        string testStr ( "bbbbac");
        cout << "For test str = " << testStr << " can evaluate to " << 'b' << "?" << endl; 
        cout << (canEvalTo('b', testStr, table) ? " Yes." : " No.") << endl;
    }
};

// 6.7
/*
A subsequence is palindromic if it is the same whether read left to right or right to left. For
instance, the sequence
A, C, G, T, G, T, C, A, A, A, A, T, C, G
has many palindromic subsequences, including A, C, G, C, A and A, A, A, A (on the other hand,
the subsequence A, C, T is not palindromic). Devise an algorithm that takes a sequence x[1 . . . n]
and returns the (length of the) longest palindromic subsequence. Its running time should be
O(n 2 ).
*/
// Longest Palindromic Sequences Ending at Each Index
// A    C    G    T    G    T    C    A    A    A    A    A    T    C    G
// 1    1    1    1    3    3    4    6    6    
//
// Whenever we see a character, what's important is when that character was seen first.
// If A[j] is the end of the palidromic subsequence, then A[i], the beginning of the subsequence,
// must == A[j]. For each character, we need to keep track of when the char was first seen,
// and also the set of locations where the character is seen.
// For any A[j], we will look up where the character was first seen. That defines the start
// of the palindromic subsequence. For each character going inward, see where it was first seen,
// past the start. Each time the condition is met, add 2 to the length. 
//

// ========================
// For each substring, find the length of the longest palindromic sub sequence
// P[i][j] = length of longest pal subseq in x[i,j]
// P[i][j] = if x[i] == x[j], 2 + P[i+1][j-1], else max(P[i+1][j-1], max(P[i][j-1], P[i+1][j])) 
//

class PalSubSeq {
public:
    int calcLen(const string& x) {
        if(x.length() < 2) return x.length();

        int len = x.length();

        vector<vector<int>> P(x.length(), vector<int>(x.length(), 0));
        for(int i = 0; i < len; ++i) 
            P[i][i] = 1;
        
        for(int i = 1; i < len; ++i) {
            if(x[i] == x[i-1])
                P[i-1][i] = 2;
            else
                P[i-1][i] = 1;
        }

        // s ends up being the size of the last sub problem solved.
        for(int s = 2; s < len; ++s) {
            for(int i = 0; i < len - s; ++i) {
                int j = i + s;
                if(x[i] == x[j]) {
                    P[i][j] = 2 + P[i+1][j-1];
                } else {
                    // Take the maximum of the three possible cases:
                    // 1 ) Neither endpoint being part of the longest palindromic subsequence
                    // 2 ) x[i] but not x[j] being part of the longest palindromic subsequence
                    // 3 ) x[j] but not x[i] being part of the longest palindromic subsequence.
                    // The text book solution only uses 2 and 3, so case 1 is covered by 2 and 3
                    P[i][j] = max(P[i+1][j-1], // This one may already be covered by the other two.
                                  max(P[i][j-1],
                                      P[i+1][j]));
                }
            }
        }
        return P[0][len-1];
    }

    void test() {
        string x = "ACGTGTCAAAATCG";
        
        cout << "Length of longest palindromic subsequence in " << x << " = " << calcLen(x) << endl;
    }
};


/*
Given two strings x = x 1 x 2 · · · x n and y = y 1 y 2 · · · y m , we wish to find 
the length of their longest common substring, that is, the largest k for which there are 
indices i and j with x i x i+1 · · · x i+k−1 = y j y j+1 · · · y j+k−1 . 
Show how to do this in time O(mn).

LCS[i][j] = length of the longest common substring ending at x[i] and y[j]
LCS[i][j] = if x[i] == y[j], 1 + LCS[i-1][j-1], 0 otherwise over x = [0, n-1] and y = [0, m-1]

*/

// Let m <= n
// x0 x1 x2 x3
// y0 y1 y2 y3 
//
// x0 x1 x2 x3
//    y0 y1 y2 y3 

class LongestCommonSubStr {
public:
    int calcLen(const string& x, const string& y) {
        int n = x.length();
        int m = y.length();
        vector<vector<int>> LCS(n, vector<int>(m, 0));
        
        int maxLen = -1;

        for(int i = 0; i < n; ++i) {
            LCS[i][0] = x[i] == y[0] ? 1 : 0;
            maxLen = max(maxLen, LCS[i][0]);
        }

        for(int j = 0; j < m; ++j) {
            LCS[0][j] = x[0] == y[j] ? 1 : 0;
            maxLen = max(maxLen, LCS[0][j]);
        }

        for(int i = 1; i < n; ++i) {
            for(int j = 1; j < m; ++j) {
                LCS[i][j] = (x[i] == y[j]) ? 1 + LCS[i-1][j-1] : 0;
                maxLen = max(maxLen, LCS[i][j]);
            }
        }

        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                cout << setw(4) << LCS[i][j];
            }
            cout << endl;
        }

        return maxLen;
    }

    void test() {
        string x = "OldSite:GeeksforGeeks.org";//"mottu loves boogie so much";
        string y = "NewSite:GeeksQuiz.com"; //"b"; //oogie mottu loves to infinity";

        cout << "x = " << x << endl;
        cout << "y = " << y << endl;
        cout << "length of longest common substring = " << calcLen(x, y) << endl;

    }
};

// 6.9
/* 
A certain string-processing language offers a primitive operation which splits a string into two
pieces. Since this operation involves copying the original string, it takes n units of time for a
string of length n, regardless of the location of the cut. Suppose, now, that you want to break a
string into many pieces. The order in which the breaks are made can affect the total running
time. For example, if you want to cut a 20-character string at positions 3 and 10, then making
the first cut at position 3 incurs a total cost of 20 + 17 = 37, while doing position 10 first has a
better cost of 20 + 10 = 30.
Give a dynamic programming algorithm that, given the locations of M cuts in a string of length
N, finds the minimum cost of breaking the string into M + 1 pieces.
*/
// 
// First cut will always be length of the string.
// We can have a table where each row is for the number of the cut (1st 2nd 3rd etc...)
// and where the column is for the cut itself. Each entry is the cost of making the cut.
//
//      3    10
// 0   10    INT_MAX
//
// 1   INT_MAX INT_MAX    
// C[i][m] = min(C[i-1][j] + substr_len for cut j) for j = 0 .. M-1
// 
// C[i][j] = length of substring + cost of cutting left substring + cost of cutting right substring.
//         = min(j - i + 1 + C[i][k-1] + C[k][j])
// cut index will always fall on the right substring.
// Picking a cut automatically makes it out of bounds in the resulting two substrings.
// So a particular cut can occur only once. So no need to mark it.
// 
constexpr int VV_INT_MAX = INT_MAX/20;
class StringSplit {
public:
        
    struct PickedCut {
        int cost;
        int k; // index of cut in string.
        PickedCut() : cost(0), k(-1) {}
        PickedCut(int argCost, int i) : cost(argCost), k(i) {};
    };

    // Use BFS to compute the order of cuts
    void computeOptOrder(int len, vector<vector<PickedCut>> & C, vector<int>& cuts,  vector<int>& orderedCuts) {
        struct Bounds {
            int i;
            int j;
            Bounds(int a, int b) : i(a), j(b) {}
        };

        queue<Bounds> q;
        int ordCutsIndex = 0;
        int numOrdCuts = orderedCuts.size();
        q.push(Bounds(0, len-1));
            
        while(!q.empty() && ordCutsIndex < numOrdCuts ) {
            Bounds curBounds = q.front();
            q.pop();

            PickedCut pc = C[curBounds.i][curBounds.j];
            cout << "curBounds.i = " << curBounds.i << " curBounds.j = " << curBounds.j 
                 << " pc.cost = " << pc.cost << " pc.k = " << pc.k << endl;
            if(pc.k != -1) {    
                Bounds left(curBounds.i, pc.k-1);
                Bounds right(pc.k, curBounds.j);
                orderedCuts[ordCutsIndex] = pc.k;
                cout << "Picked cut at " << pc.k << endl;
                ++ordCutsIndex;
                q.push(left);
                q.push(right);
            }
        }
    }

    vector<int> splitOrder(int len, vector<int>& cuts) {
        vector<int> orderedCuts(cuts.size());
        vector<vector<PickedCut>> C(len, vector<PickedCut>(len));
        
        for(int s = 1; s <= len; ++s) {
            for(int i = 0; i < len-s + 1; ++i) {
                int j = i + s - 1;
                int k = -1;
                int minCost = VV_INT_MAX;
                int indexOfPickedCut = -1;
                     
                for(unsigned z = 0; z < cuts.size(); ++z) {
                    int cost = VV_INT_MAX;

                    // We only want cuts that are inside the bounds,
                    // NOT on them.
                    if(cuts[z] <= i || cuts[z] >= j)
                    {
                        continue;
                    }
                    //cout << "cuts[z] = " << cuts[z] << endl;
                    //cout << "i = " << i << " j = " << j << endl;
 
                    k = cuts[z];
                    cost = s + C[i][k-1].cost + C[k][j].cost;
                    if(cost < minCost) {
                        minCost = cost;
                        indexOfPickedCut = z;
                    }
                }

                                
                if(indexOfPickedCut == -1)
                    continue;
                C[i][j].cost = minCost;
                C[i][j].k = cuts[indexOfPickedCut];
 
                cout << "C[" << i << "][" << j << "] = (" << minCost << ", " << indexOfPickedCut
                     << ")" << endl;

            }
        }
        cout << "Done creating matrix" << endl;
        computeOptOrder(len, C, cuts, orderedCuts);
        return orderedCuts;
    }

    void test() {
        
        int len = 20;
        vector<int> cuts {3, 10};
        vector<int> ordCuts = splitOrder(len, cuts);

        cout << "For str of len = " << len 
             << " make cuts in this order: ";
        
        for(int c : ordCuts)
            cout << c << " ";
        cout << endl;

    }
};

/*
Counting heads. Given integers n and k, along with p 1 , . . . , p n ∈ [0, 1], you want to determine the
probability of obtaining exactly k heads when n biased coins are tossed independently at random,
where p i is the probability that the ith coin comes up heads. Give an O(n 2 ) algorithm for this
task. 2 Assume you can multiply and add two numbers in [0, 1] in O(1) time.
*/
//    
//    Two ways of getting k heads with q coins:
//    k-1 coins are heads, and kth coin is a head
//    q-1 coins already have k heads and qth coin is a tail
//    P(q, k) = P(q-1, k-1)*p[k] + P(q-1, k)*(1-p[k])
//
//    Need to add up all the ways to get k heads
//    If each coin is tossed T times
//    Num heads for each coin is p1 * T, p2 * T, etc...
//    P(z tosses = k) = SUM( PRODUCT(subset with k elems
//
//    {0, 1, 2 ,3, 4, 5, 6}
//    {0, 1, 2},  
class CountingHeads {
public:
    double calcProb(int n, int k, vector<double> p) {
        vector<vector<double>> P(n, vector<double>(n, 0.0));
        
        P[0][0] = p[0];

        for(int i = 1; i < n; ++i)
            P[i][0] = p[i] + P[i-1][0]*(1-p[i]);

        for(int i = 1; i < n; ++i) {
            for(int j = 1; j < i; ++j) {
                P[i][j] = P[i-1][j-1]*p[j] + P[i-1][j]*(1-p[j]);
            }
        }

        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                cout << setw(5) << P[i][j];
            }
            cout << endl;
        }

        return P[n-1][k-1];
    }

    void test() {
        vector<double> p {0.333333, 0.5, 0.75};
        
        cout << "p = ";
        for(double z : p)
            cout << setw(5) << z;
        cout << endl;

        cout << "P(2 heads) = " << calcProb(p.size(), 2, p) << endl;
    }
};


/*
Given two strings x = x 1 x 2 · · · x n and y = y 1 y 2 · · · y m , 
we wish to find the length of their longest common subsequence, that is, 
the largest k for which there are indices i 1 < i 2 < · · · < i k and
j 1 < j 2 < · · · < j k with x i 1 x i 2 · · · x i k = y j 1 y j 2 · · · y j k . 
Show how to do this in time O(mn).
*/
//  LCS[i, j] = longest common subsequence ending at x[i] and y[j]
//  LCS[i, j] = if x[i] == y[j], 1 + LCS[i-1, j-1]
//              else max(LCS[i-1, j], max(LCS[i, j-1], LCS[i-1, j-1]))
//  LCS[0, j] = 1 if x[0] == y[j] for all j = [0, m)
//  LCS[i, 0] = 1 if y[0] == x[i] for all i = [0, n)
//
//
//
//
class LongestCommonSubSeq {
public:
    int calcLcs(const string& x, const string& y) {
        int n = x.length();
        int m = y.length();
        vector<vector<int>> LCS(x.length(), vector<int>(y.length(), 0));
        
        for(int i = 0; i < n; ++i)
            LCS[i][0] = (x[i] == y[0]) ? 1 : 0;
        for(int j = 0; j < m; ++j)
            LCS[0][j] = (x[0] == y[j]) ? 1 : 0;

        for(int i = 1; i < n; ++i) {
            for(int j = 1; j < m; ++j) {
                if(x[i] == y[j])
                    LCS[i][j] = 1 + LCS[i-1][j-1];
                else
                    LCS[i][j] = max(LCS[i-1][j],
                                    max(LCS[i][j-1], LCS[i-1][j-1]));
            }
        }
        return LCS[n-1][m-1];
    }

    void test() {
        string x = "polynomial";
        string y = "exponential";

        cout << "Length of longest common subsequence of " 
             << x << " and " << y << " is " << calcLcs(x, y) << endl;
    }
};


/*
Consider the following game. A “dealer” produces a sequence s 1 · · · s n of “cards,” face up, where
each card s i has a value v i . Then two players take turns picking a card from the sequence, but
can only pick the first or the last card of the (remaining) sequence. The goal is to collect cards of
largest total value. (For example, you can think of the cards as bills of different denominations.)
Assume n is even.

Give an O(n^2 ) algorithm to compute an optimal strategy for the first player. Given the
initial sequence, your algorithm should precompute in O(n^2 ) time some information, and
then the first player should be able to make each move optimally in O(1) time by looking
up the precomputed information.
*/
//   
//   score[i, j] = max( 
//                      min(score[i+1, j-1] + v[j], score[i+2, j] + v[i+1]) + v[i],
//                      min(score[i, j-2] + v[j-1], score[i+1, j-1] + v[i]) + v[j])
//
//
//

struct Pick {
    int scr;
    int indexChosen;
    Pick() : scr(0), indexChosen(-1) {}
    Pick(int s, int i) : scr(s), indexChosen(i) {}
};



class CardSeqGame {
    vector<vector<Pick>> score;     
    public:
    
    CardSeqGame(vector<int>& v) : score(v.size(), vector<Pick>(v.size())) {    
        int n = v.size();
        
        for(int i = 0; i < n; ++i) 
            score[i][0].scr = 0;            

        for(int j = 0; j < n; ++j)
            score[0][j].scr = 0;

        for(int z = 0; z < n; ++z) {
            score[z][z].scr = v[z];
            score[z][z].indexChosen = z;
        }

        for(int k = 0; k < n-1; ++k) {
            if(v[k] > v[k+1]) {
                score[k][k+1].scr = v[k];
                score[k][k+1].indexChosen = k;
            } else {
                score[k][k+1].scr = v[k+1];
                score[k][k+1].indexChosen = k+1;
            }
        }
        
        for(int s = 2; s < n; ++s) {
            for(int i = 0; i < n - s; ++i) {
                int j = i + s;
                
                int score_with_first = 
                    min(score[i+1][j-1].scr, score[i+2][j].scr) + v[i];
                int score_with_last = 
                    min(score[i][j-2].scr, score[i+1][j-1].scr) + v[j];

                if(score_with_first > score_with_last) {
                    score[i][j].scr = score_with_first;
                    score[i][j].indexChosen = i;
                } else {
                    score[i][j].scr = score_with_last;
                    score[i][j].indexChosen = j;
                }
            }
        }
    }

    Pick takeCard(int i, int j) {
        return score[i][j];   
    }

    void printStrategyTable() {
        cout << "Strategy Table" << endl;;
        for(int i = 0; i < score.size(); ++i) {
            for(int j = 0; j < score.size(); ++j) {
                cout << setw(4) << score[i][j].scr;
            }
            cout << endl;
        }
        cout << "----------------------------" << endl;

        for(int i = 0; i < score.size(); ++i) {
            for(int j = 0; j < score.size(); ++j) {
                cout << setw(4) << score[i][j].indexChosen;
            }
            cout << endl;
        }
    }

    static void test() {
        vector<int> v{2,10,3,20,4,1};
        CardSeqGame csg(v);

        csg.printStrategyTable();

        int i = 0;
        int j = v.size()-1;
        vector<int> player {0, 0};
        unsigned int curPlayer = 0;
        while(i < j) {
            cout << "Remaining cards (" << i << ", " << j << "): ";
            for(int k = i; k <= j; ++k)
                cout << v[k] << " ";
            cout << endl;
            Pick curPick = csg.takeCard(i, j);
            player[curPlayer] += v[curPick.indexChosen];
            if(curPick.indexChosen == i)
                ++i;
            else if (curPick.indexChosen == j)
                --j;
            else {
                ++i; --j;
                cout << "We have a problem!!!" << endl;
            }
            cout << "Player " << curPlayer << " picked card " << curPick.indexChosen 
                 << " with value " << v[curPick.indexChosen] << endl;;
            cout << "Player 0 score = " << player[0] << endl;
            cout << "Player 1 score = " << player[1] << endl;
            curPlayer = (curPlayer+1) % 2;
        }

        cout << "Player 0 Score = " << player[0] 
             << " Predicted = " << csg.takeCard(0, v.size()-1).scr << endl;
    }
};


/*
Cutting cloth. You are given a rectangular piece of cloth with dimensions X × Y , where X and
Y are positive integers, and a list of n products that can be made using the cloth. For each
product i ∈ [1, n] you know that a rectangle of cloth of dimensions a i × b i is needed and that the
final selling price of the product is c i . Assume the a i , b i , and c i are all positive 
integers. You have a machine that can cut any rectangular piece of cloth into two pieces either 
horizontally or vertically. Design an algorithm that determines the best return on the X × Y 
piece of cloth, that is, a strategy for cutting the cloth so that the products made from the 
resulting pieces give the maximum sum of selling prices. You are free to make as many copies 
of a given product as you wish, or none if desired.
*/
//
//    MR[y, x] = max ( MR[y - b, x - a] + c) 
//
//
class CuttingCloth {
public:
   struct Product {
        string name;
        int a;
        int b;
        int c;
    };

    struct Choice {
        int val;
        int index;
        int way;  // 0 if choice fits a x b, 1 if it fits b x a

        Choice() : val(0), index(-1), way(0) {}
    };

private:
    vector<vector<Choice>> MR;
    vector<Product> prods;
public:
    CuttingCloth(int X, int Y, vector<Product>& prodListArg) : 
        MR(Y+1, vector<Choice>(X+1)), prods(prodListArg) {
        int numProds = prods.size();

        // TODO: Finish population of defaults
        for(int k = 0; k < numProds; ++k) {
            Choice ch1;
            Choice ch2;
            ch2.way = 1;

            ch1.val = prods[k].c;
            ch1.index = k;
            ch2.val = prods[k].c;
            ch2.index = k;

            if(prods[k].b <= Y && prods[k].a <= X)
                MR[prods[k].b][prods[k].a] = ch1;

            if(prods[k].a <= Y && prods[k].b <= X)
                MR[prods[k].a][prods[k].b] = ch2;

        }

        for(int i = 1; i <= Y; ++i) {
            for(int j = 1; j <= X; ++j) {
                Choice optChoiceWay1;
                Choice optChoiceWay2;
                optChoiceWay2.way = 1;
                for(int k = 0; k < numProds; ++k) {

                    if((i > prods[k].b && j > prods[k].a) ||
                       (i == prods[k].b && j > prods[k].a) || 
                       (i > prods[k].b && j == prods[k].a)) {
                        int curTotPrice = MR[i - prods[k].b][j - prods[k].a].val + 
                                          MR[prods[k].b][j - prods[k].a].val + 
                                          MR[i - prods[k].b][prods[k].a].val + 
                                          prods[k].c;
                        if(curTotPrice > optChoiceWay1.val) {
                            optChoiceWay1.val = curTotPrice;
                            optChoiceWay1.index = k;
                        }
                    }

                    if((i > prods[k].a && j > prods[k].b) ||
                       (i == prods[k].a && j > prods[k].b) ||
                       (i > prods[k].a && j == prods[k].b)) {
                        int curTotPrice = MR[i - prods[k].a][j - prods[k].b].val + 
                                          MR[prods[k].a][j - prods[k].b].val + 
                                          MR[i - prods[k].a][prods[k].b].val + 
                                          prods[k].c;
                        if(curTotPrice > optChoiceWay2.val) {
                            optChoiceWay2.val = curTotPrice;
                            optChoiceWay2.index = k;
                        }
                    }
                }
                if(optChoiceWay1.val > optChoiceWay2.val)
                    MR[i][j] = optChoiceWay1;
                else
                    MR[i][j] = optChoiceWay2;
            }
        }
    }

    vector<int> maxReturn(int x, int y, int& totalSellPrice) {
        vector<int> prodList;
        totalSellPrice = 0;
        int i = y;
        int j = x;

        while ( i > 0 && j > 0) {
            totalSellPrice += MR[i][j].val;
            prodList.push_back(MR[i][j].index);

            i -= prods[MR[i][j].index].b;
            j -= prods[MR[i][j].index].a;
        }
        return prodList;
    }

    void printTable() {
        cout << "Total Returns" << endl;
        int Y = MR.size()-1;
        int X = MR[0].size()-1;
        for(int i = 0; i <= Y; ++i) {
            for(int j = 0; j <= X; ++j) {
                cout << setw(6) << MR[i][j].val;
            }
            cout << endl;
        }
        cout << "Optimal Choices" << endl;
        for(int i = 0; i <= Y; ++i) {
            for(int j = 0; j <= X; ++j) {
                cout << setw(6) << MR[i][j].index;
            }
            cout << endl;
        }
    }

    static void test() {
        vector<Product> prodList { {"A", 1, 1, 1}, {"B", 3, 1, 3}, {"C",3, 3, 9}, {"D", 2, 5, 10}};
        CuttingCloth cc(5, 6, prodList);
       
        cout << "Possible Products:" << endl;
        for(Product p: prodList) {
            cout << p.name << " = (" << p.a << ", " << p.b << ", " << p.c << ")" << endl;
        }
 
        cc.printTable();
        
        int totPrice = 0;
        vector<int> optProdChoices = cc.maxReturn(3, 4, totPrice);

        cout << "Optimal Product List (3, 4): ";
        for(int index : optProdChoices)
            cout << setw(4) << prodList[index].name;
        cout << endl;
        cout << "Total Price = " << totPrice << endl;
    }
};


class CoinChangeProblems {
public:

    /*
        Given an unlimited supply of coins of denominations x 1 , x 2 , . . . , x n , we wish to make change for
        a value v; that is, we wish to find a set of coins whose total value is v. This might not be possible:
        for instance, if the denominations are 5 and 10 then we can make change for 15 but not for 12.
        Give an O(nv) dynamic-programming algorithm for the following problem.
        Input: x 1 , . . . , x n ; v.
        Question: Is it possible to make change for v using coins of denominations x 1 , . . . , x n ?
    */
    // change is the amount of each coin needed to make change.
    // true is returned if it is possible to make change with the given coins
    bool makeChangeUnlimitedSupplyOfCoins(int val, vector<int>& x, vector<int>& change ) {
        vector<vector<int>> D (val + 1, vector<int>(x.size(), 0));
        vector<int> numCoins(val+1, -1);

        for(int k = 0; k < x.size(); ++k) {
            if(val >= x[k]) {
                numCoins[x[k]] = 1;
                D[x[k]][k] = 1;
            }
        }
        cout << "Finisihed initializing amounts for 1 coin" << endl;
        for(int i = 1; i <= val; ++i) {
            cout << "i = " << i << endl;;
            for(int k = 0; k < x.size() ; ++k) {
                cout << "   k = " << k << endl;
                if(i - x[k] > 0 && numCoins[i - x[k]] != -1){
                    numCoins[i] = numCoins[i - x[k]] + 1;
                    for(int z = 0; z < x.size(); ++z) 
                        D[i][z] = D[i - x[k]][z];
                    D[i][k] += 1;
                }
            }
        }
        cout << "Check that change was made for val" << endl;
        if(numCoins[val] != -1) {
            for(int k = 0; k < x.size(); ++k) 
                change[k] = D[val][k];
            return true;
        }
        else
            return false;
    }

    static void testMakeChangeUnlimitedSupply() {
        CoinChangeProblems ccp;
        
        vector<int> x {3, 4, 9};
        int val;

        cout << "Making change with denominations: ";
        for(int d : x)
            cout << setw(4) << d;
        cout << endl;

        cout << "Enter value to make change for or non-number to quit: ";
  
        while(cin >> val) {
            vector<int> change(x.size(), 0);
            cout << "Entered " << val << endl;;
            if(ccp.makeChangeUnlimitedSupplyOfCoins(val, x, change)) {
                cout << "Your change is: "; 
                for(int k = 0; k < change.size(); ++k)
                    if(change[k] != 0)
                        cout << setw(4) << x[k] << ":" << setw(4) << change[k];
                cout << endl;
            }
            else
                cout << "Change illah!!!" << endl;
            

            cout << "Enter value to make change for or non-number to quit: ";
        }
    }

    /*
        Consider the following variation on the change-making problem (Exercise 6.17): you are given
        denominations x 1 , x 2 , . . . , x n , and you want to make change for a value v, but you are allowed to
        use each denomination at most once. For instance, if the denominations are 1, 5, 10, 20, then you
        can make change for 16 = 1 + 15 and for 31 = 1 + 10 + 20 but not for 40 (because you can’t use 20
        twice).
        Input: Positive integers x 1 , x 2 , . . . , x n ; another integer v.
        Output: Can you make change for v, using each denomination x i at most once?
        Show how to solve this problem in time O(nv).
    */
    bool makeChangeUsingOneOfEach(int val, vector<int>& x, vector<bool>& change )     {
        // K is our knapsack :)
        // coins 1 .. j
        // With coins, the weight and value are one and the same.
        // K(w, j) = max value produced by coins 1 .. j for given cap w
        // K(w, j) = max(K[w - x[j], j-1] + x[k], K(w, j-1))
        vector<vector<int>> K(val+1, vector<int>(x.size()+1, 0));

        // Maps a value to the coin set making up that value.
        vector<vector<bool>> coinSetsForEachValue (val+1, vector<bool>(x.size(), false));

        for(unsigned j = 0; j <= x.size(); ++j) 
            K[0][j] = 0;
        for(int w = 0; w <= val; ++w)
            K[w][0] = 0;

        // Go through each coin one by one.
        for(unsigned i = 0; i < x.size(); ++i) {
            unsigned j = i + 1;

            // Try picking the coin for each value up to and including val
            for(int w = 1; w <= val; ++w) {
                int pickCoinVal = 0;
                int notPickCoinVal = K[w][j-1];
                if(w - x[i] >= 0) {
                    pickCoinVal =  K[w - x[i]][j-1] + x[i];
                    for(int z = 0; z < i; ++z)
                        coinSetsForEachValue[w][z] = coinSetsForEachValue[w - x[i]][z];
                }

                if(pickCoinVal > notPickCoinVal) { // We're picking the coin
                    coinSetsForEachValue[w][i] = true;
                    K[w][j] = pickCoinVal;
                } else {
                    coinSetsForEachValue[w][i] = false;
                    K[w][j] = notPickCoinVal;
                }    
            }
        }
        if(K[val][x.size()] == val) {
            for(int i = 0; i < change.size(); ++i)
                change[i] = coinSetsForEachValue[val][i];
            return true;
        } else
            return false;
        
    }


    static void testMakeChangeOneOfEach() {
        CoinChangeProblems ccp;
        
        vector<int> x {1, 5, 10, 20 , 25};
        int val;

        cout << "Making change with one of each denominations: ";
        for(int d : x)
            cout << setw(4) << d;
        cout << endl;

        cout << "Enter value to make change for or non-number to quit: ";
  
        while(cin >> val) {
            vector<bool> change(x.size(), false);
            cout << "Entered " << val << endl;;
            if(ccp.makeChangeUsingOneOfEach(val, x, change)) {
                cout << "Your change is: "; 
                for(int k = 0; k < change.size(); ++k)
                    if(change[k] )
                        cout << setw(4) << x[k];
                cout << endl;
            }
            else
                cout << "Change illah!!!" << endl;
            

            cout << "Enter value to make change for or non-number to quit: ";
        }
    }

    /*
     Give an efficient algorithm for the following task.
    Input: n words (in sorted order); frequencies of these words: p 1 , p 2 , . . . , p n .
    Output: The binary search tree of lowest cost (defined above as the expected number
            of comparisons in looking up a word). 
    */
    // C(i, j) = min cost of binary tree constructed from words pi .. pj
    // This bin tree can be divided into three parts consisting of words
    // from pi .. pk-1, node for pk, and subtree made from pk+1 .. pj
    // C(i, j) = min((depth + 1) C(i, k-1)
    // Cost of node = depth * frequency
    //
    //             D
    //            / \
    //           /   \
    //          B     F
    //         / \   / \
    //        /   \ /   \
    //       A    C E    G
    //
    // Cost of Subtree B A C = B + 2A + 2C
    // Cost of Subtree F E G = F + 2E + 2G
    // Cost of tree = D + 2B + 2F + 3A + 3C + 3E + 3G
    //              = D + 2B + 3A + 3C + 2F + 3E + 3G
    //              = D + (B + A + C) + B + 2A + 2C + (F + E + G) + F + 2E + 2G
    //              = D + Sum(BAC) + Cost of subtree BAC + Sum(FEG) + Cost of subtree FEG
    // C(i, j) = min(pk + Sum(pi .. pk-1) + C(i, k-1) + Sum(pk+1 .. pj) + C(k+1, j))
    // 1 1 1 1 1 1 1
    // 1 2 3 4 5 6 7
    //  
/*    struct TreeNode {
        int wordIndex;
        TreeNode* left;
        TreeNode* right;

        TreeNode() : wordIndex(-1), left(NULL), right(NULL) {} 
    };
    struct Item {
        int k;
        int cost;
        Item() : k(-1), cost(0) {}
        Item(int ak, int acost) : k(ak), cost(acost) {}
    };
    class OptimalBinarySearchTrees {
        vector<vector<Item>> C;
        vector<int> sums;
        vector<TreeNode> nodes;
    public:
        OptimalBinarySearchTrees(vector<int>& p) : C(p.size(), vector<Item>(p.size())), 
                                                   sums(p.size(), 0),
                                                   nodes(p.size()) {
            int sum = 0;
            for(int i = 0; i < p.size(); ++i) {
                sum += p[i];
                sums[i] = sum;
            }

            // Solve the problem for subproblems of size 1
            for(int i = 0; i < p.size(); ++i) 
                C[i][i] = Item(i, p[i]);
            // Solve for size 2 as well.
            for(int i = 0; i < p.size()-1; ++i) {
                if(p[i] < p[i+1])
                    C[i][i+1] = Item(i, p[i+1] + 2*p[i]);
                else
                    C[i][i+1] = Item(i+1, p[i] + 2*p[i+1]);
            }
            
            for(int s = 2; s < p.size(); ++s) {
                for(int i = 0; i < p.size(); ++i) {
                    int j = i + s;
                    Item minItem;
                    minItem.cost = MOTTU_INF;
                    for(int k = i+1; k < j; ++k) {
                        int curCost = p[k] + (sums[k-1] - sums[i] + p[i])   + C[i][k-1] 
                                           + (sums[j] - sums[k+1] + p[k+1]) + C[k+1][j];
                        if(curCost < minItem.cost) {
                            minItem.cost = curCost;
                            minItem.k = k;
                        }
                    }
                }
            }
        }
        
        struct QueueItem {
            int beg;
            int end;
            QueueItem(int b, int e) : beg(b), end(e) {}
        };
        TreeNode* makeTree(int& cost) {
            TreeNode* root;
            queue<QueueItem> nodeIndexQueue;
            int k = C[0][p.size()-1].k;

            root = nodes[C[0][p.size()-1].k].wordIndex = k;
            cost = C[0][p.size()-1].cost;
            nodeIndexQueue.push(QueueItem(0, k-1));
            nodeIndexQueue.push(QueueItem(k+1, p.size()-1));

            while(!nodeIndexQueue.empty()) {
                QueueItem item = nodeIndexQueue.
            }

            for(int z = 1; z < p.size(); ++z) {
                
            }


        }
        
    };
    */
 
};



int main() {
    //ContigSubseq cs;
    //cs.test();
    //ReconstructDoc rd;
    //rd.test();
    //SpecialMult sm;
    //sm.test();
    //PalSubSeq pss;
    //pss.test();
    //LongestCommonSubStr lcs;
    //lcs.test();
    //StringSplit ss;
    //ss.test();
    //CountingHeads ch;
    //ch.test();
    //LongestCommonSubSeq lcs;
    //lcs.test();
    //CardSeqGame::test();
    //CuttingCloth::test();
    //CoinChangeProblems::testMakeChangeUnlimitedSupply();
    CoinChangeProblems::testMakeChangeOneOfEach();
    return 0;
}
