#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <climits>
#include "LookUpMatrix.h"


using namespace std;

/*
Given an integer representing a given amount of change, write a
function to compute the total number of coins required to make
that amount of change. You can assume that there is always a
1¢ coin.

Given:
Total C to make change for.
coins 

Come up with
min number of coins to make the change
*/
class CoinChangeProb {

public:
    
    int rec(int C, vector<int>& coins) {
        if(C == 0) {
            return 0;
        } else {
            int minCoins = INT_MAX;
            for(unsigned i  = 0; i < coins.size(); ++i) {
                int temp = INT_MAX;
                if(C - coins[i] >= 0)
                    temp = 1 + rec(C - coins[i], coins);
                minCoins = min(minCoins, temp);
            }
            return minCoins;
        }        
    }

    int recMemoHelper(int C, vector<int>& coins, vector<int>& table) {
        int minCoins = INT_MAX;

        if(table[C] != INT_MAX)
            return table[C];

        if(C == 0) {
            minCoins = 0;
        } else {
            
            for(unsigned i = 0; i < coins.size(); ++i) {
                int temp = INT_MAX;
                if(C - coins[i] >= 0)
                    temp = 1 + recMemoHelper(C - coins[i], coins, table);
                // This checks for overflow in case there is no way to make
                // change.
                if(temp - 1 < temp)
                    minCoins = min(minCoins, temp);
            }
        }

        table[C] = minCoins;
        return minCoins;
    }
    int recMemo(int C, vector<int>& coins) {
        vector<int> table(C+1, INT_MAX);
        return recMemoHelper(C, coins, table);
    }

    int iter(int C, vector<int>& coins) {
        vector<int> table(C+1, INT_MAX);

        table[0] = 0;
        for(int amt = 1; amt <= C; ++amt) {
            int minCoins = INT_MAX;
            for(unsigned i = 0; i < coins.size(); ++i) { 
                int temp = INT_MAX;
                if(C - coins[i] >= 0)
                    temp = 1 + table[amt - coins[i]];
                if(temp - 1 < temp)
                    minCoins = min(minCoins, temp);

            }
            table[amt] = minCoins;
        }

        return table[C];

    }
    
};


/*
Given a 2D boolean array, find the largest square subarray of
true values. The return value should be the side length of the
largest square subarray subarray.
*/

class LargestSquare {
public:
    // For the recursive solutions, the sub problems are the 
    // subsquares whose UPPER LEFT corners are i,j.
    // FOr the iter solution, the subproblems are
    // the subsquares whose BOTTOM RIGHT corners are i, j


    
    int recForOneSq(vector<vector<bool>>& mat, int i, int j) {
        int rows = mat.size();
        int cols = mat[0].size();
        if(i >= rows || j >= cols || !mat[i][j])
            return 0;
        else {
            return 1 + min(min(recForOneSq(mat, i, j+1), 
                               recForOneSq(mat, i+1, j+1)),
                           recForOneSq(mat, i+1, j));
        }
    }

    int rec(vector<vector<bool>>& mat) {
        int rows = mat.size();
        if(rows == 0)
            return 0;

        int cols = mat[0].size(); 
        if(cols == 0)
            return 0;

        int sideLength = 0;
        for(int i = 0; i < rows; ++i) {
            for(int j = 0; j < cols; ++j) {
                if(!mat[i][j])
                    continue;
                sideLength = max(sideLength, 
                                 recForOneSq(mat, i, j));
            }
        }

        return sideLength;        
    }



    int recForOneSqMemo(vector<vector<bool>> & mat, 
                        int i, 
                        int j, 
                        vector<vector<int>>& memo)
    {
        int rows = mat.size();
        int cols = mat[0].size(); 
        int retval = 0;        
        if(i == rows || j == cols || !mat[i][j])
            retval = 0;
        else if (memo[i][j] > -1)
            retval = memo[i][j];
        else {
            retval = 1 + min(min(recForOneSq(mat, i, j+1), 
                               recForOneSq(mat, i+1, j+1)),
                           recForOneSq(mat, i+1, j));
        } 
        memo[i][j] = retval;
        return retval;
    }

    int recMemo(vector<vector<bool>>& mat) {
        vector<vector<int>> memo (mat.size()+1, vector<int>(mat[0].size()+1, -1));
        int rows = mat.size();
        if(rows == 0)
            return 0;

        int cols = mat[0].size(); 
        if(cols == 0)
            return 0;

        int sideLength = 0;
        for(int i = 0; i < rows; ++i) {
            for(int j = 0; j < cols; ++j) {
                if(!mat[i][j])
                    continue;
                sideLength = max(sideLength, 
                                 recForOneSqMemo(mat, i, j, memo));
            }
        }

        return sideLength;        
    }


    


    int iter(vector<vector<bool>> & mat) {
        vector<vector<int>> memo(mat.size(), 
                                 vector<int>(mat[0].size(), 
                                             INT_MAX));

        int rows = mat.size();
        int cols = mat[0].size(); 
        
        int sideLength = 0;
        for(int i = 0; i < rows; ++i) {
            for(int j = 0; j < cols; ++j) {
                if(i == 0 || j == 0)
                    memo[i][j] = mat[i][j] ? 1 : 0;
                else {
                    memo[i][j] = 1 + min(memo[i][j-1],
                                         min(memo[i-1][j], 
                                             memo[i-1][j-1]));
                }

                // In the recursive cases, we maximized over the 
                // largest subsquare starting at each i,j.
                // Here we maximize over all subsquares ending at i,j                
                sideLength = max(sideLength, memo[i][j]);                
            }
        }
        return sideLength;
    }
};


/*
Given a list of items with values and weights, as well as a max
weight, find the maximum value you can generate from items,
where the sum of the weights is less than or equal to the max.
eg.
items = {(w:2, v:6), (w:2, v:10), (w:3, v:12)}
max weight = 5
knapsack(items, max weight) = 22
ITEMS CAN'T BE REPEATED!!!
*/
class Knapsack01 {
public:
    struct Item {
        int weight;
        int value;
    };
    // Returns the max value of the items that will fit in 
    // the cap.

    // i keeps track of the current item.
    int recHelper(int cap, vector<Item>& items, int i) {
        if(cap <= 0)
            return 0;
        else if (i == (int)items.size())
            return 0;
        else {
            int maxValueWithItem = 0;
            int maxValueWithoutItem = 0;
            if(cap - items[i].weight >= 0) {
                maxValueWithItem = items[i].value + 
                                    recHelper(cap - items[i].weight, 
                                              items, 
                                              i+1);
            }
            maxValueWithoutItem = recHelper(cap, items, i+1);
            return max(maxValueWithItem, maxValueWithoutItem);
        }
    }

    int rec(int cap, vector<Item>& items) {
        return recHelper(cap, items, 0);
    }

    
    int recMemoHelper(int cap, 
                      vector<Item>& items, 
                      int i, 
                      unordered_map<int, vector<int>>& memo)
    {
        if(cap <= 0)
            return 0;
        else if(i == (int)items.size())
            return 0;
        else {
            unordered_map<int, vector<int>>::iterator it = memo.find(cap);
            if(it != memo.end() && it->second[i] != -1) 
                return it->second[i];
            
            int maxValueWithItem = 0;
            int maxValueWithoutItem = 0;
            if(cap - items[i].weight >= 0) {
                maxValueWithItem = items[i].value +
                                   recMemoHelper(cap - items[i].weight,
                                             items,
                                             i+1,
                                             memo);            
            }
            maxValueWithoutItem = recMemoHelper(cap, 
                                            items,
                                            i+1,
                                            memo);
            it = memo.find(cap);
            if(it == memo.end())
                memo[cap] = vector<int>(items.size(), -1);

            memo[cap][i] = max(maxValueWithItem, maxValueWithoutItem);

            return memo[cap][i];
        }
    }

    int recMemo(int cap, vector<Item>& items) {
        unordered_map<int, vector<int>> memo;
        return recMemoHelper(cap, items, 0, memo);
    }

    int iter(int cap, vector<Item>& items) {
        vector<vector<int>> memo(cap+1, vector<int>(items.size(), 0));
        for(int i = 0; i < (int)items.size(); ++i) {
            
            for(int w = 0; w <= cap; ++w) {
                int maxValueWithItem = 0;
                int maxValueWithoutItem = 0;
                if(w - items[i].weight >= 0) {
                    maxValueWithItem = max(items[i].value + memo[w - items[i].weight][i],
                                           maxValueWithItem);
                }
                if(i > 0)
                    maxValueWithoutItem = memo[w][i-1];
                memo[w][i] = max(maxValueWithItem, maxValueWithoutItem);
            }
        }
        return memo[cap][items.size()-1];
    }

    void test() {
        /*
           where W =
5, i = 0 and items={(w:2, v:6), (w:2, v:10), (w:3,
v:12)} 
        */
        vector<Item> items {{2, 6}, {2, 10}, {3, 12} };
        int cap = 5;

        int maxValueIter = iter(cap, items);
        int maxValueRecMemo = recMemo(cap, items);

        cout << "For cap = " << cap << " with items: " << endl;
        for(Item& it : items)
            cout << "weight = " << it.weight << " value = " << it.value << endl;
        cout << "maxValueIter = " << maxValueIter << " maxValueRecMemo = "  
             << maxValueRecMemo << endl;
    }
};


/*
Given an array of integers, nums and a target value T, find the
number of ways that you can add and subtract the values in
nums to add up to T.
eg.
nums = {1, 1, 1, 1, 1}
T = 3
1 + 1 + 1 + 1 - 1
1 + 1 + 1 - 1 + 1
1 + 1 - 1 + 1 + 1
1 - 1 + 1 + 1 + 1
-1 + 1 + 1 + 1 + 1
targetSum(nums, T) = 5
*/

class TargetSum {

    public:

    /*
        We can, therefore, define our subproblem as follows:
        targetSum(nums, T, i, sum) is the number of possible
        combinations of adding and subtracting the numbers at or after
        index i, where the sum of those numbers plus the sum equals
        T. It can also be stated as follows: The number of combinations
        where the sum equals T - sum.  
    */
    
    int recHelper(vector<int>& nums, int T, int i, int sum){
        if(i == (int)nums.size())
            return (sum == T) ? 1 : 0;
        else { 
            return recHelper(nums, T, i+1, sum + nums[i]) + 
                   recHelper(nums, T, i+1, sum - nums[i]);
        }
    }
    int rec(vector<int>& nums, int T){
        // Returning the number of ways to add up numbers starting at index 0
        // so that they equal T-0
        return recHelper(nums, T, 0, 0);
    }


    int recMemoHelper(vector<int>& nums, int T, int i, int sum, LookUpMatrix<int,int,int>& memo) {
        int retval = 0;
        if(i == (int)nums.size())
            retval = (sum == T) ? 1 : 0;
        else {
            if(memo.exists(sum, i))
                retval = memo.get(sum,i);
            else {                
                retval = recMemoHelper(nums, T, i+1, sum+nums[i], memo) +
                               recMemoHelper(nums, T, i+1, sum-nums[i], memo);
                memo.put(retval, sum, i);
            }
        }
        return retval;
    }

    int recMemo(vector<int>& nums, int T) {
        LookUpMatrix<int, int, int> memo;
        return recMemoHelper(nums, T, 0, 0, memo);
    }

    int iter(vector<int>& nums, int T) {
        //LookUpMatrix<int, int, int> mat;   
        int sum = 0;
        for(int n : nums)
            sum += n;
        
        // biasing sums by sum.
        // total number of possible sums is 2*sum+1 [-sum, sum]
        vector<vector<int>> cache(nums.size()+1, vector<int>(2*sum+1, 0));
        if(sum == 0) return 0;
        // Initilaze i = 0; T=0;
        cache[0][sum] = 1;
        // iterate over previous row and current
        for(int i = 1; i <= nums.size(); ++i) {
            // This loop iterates through all possible sums.
            for(int s = 0; s < 2*sum+1; ++s) {
                int prev = cache[i-1][s];
                if(prev != 0) {
                    cache[i][s - nums[i-1]] += prev;
                    cache[i][s + nums[i-1]] += prev;
                }
            }
        }



        for(unsigned k = 0; k < cache.size(); ++k) {
            for(unsigned m = 0; m < cache[0].size(); ++m) {
                cout << cache[k][m] << " ";
            }
            cout << endl;
        }

        
        return cache[nums.size()][sum+T];
    }

    void test() {
        vector<int> nums {1, 1, 1, 1, 1};
        int T = 3 ;
        
        //iter(nums, T);
        cout << recMemo(nums, T);

    }
};



int main() {

    //Knapsack01 ks;
    //ks.test();

    TargetSum ts;
    ts.test();
    
    return 0;
}




