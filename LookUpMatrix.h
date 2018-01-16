#include <unordered_map>
using namespace std;

// a = lookup[b][c]
template <typename A, typename B, typename C>
class LookUpMatrix {
    private:
        
        unordered_map<B, unordered_map<C, A>> _mat;
    public:
        bool exists(const B& b, const C& c);

        A& get(const B& b, const C& c); 

        void put(const A& a, const B& b, const C& c); 
};

template<typename A, typename B, typename C>
bool LookUpMatrix<A, B, C>::exists(const B& b, const C& c) {
    typename unordered_map<B, unordered_map<C, A>>::iterator 
        outIt = _mat.find(b);
    if(outIt != _mat.end()) {
        typename unordered_map<C, A>::iterator inIt = 
            outIt->second.find(c);
        if(inIt != outIt->second.end())
            return true;
        
    }
    return false;
}


template<typename A, typename B, typename C>
A& LookUpMatrix<A, B, C>::get(const B& b, const C& c) {
    return _mat[b][c];
}

template<typename A, typename B, typename C>
void LookUpMatrix<A, B, C>::put(const A& a, const B& b, const C& c) {
    typename unordered_map<B, unordered_map<C, A>>::iterator outIt = 
        _mat.find(b);
    if(outIt == _mat.end())
        _mat[b] = unordered_map<C, A>();
    _mat[b][c] = a;
}

