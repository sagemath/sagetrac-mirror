#ifndef EXAM_TUPLE_H
#define EXAM_TUPLE_H
#include <string>
#include <cstring>
#include <stdexcept>
#include <sstream>
#include <iostream>
template<typename T>
inline void hash_combine(std::size_t &seed, const T &v, bool unordered = false);
template<typename T>
class Tuple {
private:
	int numOfElements;
	T *content;
	size_t hash;
public:
	size_t getHash() const;

	void rehash();

	int size() const;

	bool operator==(const Tuple &rhs) const;

	Tuple();
        ~Tuple();
        Tuple(Tuple<T>&& b){
                //std::cerr << "Moving a tuple" << std::endl;
                this->numOfElements = b.numOfElements;
                this->content = b.content;
                this->hash = b.hash;
                b.content = nullptr;
        }
        Tuple(const Tuple<T>& b){
                //std::cerr << "Copying a tuple" << std::endl;
                this->numOfElements = b.numOfElements;
                this->content = new T[numOfElements];
                memcpy(this->content, b.content, sizeof(T)*numOfElements);
                this->hash = b.hash;
        }
        Tuple<T>& operator=(const Tuple<T>& b){
                //std::cerr << "Copying a tuple with assignment" << std::endl;
                this->numOfElements = b.numOfElements;
                this->content = new T[numOfElements];
                memcpy(this->content, b.content, sizeof(T)*numOfElements);
                this->hash = b.hash;
                return *this;
        }
        operator std::vector<int>() const{
                std::vector<int> res(numOfElements);
                for(int i = 0; i < numOfElements; i++){
                        res[i] = content[i];
                }
                return res;
        }
	explicit Tuple(long unsigned int size);

	explicit Tuple(long int size);
        Tuple(int size);
	/*
	 * Constructor that allocates the needed space for the tuple,
	 * derived from the size of the argument, and copies each
	 * element of the container given as argument inside the tuple, in order.
	 * Lastly, it also computes the hash of the Tuple at the end of the creation process.
	 */
	template<typename container>
	Tuple(const container &arr): numOfElements(arr.size()) {
		content = new T[numOfElements];
		hash = 0;
		int i = 0;
		for (const auto &elem: arr) {
			content[i] = elem;
			hash_combine(hash, content[i]);
			i++;
		}
	}
	/*
	 * Helpful constructor that does exactly the same thing the previous one does, but on 2 containers,
	 * avoiding the expense of time and space in joining them before calling the construct
	 */
	template<typename container>
	Tuple(const container &arr1, const container &arr2): numOfElements(arr1.size() + arr2.size()) {
		content = new T[numOfElements];
		hash = 0;
		int i = 0;
		for (auto &elem: arr1) {
			content[i] = elem;
			hash_combine(hash, content[i]);
			i++;
		}
		for (auto &elem: arr2) {
			content[i] = elem;
			hash_combine(hash, content[i]);
			i++;
		}
	}
	/* Constructor that allocates the needed space for the tuple,
	 * derived from the size of the argument, and copies each
	 * element of the initializer_list given as argument inside the tuple, in order.
	 * Lastly, it also computes the hash of the Tuple at the end of the creation process.
	 */
	Tuple(const std::initializer_list<T> &list) : numOfElements(list.size()) {
		content = new T[list.size()];
		int cnt = 0;
		hash = 0;
		for (auto &el: list) {
			content[cnt] = el;
			hash_combine(hash, content[cnt++]);
		}
	}

	T operator[](int idx) const;

	T &operator[](int idx);
        
        Tuple<T> modify(int idx, T value) const;
	
        template<typename U>
	friend std::ostream &operator<<(std::ostream &, const Tuple<U> &);

	bool operator<(const Tuple &other) const;

	typedef T* iterator;

	const Tuple<T>::iterator begin() const;

	const Tuple<T>::iterator end() const;
     
        __attribute__ ((noinline)) std::string to_string() const{
                if (this->size() < 1) return "";
                std::stringstream ss;
                ss << "(" << this->content[0];
                for (int i = 1; i < this->size(); i++) {
                        ss << "," << this->content[i];
                }
                ss << ")";
                return ss.str();
        }
};

/*
 * Overloaded method needed to easily print a Tuple
 */
template<typename U>
std::ostream &operator<<(std::ostream &stream, const Tuple<U> &t) {
	if (t.size() < 1) return stream;
	stream << "(" << t.content[0];
	for (int i = 1; i < t.size(); i++) {
		stream << "," << t.content[i];
	}
	stream << ")";
	return stream;
}

using namespace std;
template<typename T>
inline void hash_combine(std::size_t &seed, const T &v, bool unordered) {
        std::hash <T> hasher;
        auto vhash = hasher(v);
        if(unordered) seed ^= vhash + 0x9e3779b9 + (vhash << 6) + (vhash >> 2);
        else seed ^= vhash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/*
 * Method that returns the hash of the tuple.
 * Notice that it only returns the pre-computed hash (to save time),
 * so it's important to manually order to recompute the hash when
 * one finished modifying the tuple
 */
template<typename T>
size_t Tuple<T>::getHash() const {
	return hash;
}
template<typename T>
Tuple<T>::~Tuple(){
		if(numOfElements >= 0) delete[] content;
}
/*
 * Method to manually re-compute the hash, which will be saved for later use
 */
template<typename T>
void Tuple<T>::rehash() {
	hash = 0;
	for (int i = 0; i < numOfElements; i++) {
		hash_combine(hash, content[i]);
	}
}

template<typename T>
int Tuple<T>::size() const {
	return numOfElements;
}
/*
 * Simple equality operator overload.
 * Equality is based on equality of size first, then of the hash, then on equality of the single elements.
 * This is another reason why it's important to keep the pre-computed hash updated
 * when the tuple is in a stable state
 */
template<typename T>
bool Tuple<T>::operator==(const Tuple &rhs) const {
	if (numOfElements != rhs.numOfElements) return false;
	if(hash != rhs.hash) return false;
	for (int i = 0; i < numOfElements; i++) {
		if (content[i] != rhs.content[i]) return false;
	}
	return true;
}

/*
 * Default constructor that creates a null tuple.
 * Used, for example, to store the results of a query with no group by clause.
 * There wouldn't be any tuple under which to store the result, otherwise.
 */
template<typename T>
Tuple<T>::Tuple() {
	numOfElements = -1;
	hash = 0;
}

/*
 * Constructors that only allocate the space needed for the tuple, without assigning anything to it
 * or computing the hash
 */
template<typename T>
Tuple<T>::Tuple(long unsigned int size) : numOfElements(size) {
	content = new T[size];
}
template<typename T>
Tuple<T>::Tuple(int size) : numOfElements(size) {
	content = new T[size];
}
template<typename T>
Tuple<T>::Tuple(long int size) : numOfElements(size) {
	content = new T[size];
}

/*
 * Accessor methods that allow to read and modify the values in the tuple
 * Unsafe, it doesn't rehash
 */
template<typename T>
T Tuple<T>::operator[](int idx) const {
	return content[idx];
}

template<typename T>
T &Tuple<T>::operator[](int idx) {
	return content[idx];
}

template<typename T>
Tuple<T> Tuple<T>::modify(int idx, T value) const{
	if(idx >= numOfElements) throw out_of_range("Out of range tuple modification");
        Tuple<T> newTuple(numOfElements);
        for(int i = 0; i < numOfElements; i++){
                if(idx == i) newTuple[i] = value;
                else newTuple[i] = content[i];
        }
        newTuple.rehash();
        return newTuple;
}


/*
 * Simple overload of the less than operator.
 * A tuple is considered smaller than another if its number of elements is smaller,
 * or if of the first non-equal corresponding elements of the tuples, it has the smallest one.
 */
template<typename T>
bool Tuple<T>::operator<(const Tuple<T> &other) const {
	if (numOfElements < other.numOfElements) return true;
	if (numOfElements > other.numOfElements) return false;
	for (int i = 0; i < numOfElements; i++) {
		if (content[i] < other.content[i]) return true;
		else if (content[i] > other.content[i]) return false;
	}
	return false;
}
/*
 * .begin() and .end() methods, they simply return the content array with the appropriate offset
 */
template<typename T>
const typename Tuple<T>::iterator Tuple<T>::begin() const {
	return content;
}

template<typename T>
const typename Tuple<T>::iterator Tuple<T>::end() const {
	return content + numOfElements;
}



namespace std {
	template<typename T>
	//Hash functions for tuples, implemented in the header due to problems between the compiler and the templates
	struct hash<Tuple<T>> {
		std::size_t operator()(const Tuple<T> &f) const {
			return f.getHash();
		}
	};
}
#endif //EXAM_TUPLE_H
