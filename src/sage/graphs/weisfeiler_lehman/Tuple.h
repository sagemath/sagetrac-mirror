#ifndef EXAM_TUPLE_H
#define EXAM_TUPLE_H
#include <string>

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

	explicit Tuple(long unsigned int size);

	explicit Tuple(long int size);

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

	typedef T *iterator;

	const Tuple<T>::iterator begin() const;

	const Tuple<T>::iterator end() const;

};

#include "Tuple.tpp"
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
