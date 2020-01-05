/**
 * File: KDTree.h
 * Author: Sudarshana Panditha
 * ------------------------
 *
 */

#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

#include "Point.h"
#include "BoundedPQueue.h"
#include <stdexcept>
#include <cmath>
#include <memory>
#include <unordered_map>

#include <iostream>

// "using namespace" in a header file is conventionally frowned upon, but I'm
// including it here so that you may use things like size_t without having to
// type std::size_t every time.
using namespace std;

template <size_t N, typename ElemType>
class KDNode;

template <size_t N, typename ElemType>
bool IsLessThanOther (const KDNode<N, ElemType> & n1,const KDNode<N, ElemType> & n2);

template <size_t N, typename ElemType>
class KDTree {
public:
    // Constructor: KDTree();
    // Usage: KDTree<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTree.
    KDTree();

    // Destructor: ~KDTree()
    // Usage: (implicit)
    // ----------------------------------------------------
    // Cleans up all resources used by the KDTree.
    ~KDTree();

    // KDTree(const KDTree& rhs);
    // KDTree& operator=(const KDTree& rhs);
    // Usage: KDTree<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Deep-copies the contents of another KDTree into this one.
    KDTree(const KDTree& other);
    KDTree& operator=(const KDTree& rhs);

    // KDTree(KDTree&& rhs);
    // KDTree& operator=(KDTree&& rhs);
    // Usage: KDTree<3, int> one = std::move(two);
    // Usage: one = std::move(two);
    // -----------------------------------------------------
    // Moves the contents of another KDTree into this one,
    // leaving the other KDTree in a valid but undefined state.
    KDTree(KDTree&& other);
    KDTree& operator=(KDTree&& rhs);

    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTree.
    size_t dimension() const;

    // size_t size() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree and whether the tree is
    // empty.
    size_t size() const;
    bool empty() const;

    // bool contains(const Point<N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTree.
    bool contains(const Point<N>& pt) const;

    // void insert(const Point<N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTree, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<N>& pt, const ElemType& value);

    // ElemType& operator[](const Point<N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTree.
    // If the point does not exist, then it is added to the KDTree using the
    // default value of ElemType as its key.
    ElemType& operator[](const Point<N>& pt);

    // ElemType& at(const Point<N>& pt);
    // const ElemType& at(const Point<N>& pt) const;
    // Usage: cout << kd.at(v) << endl;
    // ----------------------------------------------------
    // Returns a reference to the key associated with the point pt. If the point
    // is not in the tree, this function throws an out_of_range exception.
    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;

    // ElemType kNNValue(const Point<N>& key, size_t k) const
    // Usage: cout << kd.kNNValue(v, 3) << endl;
    // ----------------------------------------------------
    // Given a point v and an integer k, finds the k points in the KDTree
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<N>& key, size_t k) const;

    /* TODO: Remove this temporary function */
    void traverse() const{
      inorderTraverse(m_head);
    }

private:
    shared_ptr<KDNode<N, ElemType>> m_head; // root of the tree
    size_t m_size; // size of the tree

    // Recursive kNNValue search
    void RecursiveKNNLookup(shared_ptr<KDNode<N, ElemType>> current ,
      const Point<N>& key,  BoundedPQueue <Point<N>> & queue) const;

};

/* Node class for KDTree Structure */
template <size_t N, typename ElemType>
class KDNode {
public:
    /* Friends */
    friend class KDTree<N, ElemType>;
    friend bool IsLessThanOther <N, ElemType>(const KDNode<N, ElemType> & n1,const KDNode<N, ElemType> & n2);

    KDNode(const Point<N> &, const ElemType & ); // Constructor

    KDNode(const KDNode &);
    //KDNode(KDNode &&);
    KDNode & operator =(const KDNode &);
    //KDNode & operator =(KDNode &&);

private:
    Point<N> m_key;
    ElemType m_value;
    size_t m_split_index;
    shared_ptr<KDNode<N, ElemType>> m_left;
    shared_ptr<KDNode<N, ElemType>> m_right;
    shared_ptr<KDNode<N, ElemType>> m_parent;
};


/** KDNode Implementation **/

/* Constructor of KDNode */
template <size_t N, typename ElemType>
KDNode<N, ElemType>::KDNode(const Point<N> & key, const ElemType & value)
    :m_key{key}, m_value{value}, m_split_index{0}, m_left{nullptr}, m_right{nullptr}, m_parent{nullptr} {}

    template <size_t N, typename ElemType>
    KDNode<N, ElemType>::KDNode(const KDNode & other)
      : m_key{other.m_key}, m_value{other.m_value}, m_split_index{other.m_split_index} {
        if(other.m_left.get()) m_left = make_shared<KDNode<N, ElemType>>(KDNode<N, ElemType>{(*other.m_left)});
        if(other.m_right.get()) m_right = make_shared<KDNode<N, ElemType>>(KDNode<N, ElemType>{(*other.m_right)});
        //if(other.m_parent.get()) m_parent = make_shared<KDNode<N, ElemType>>(KDNode<N, ElemType>{(*other.m_parent)});
    }

    template <size_t N, typename ElemType>
    KDNode<N, ElemType> & KDNode<N, ElemType>::operator =(const KDNode & other) {
      m_key = other.m_key;
      m_value = other.m_value;
      m_split_index = other.m_split_index;
      *m_left = *other.m_left;
      *m_right = *other.m_right;
      //*m_parent = *other.m_parent;
    };

/*template <size_t N, typename ElemType>
KDNode<N, ElemType>::KDNode(KDNode && other)
  :m_key{other.m_key}, m_value{other.m_value}, m_split_index{other.m_split_index},
    m_left{other.m_left}, m_right{other.m_right}, m_parent{other.m_parent} {
  other.m_left = nullptr;
  other.m_right = nullptr;
  other.m_parent = nullptr;
}

template <size_t N, typename ElemType>
KDNode<N, ElemType> & KDNode<N, ElemType>::operator =(KDNode && other) {
  m_key = other.m_key;
  m_value = other.m_value;
  m_split_index = other.m_split_index;
  m_left = other.m_left;
  m_right = other.m_right;
  m_parent = other.m_parent;
  other.m_left = nullptr;
  other.m_right = nullptr;
  other.m_parent = nullptr;
};*/

/* This template function is used to compare two nodes */
template <size_t N, typename ElemType>
bool IsLessThanOther (const KDNode<N, ElemType> & n1,const KDNode<N, ElemType> & n2) {
    const size_t splitIndex = n2.m_split_index;
    auto key1 = n1.m_key;
    auto key2 = n2.m_key;
    return key1[splitIndex] < key2[splitIndex];
}

/** KDTree class implementation details */

/* KDTree Implementation */
/* Constructor */
template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() : m_size{0}, m_head{nullptr}{ // tree with no nodes

}

/* Destructor */
template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
}

/* This function returns the dimensions of the KDTree */
template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
    return N;
}

/* This function returns the number of points in the KDTree */
template <size_t N, typename ElemType >
size_t KDTree<N, ElemType>::size() const {
    return m_size;
}

/* This function returns whether the KDTree is empty or not */
template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
    return m_size == 0;
}

/*
 * This function is used to insert a new point with its
 * associated value to the KDTree
*/
template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& point, const ElemType & value){
    shared_ptr<KDNode<N, ElemType>> node = make_shared<KDNode<N, ElemType>>(point, value);
    shared_ptr<KDNode<N, ElemType>> y {nullptr};
    shared_ptr<KDNode<N, ElemType>> x = m_head;

    size_t splitIndex{0};

    // TODO: Set Split Index Correctly
    while(x.get()){
        if((*node).m_key == (*x).m_key) { // previously added node, update just the value
          auto & oldNode = (*x);
          oldNode.m_value = value;
          return;
        }
        y = x;
        if( IsLessThanOther<N,ElemType>(*node, *x)){
            x = (*x).m_left;
        } else {
            x = (*x).m_right;
        }
        ++splitIndex;
    }
    (*node).m_parent = y;
    (*node).m_split_index = (splitIndex % N);
    if(!y.get()){
        m_head = node;
    } else if (IsLessThanOther<N,ElemType>(*node,*y)){
        (*y).m_left = node;
    } else {
        (*y).m_right = node;
    }
    ++m_size;
}

/* This function returns whether the tree contains the given point or not */
template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& point) const {
  shared_ptr<KDNode<N, ElemType>> node = make_shared<KDNode<N, ElemType>>(point, ElemType{});
  shared_ptr<KDNode<N, ElemType>> x = m_head;

  while(x.get()) {
    if((*x).m_key == point) return true;
    if( IsLessThanOther<N,ElemType>(*node, *x)){
        x = (*x).m_left;
    } else {
        x = (*x).m_right;
    }
  }

  return false;
}

/* This function returns the value of a given point */
template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& point) {
  shared_ptr<KDNode<N, ElemType>> node = make_shared<KDNode<N, ElemType>>(point, ElemType{});
  shared_ptr<KDNode<N, ElemType>> y {nullptr};
  shared_ptr<KDNode<N, ElemType>> x = m_head;

  size_t splitIndex{0};

  // TODO: Set Split Index Correctly
  while(x.get()){
      if((*x).m_key == (*node).m_key) {
        return (*x).m_value;
      }
      y = x;
      if( IsLessThanOther<N,ElemType>(*node, *x)){
          x = (*x).m_left;
      } else {
          x = (*x).m_right;
      }
      ++splitIndex;
  }
  (*node).m_parent = y;
  (*node).m_split_index = (splitIndex % N);
  if(!y.get()){
      m_head = node;
  } else if (IsLessThanOther<N,ElemType>(*node,*y)){
      (*y).m_left = node;
  } else {
      (*y).m_right = node;
  }
  ++m_size;

  return (*node).m_value;
}

/*
 * This functions returns the value of the given point or throws an out_of_range
 * exception if the point doesn't exist in the tree
*/
template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& point) {
  shared_ptr<KDNode<N, ElemType>> node = make_shared<KDNode<N, ElemType>>(point, ElemType{});
  shared_ptr<KDNode<N, ElemType>> x = m_head;

  while(x.get()) {
    if((*x).m_key == point) return (*x).m_value;
    if( IsLessThanOther<N,ElemType>(*node, *x)){
        x = (*x).m_left;
    } else {
        x = (*x).m_right;
    }
  }

  throw std::out_of_range{"Point not found in the kd-tree"};
}

/*
 * This functions returns the value of the given point or throws an out_of_range
 * exception if the point doesn't exist in the tree
*/
template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& point) const {
  shared_ptr<KDNode<N, ElemType>> node = make_shared<KDNode<N, ElemType>>(point, ElemType{});
  shared_ptr<KDNode<N, ElemType>> x = m_head;

  while(x.get()) {
    if((*x).m_key == point) return (*x).m_value;
    if( IsLessThanOther<N,ElemType>(*node, *x)){
        x = (*x).m_left;
    } else {
        x = (*x).m_right;
    }
  }

  throw std::out_of_range{"Point not found in the kd-tree"};
}

/* Recursive k-NN lookup helper member function */
template <size_t N, typename ElemType>
void KDTree<N, ElemType>::RecursiveKNNLookup(shared_ptr<KDNode<N, ElemType>> current ,
  const Point<N>& key, BoundedPQueue <Point<N>>& queue) const {
    if(!current.get()) return;
    auto & currentNode = (*current);
    queue.enqueue(currentNode.m_key, Distance(currentNode.m_key, key));
    shared_ptr<KDNode<N, ElemType>> otherSubTree {nullptr};
    auto splitIndex = currentNode.m_split_index;
    if(key[splitIndex] < currentNode.m_key[splitIndex]) {
      otherSubTree = currentNode.m_right;
      RecursiveKNNLookup(currentNode.m_left, key, queue);
    } else {
      otherSubTree = currentNode.m_left;
      RecursiveKNNLookup(currentNode.m_right, key, queue);
    }

    if((queue.size() != queue.maxSize())
      || fabs((currentNode.m_key[splitIndex] - key[splitIndex]))) {
        RecursiveKNNLookup(otherSubTree, key, queue);
    }

}

/* This function implements k-nearest neighbour lookup for KDTree */
template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::kNNValue(const Point<N>& key, size_t k) const {
  BoundedPQueue <Point<N>> queue {k};
  if(k == 1 && this->contains(key)) return this->at(key);

  RecursiveKNNLookup(m_head, key, queue);

  unordered_map<ElemType, size_t> valueOccurences {};

  ElemType hfv {};
  size_t hfvFrequency {0};
  while(queue.size()){
    auto point = queue.dequeueMin();
    auto & value = this->at(point);
    if(!valueOccurences.count(value)) valueOccurences[value] = 0; // initialize value occurences
    valueOccurences[value] += 1;
    if(hfvFrequency < valueOccurences[value]) {
      hfv = value;
      hfvFrequency = valueOccurences[value];
    }
  }

  return hfv;
}

/** Copy semantics for KDTree */
template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& other) : m_head{make_shared<KDNode<N, ElemType>>(*other.m_head)}, m_size{other.m_size} {

}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs) {
  if(m_head == rhs.m_head && m_size == rhs.m_size) return *this;
  m_head = make_shared<KDNode<N, ElemType>>(*rhs.m_head);
  m_size = rhs.m_size;
  return *this;
}

/** Move semantics for KDTree  **/
template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(KDTree&& other) : m_head{other.m_head}, m_size{other.m_size} {
  other.m_head = nullptr;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(KDTree&& other) {
  m_head = other.m_head;
  m_size = other.m_size;
  other.m_head = nullptr;
}

#endif // KDTREE_INCLUDED
