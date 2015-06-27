#ifndef KDTREE_H
#define KDTREE_H
#include "../Utils/Utils.h"
#include "../Utils/Debug.h"
#include "../Utils/Vector.h"
#include "../Heaps/Heap.h"
#include "../HashTable/ChainingHashTable.h"
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/GCFreelist.h"
#include <cmath>
namespace igmdk{

template<typename NODE> struct QNode
{
    NODE* node;
    double d;
    bool operator<(QNode const& rhs)const{return d > rhs.d;}
    static double dh(Heap<QNode>& heap, int k)
    {
        return heap.getSize() < k ?
            numeric_limits<double>::max() : heap.getMin().d;
    }
};

template<typename KEY, typename VALUE, int D,
    typename INDEXED_COMPARATOR = LexicographicComparator<KEY> > class KDTree
{
    INDEXED_COMPARATOR c;
    struct Node
    {
        KEY key;
        VALUE value;
        Node *left, *right;
        Node(KEY const& theKey, VALUE const& theValue):
            key(theKey), value(theValue), left(0), right(0){}
    }* root;
    Freelist<Node> freelist;
public:
    typedef Node NodeType;
    bool isEmpty(){return !root;}
    KDTree(INDEXED_COMPARATOR theComparator = INDEXED_COMPARATOR()):
        root(0), c(theComparator){}
    Node* constructFrom(Node* node)
    {
        Node* tree = 0;
        if(node)
        {
            tree = new(freelist.allocate())Node(node->key, node->value);
            tree->left = constructFrom(node->left);
            tree->right = constructFrom(node->right);
        }
        return tree;
    }
    KDTree(KDTree const& other): c(other.c)
        {root = constructFrom(other.root);}
    KDTree& operator=(KDTree const& rhs){return genericAssign(*this, rhs);}

    Node** findPointer(KEY const& key, Node*& parent)
    {
        Node* node, **pointer = &root;
        parent = 0;
        for(int i = 0; (node = *pointer) && !c.isEqual(key, node->key);
            i = (i + 1) % D)
        {
            parent = node;
            pointer = &(c.isLess(key, node->key, i) ?
                node->left : node->right);
        }
        return pointer;
    }
    VALUE* find(KEY const& key)
    {
        Node *node = *findPointer(key, node);
        return node ? &node->value : 0;
    }
    void insert(KEY const& key, VALUE const& value)
    {
        Node *dummy, **pointer = findPointer(key, dummy);
        if(*pointer) (*pointer)->value = value;
        else *pointer = new(freelist.allocate())Node(key, value);
    }

    void rangeQuery(KEY const& l, KEY const& u, bool* dimensions,
        Vector<Node*>& result, Node* node, int i)
    {
        if(!node) return;
        bool inRange = true;
        for(int j = 0; j < D; ++j)
            if(dimensions[j] && (c.isLess(node->key, l, j) ||
                c.isLess(u, node->key, i))) inRange = false;
        if(inRange) result.append(node);
        int j = (i + 1) % D;
        if(!(dimensions[i] && c.isLess(node->key, l, i)))
            rangeQuery(l, u, dimensions, result, node->left, j);
        if(!(dimensions[i] && c.isLess(u, node->key, i)))
            rangeQuery(l, u, dimensions, result, node->right, j);
    }
    void rangeQuery(KEY const& l, KEY const& u, bool* dimensions,
        Vector<Node*>& result)
        {rangeQuery(l, u, dimensions, result, root, 0);}

    template<typename DISTANCE> void distanceQuery(KEY const& x,
        double radius, Vector<Node*>& result, Node* node, int i,
        DISTANCE const& distance, bool cut)
    {
        if(!node) return;
        if(distance(node->key, x) <= radius) result.append(node);
        else if(cut) return;
        else cut = true;
        i = (i + 1) % D;
        Node* nodes[] = {node->left, node->right};
        for(int j = 0; j < 2; ++j)
            distanceQuery(x, radius, result, nodes[j], i, distance, cut);
    }
    template<typename DISTANCE> void distanceQuery(KEY const& x,
        double radius, Vector<Node*>& result, DISTANCE const& distance)
        {distanceQuery(x, radius, result, root, 0, distance, false);}

    typedef QNode<Node> HEAP_ITEM;
    template<typename DISTANCE> void kNN(Node* node, KEY const& key,
        Heap<HEAP_ITEM>& heap, int k, int i, KEY& partial,
        double partialDistance, DISTANCE const& distance)
    {
        double best = HEAP_ITEM::dh(heap, k);
        if(node && partialDistance < best)
        {//update partial distance
            double newPartialDistance = distance(key, node->key, i) -
                distance(key, partial, i);
            if(heap.getSize() < k)
            {
                HEAP_ITEM x = {node, distance(key, node->key)};
                heap.insert(x);
            }
            //use new partial distance to check for a cut again
            else if(newPartialDistance < best)
            {//incremental calculate-compare
                double d = distance(best, key, node->key);
                if(d < best)
                {
                    HEAP_ITEM x = {node, d};
                    heap.changeKey(0, x);
                }
            }
            int j = (i + 1) % D;
            //swap children for best order
            Node *l = node->left, *r = node->right;
            if(!c.isLess(key, node->key, i)) swap(l, r);
            kNN(l, key, heap, k, j, partial, partialDistance, distance);
            //set partial component to the node component, use the node
            //as temporary storage
            swap(partial[i], node->key[i]);
            kNN(r, key, heap, k, j, partial, newPartialDistance, distance);
            swap(partial[i], node->key[i]);
        }
    }
    template<typename DISTANCE> Vector<NodeType*> kNN(KEY const& key, int k,
        DISTANCE const& distance)
    {
        Heap<HEAP_ITEM> heap;
        KEY partial = key;
        kNN(root, key, heap, k, 0, partial, 0, distance);
        Vector<Node*> result(k);
        while(!heap.isEmpty()) result.append(heap.deleteMin().node);
        result.reverse();
        return result;
    }
    template<typename DISTANCE> NodeType* nearestNeighbor(KEY const& key,
        DISTANCE const& distance)
    {
        assert(!isEmpty());
        Node* parent, *result = *findPointer(key, parent);
        if(result) return result;
        Heap<HEAP_ITEM> heap;
        HEAP_ITEM x = {parent, distance(key, parent->key)};
        heap.insert(x);
        KEY partial = key;
        kNN(root, key, heap, 1, 0, partial, 0, distance);
        return heap.getMin().node;
    }
};

template<typename KEY, typename VALUE, typename DISTANCE> class VpTree
{
    DISTANCE distance;
    static double bound(double keyDistance, double rLow, double rHigh)
        {return max(0., max(keyDistance - rHigh, rLow - keyDistance));}
    struct Node
    {
        KEY key;
        VALUE value;
        double leftChildDistance, radius;
        Node *left, *right;
        Node(KEY const& theKey, VALUE const& theValue): key(theKey), left(0),
            right(0), value(theValue), leftChildDistance(0), radius(0){}
        double leftChildBound(double keyDistance)
            {return bound(keyDistance, 0, leftChildDistance);}
        double rightChildBound(double keyDistance)
            {return bound(keyDistance, leftChildDistance, radius);}
    }* root;
    Freelist<Node> freelist;
public:
    typedef Node NodeType;
    bool isEmpty(){return !root;}

    VpTree(DISTANCE const& theDistance = DISTANCE()):
        root(0), distance(theDistance){}
    Node* constructFrom(Node* node)
    {
        Node* tree = 0;
        if(node)
        {
            tree = new(freelist.allocate())Node(node->key, node->value);
            tree->leftChildDistance = node->leftChildDistance;
            tree->radius = node->radius;
            tree->left = constructFrom(node->left);
            tree->right = constructFrom(node->right);
        }
        return tree;
    }
    VpTree(VpTree const& other): distance(other.distance)
        {root = constructFrom(other.root);}
    VpTree& operator=(VpTree const& rhs){return genericAssign(*this, rhs);}

    void distanceQuery(KEY const& key, double radius,
        Vector<Node*>& result, Node* node)
    {
        if(!node) return;
        double d = distance(node->key, key);
        if(d <= radius) result.append(node);
        if(node->leftChildBound(d) <= radius)
            distanceQuery(key, radius, result, node->left);
        if(node->rightChildBound(d) <= radius)
            distanceQuery(key, radius, result, node->right);
    }
    Vector<NodeType*> distanceQuery(KEY const& key, double radius)
    {
        Vector<NodeType*> result;
        distanceQuery(key, radius, result, root);
        return result;
    }

    typedef QNode<Node> HEAP_ITEM;
    void kNN(Node* node, KEY const& key, Heap<HEAP_ITEM>& heap, int k)
    {
        if(!node) return;
        //replace furthest node in heap with the current node if it's closer
        HEAP_ITEM x = {node, distance(key, node->key)};
        if(heap.getSize() < k) heap.insert(x);
        else if(x.d < HEAP_ITEM::dh(heap, k)) heap.changeKey(0, x);
        //expand closer child first
        double lb = node->leftChildBound(x.d),
            rb = node->rightChildBound(x.d);
        Node* l = node->left, *r = node->right;
        if(lb > rb)
        {
            swap(lb, rb);
            swap(l, r);
        }
        if(lb <= HEAP_ITEM::dh(heap, k)) kNN(l, key, heap, k);
        if(rb <= HEAP_ITEM::dh(heap, k)) kNN(r, key, heap, k);
    }
    Vector<NodeType*> kNN(KEY const& key, int k)
    {
        Heap<HEAP_ITEM> heap;
        kNN(root, key, heap, k);
        Vector<Node*> result(k);
        while(!heap.isEmpty()) result.append(heap.deleteMin().node);
        result.reverse();
        return result;
    }
    NodeType* nearestNeighbor(KEY const& key)
    {
        assert(!isEmpty());
        return kNN(key, 1)[0];
    }

    VALUE* find(KEY const& key)
    {
        Node* node = root;
        while(node && key != node->key)
            node = (distance(key, node->key) <= node->leftChildDistance) ?
                node->left : node->right;
        return node ? &node->value : 0;
    }
    void insert(KEY const& key, VALUE const& value)
    {
        Node **pointer = &root, *node;
        while((node = *pointer) && key != node->key)
        {
            double d = distance(key, node->key);
            node->radius = max(node->radius, d);
            if(!node->left) node->leftChildDistance = d;
            pointer = &(d <= node->leftChildDistance ?
                node->left : node->right);
        }
        if(node) node->value = value;
        else *pointer = new(freelist.allocate())Node(key, value);
    }
};

template<typename KEY, typename VALUE, typename DISTANCE> class KNNBruteForce
{
    DISTANCE distance;
    typedef KVPair<KEY, VALUE> Node;
    Vector<Node> nodes;
    struct QNode
    {
        double distance;
        int result;
        bool operator<(QNode const& rhs)const
            {return distance > rhs.distance;}
    };
public:
    KNNBruteForce(DISTANCE const& theDistance = DISTANCE()):
        distance(theDistance){}
    typedef Node NodeType;
    void insert(KEY const& key, VALUE const& value)
        {nodes.append(Node(key, value));}
    Vector<NodeType*> kNN(KEY const& key, int k)
    {
        Heap<QNode> q;
        for(int i = 0; i < nodes.getSize(); ++i)
        {
            QNode node = {distance(key, nodes[i].key), i};
            if(q.getSize() < k) q.insert(node);
            else if(node.distance < q.getMin().distance)
                q.changeKey(0, node);
        }
        Vector<NodeType*> result;
        while(!q.isEmpty()) result.append(&nodes[q.deleteMin().result]);
        result.reverse();
        return result;
    }
    NodeType* nearestNeighbor(KEY const& key){return kNN(key, 1)[0];}
};

}//end namespace
#endif
