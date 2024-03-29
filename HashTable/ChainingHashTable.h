#ifndef CHAINING_HASH_TABLE_H
#define CHAINING_HASH_TABLE_H
#include "HashFunction.h"
#include "../Utils/GCFreeList.h"
namespace igmdk{

template<typename KEY, typename VALUE, typename HASHER = EHash32<BUHash>,
typename COMPARATOR = DefaultComparator<KEY> >class ChainingHashTable
{
    int capacity, size;
    struct Node
    {
        KEY key;
        VALUE value;
        Node* next;
        Node(KEY const& theKey, VALUE const& theValue):
            key(theKey), value(theValue), next(0){}
    }** table;
    Freelist<Node> freelist;
    HASHER h;
    COMPARATOR c;
    void allocateTable(int requestedSize)
    {
        int bits = lgCeiling(max(requestedSize, 8));
        capacity = twoPower(bits);
        h = HASHER(bits);
        size = 0;
        table = new Node*[capacity];
        for(int i = 0; i < capacity; ++i) table[i] = 0;
    }
    void resize()
    {
        int oldCapacity = capacity;
        Node** oldTable = table;
        allocateTable(size * 2);
        for(int i = 0; i < oldCapacity; ++i)
            for(Node* j = oldTable[i], *tail; j; j = tail)
            {
                tail = j->next;
                j->next = 0;
                insertNode(j);
            }
        delete[] oldTable;
    }
public:
    typedef Node NodeType;
    void getSize(){return size;}
    ChainingHashTable(int initialCapacity = 8, COMPARATOR const&
        theComparator = COMPARATOR()): c(theComparator)
        {allocateTable(initialCapacity);}
    ChainingHashTable(ChainingHashTable const& rhs): capacity(rhs.capacity),
        size(rhs.size), h(rhs.h), table(new Node*[capacity]), c(rhs.c)
    {
        for(int i = 0; i < capacity; ++i)
        {
            table[i] = 0;
            Node** target = &table[i];
            for(Node* j = rhs.table[i]; j; j = j->next)
            {
                *target = new(freelist.allocate())Node(*j);
                target = &(*target)->next;
            }
        }
    }
    ChainingHashTable& operator=(ChainingHashTable const& rhs)
        {return genericAssign(*this, rhs);}
    ~ChainingHashTable(){delete[] table;}
    Node* insertNode(Node* node)
    {
        Node** pointer = findPointer(node->key);
        if(*pointer)
        {
            (*pointer)->value = node->value;
            freelist.remove(node);
            node = *pointer;
        }
        else
        {
            *pointer = node;
            if(++size > capacity) resize();
        }
        return node;
    }
    Node* insert(KEY const& key, VALUE const& value)
        {return insertNode(new(freelist.allocate())Node(key, value));}
    Node** findPointer(KEY const& key)
    {
        Node** pointer = &table[h.hash(key)];
        for(;*pointer && !c.isEqual((*pointer)->key, key);
            pointer = &(*pointer)->next);
        return pointer;
    }
    Node* findNode(KEY const& key){return *findPointer(key);}
    VALUE* find(KEY const& key)
    {
        Node* next = findNode(key);
        return next ? &next->value : 0;
    }
    void remove(KEY const& key)
    {
        Node** pointer = findPointer(key);
        Node* i = *pointer;
        if(i)
        {
            *pointer = i->next;
            freelist.remove(i);
            if(--size < capacity * 0.1) resize();
        }
    }
    struct Iterator
    {
        int i;
        Node* nextLink;
        ChainingHashTable& t;
        void advance(){while(i < t.capacity && !(nextLink = t.table[i++]));}
    public:
        Iterator(ChainingHashTable& theHashTable):
            i(0), nextLink(0), t(theHashTable){advance();}
        Iterator& operator++()
        {
            if(nextLink) nextLink = nextLink->next;
            else ++i;
            advance();
            return *this;
        }
        NodeType& operator*(){assert(nextLink); return *nextLink;}
        NodeType* operator->(){assert(nextLink); return nextLink;}
        bool operator!=(Iterator const& rhs)
            {return nextLink != rhs.nextLink;}
    };
    Iterator begin(){return Iterator(*this);}
    Iterator end()
    {
        Iterator result(*this);
        result.i = capacity;
        result.nextLink = 0;
        return result;
    }
};

}
#endif
