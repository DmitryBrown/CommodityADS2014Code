#ifndef EMVECTOR_H
#define EMVECTOR_H
#include "File.h"
#include "../Utils/Vector.h"
#include "../Utils/Stack.h"
#include "../Utils/Sort.h"
#include "../Heaps/Heap.h"
#include "../Utils/Queue.h"
#include <cmath>
using namespace std;

namespace igmdk{

class SingleBlockBuffer
{
    long long loadedBlock;
    bool changed;
    void loadBlock(BlockFile& blockFile, long long block)
    {
        if(block != loadedBlock)
        {
            flush(blockFile);
            blockFile.readBlock(block);
            loadedBlock = block;
        }
    }
public:
    SingleBlockBuffer(): loadedBlock(-1), changed(false){}
    void flush(BlockFile& blockFile)
    {
        if(changed)
        {
            blockFile.writeBlock(loadedBlock);
            changed = false;
        }
    }
    void get(char* data, BlockFile& blockFile,
        long long block, int start, int n)
    {
        loadBlock(blockFile, block);
        for(int i = 0; i < n; ++i) data[i] = blockFile.block[start + i];
    }
    void set(char* data, BlockFile& blockFile,
        long long block, int start, int n)
    {
        loadBlock(blockFile, block);
        for(int i = 0; i < n; ++i) blockFile.block[start + i] = data[i];
        changed = true;
    }
};

template<typename POD, typename BUFFER = SingleBlockBuffer> class EMVector
{
    BlockFile blockFile;
    long long size;
    int itemsPerBlock(){return blockFile.getBlockSize()/sizeof(POD);}
    long long block(long long i){return i/itemsPerBlock();}
    long long index(long long i){return i % itemsPerBlock();}
    BUFFER buffer;
public:
    long long getSize(){return size;}
    EMVector(string const&filename, int blockSize = 2048, int extraItems = 0)
        : blockFile(filename, blockSize), size(blockFile.getSize() *
        itemsPerBlock() - extraItems){assert(blockSize % sizeof(POD) == 0);}
    long long extraItems()
        {return blockFile.getSize() * itemsPerBlock() - size;}
    void append(POD const& item)
    {
        ++size;
        if(extraItems() < 0) blockFile.appendBlock();
        set(item, size - 1);
    }
    void set(POD const& item, long long i)
    {
        assert(i >= 0 && i < size);
        char* data = (char*)&item;
        buffer.set(data, blockFile, block(i), index(i) * sizeof(POD),
            sizeof(POD));
    }
    POD operator[](long long i)
    {
        assert(i >= 0 && i < size);
        POD result;
        char* data = (char*)&result;
        buffer.get(data, blockFile, block(i), index(i) * sizeof(POD),
            sizeof(POD));
        return result;
    }
    void removeLast()
    {
        assert(size > 0);
        --size;
    }
    ~EMVector(){buffer.flush(blockFile);}

    friend void IOSort(EMVector& vector)
    {
        {
            long long C = sqrt(vector.getSize() * vector.itemsPerBlock()),
                Q = vector.getSize()/C, lastQSize = vector.getSize() % C;
            EMVector temp("IOSortTempFile.igmdk");
            typedef KVPair<POD, long long> HeapItem;
            Heap<HeapItem, KVComparator<POD, long long> > merger;
            Vector<pair<Queue<POD>, long long> > buffers(Q + 1, Q + 1);
            for(long long i = 0, k = 0; i < Q + 1; ++i)
            {
                long long n = i == Q ? lastQSize : C;
                if(n > 0)
                {//sort each block and write the result to a temp vector
                    Vector<POD> buffer;
                    for(long long j = 0; j < n; ++j)
                        buffer.append(vector[k++]);
                    quickSort(buffer.getArray(), 0, buffer.getSize() - 1);
                    for(long long j = 0; j < n; ++j) temp.append(buffer[j]);
                    //record number of unmerged items in the block
                    buffers[i].second = n - 1;
                    //put smallest item of each block on the heap
                    merger.insert(HeapItem(temp[i * n], i));
                }
            }
            for(long long i = 0; i < vector.getSize(); ++i)
            {//merge
                long long q = merger.getMin().value;
                vector.set(merger.deleteMin().key, i);
                bool bufferIsEmpty = buffers[q].first.isEmpty();
                if(!bufferIsEmpty || buffers[q].second > 0)
                {
                    if(bufferIsEmpty)
                    {//refill
                        long long j = 0, next = q * C +
                            (q == Q ? lastQSize : C) - buffers[q].second;
                        while(j < vector.itemsPerBlock() &&
                            buffers[q].second-- > 0)
                            buffers[q].first.push(temp[next + j++]);
                    }
                    merger.insert(HeapItem(buffers[q].first.pop(), q));
                }
            }
        }
        File::remove("IOSortTempFile.igmdk");
    }
};

}
#endif
