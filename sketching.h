#pragma once

#include "cs225/PNG.h"

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <set>
#include <tuple>
#include "utils.h"

using namespace cs225;

typedef uint64_t (*hashFunction)(int); 
const uint64_t GLOBAL_MAX_INT = ~(0);

std::vector<uint64_t> kminhash(std::vector<int> inList, unsigned k, hashFunction h);


std::vector<uint64_t> khash_minhash(std::vector<int> inList, std::vector<hashFunction> hv);


std::vector<uint64_t> kpartition_minhash(std::vector<int> inList, int part_bits, hashFunction h);


float minhash_jaccard(std::vector<uint64_t> mh1, std::vector<uint64_t> mh2);

int minhash_cardinality(std::vector<uint64_t> mh, unsigned k);

float exact_jaccard(std::vector<int> raw1, std::vector<int> raw2);


int exact_cardinality(std::vector<int> raw);

class MM
{
    public:
        
        MM(const cs225::PNG& input, unsigned numTiles, unsigned k, hashFunction h);

        
        std::vector<uint64_t> getMinHash(unsigned width, unsigned height) const;

        
        int countMatchTiles(const MM & other, float threshold) const;


    private:
        
        std::vector<std::vector<std::vector<uint64_t>>> min_hash_sketches_;
        std::vector<uint64_t> computeTileMinHash(unsigned int tileRow, unsigned int tileCol, unsigned int tileWidth, unsigned int tileHeight);
        void processRowTiles(unsigned int tileRow, unsigned int tileWidth, unsigned int tileHeight);
        const cs225::PNG& input_;
        unsigned numTiles_;
        unsigned k_;
        hashFunction h_;
      

};

std::vector<std::tuple<int, int, int>> build_minhash_graph(std::vector<std::string> flist, unsigned numTiles, unsigned k, hashFunction h, float threshold);
