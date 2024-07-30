#include "sketching.h"
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <cmath>

std::vector<uint64_t> kminhash(std::vector<int> inList, unsigned k, hashFunction h) {
    // code
    std::set<uint64_t> unique_hash_values;
    for (int number : inList) {
        unique_hash_values.insert(h(number));
    }
    std::vector<uint64_t> sorted_hash_values(unique_hash_values.begin(), unique_hash_values.end());
    std::sort(sorted_hash_values.begin(), sorted_hash_values.end());

    std::vector<uint64_t> result;
    result.reserve(k);

    for (size_t i = 0; i < k && i < sorted_hash_values.size(); ++i) {
        result.push_back(sorted_hash_values[i]);
    }
    while (result.size() < k) {
        result.push_back(UINT64_MAX);
    }
    return result;
}

std::vector<uint64_t> khash_minhash(std::vector<int> inList, std::vector<hashFunction> hv) {
    // code
    std::vector<uint64_t> minHashes(hv.size(), UINT64_MAX);
    for (size_t i = 0; i < hv.size(); ++i) {
        for (int num : inList) {
            uint64_t hash_value = hv[i](num);
            minHashes[i] = std::min(minHashes[i], hash_value);
        }
    }
    return minHashes;
}

std::vector<uint64_t> kpartition_minhash(std::vector<int> inList, int part_bits, hashFunction h) {
    // code
    size_t num_partitions = 1ULL << part_bits;
    std::vector<uint64_t> minHashes(num_partitions, UINT64_MAX);
    for (int num : inList) {
        uint64_t hash_value = h(num);
        size_t partition_index = hash_value >> (64 - part_bits);
        minHashes[partition_index] = std::min(minHashes[partition_index], hash_value);
    }
    return minHashes;
}

float minhash_jaccard(std::vector<uint64_t> mh1, std::vector<uint64_t> mh2) {
    // code
    const uint64_t GLOBAL_MAX_INT = UINT64_MAX;
    size_t intersection_count = 0;
    std::unordered_set<uint64_t> union_set;
    for (size_t i = 0; i < mh1.size(); ++i) {
        if (mh1[i] != GLOBAL_MAX_INT) {
            union_set.insert(mh1[i]);
            if (std::find(mh2.begin(), mh2.end(), mh1[i])!=mh2.end()) {
                intersection_count++;
            }
        }
    }
    for (size_t i = 0; i < mh2.size(); ++i) {
        if (mh2[i] != GLOBAL_MAX_INT) {
            union_set.insert(mh2[i]);
        }
    }
    size_t union_count = union_set.size();
    if (union_count == 0) {
        return 0.0f;
    }
    return static_cast<float>(intersection_count) / static_cast<float>(union_count);
}

int minhash_cardinality(std::vector<uint64_t> mh, unsigned k) {
    // code
    if (k > mh.size() || mh.empty()) {
        return 0;
    }
    float normalized_kth_min = static_cast<float>(mh[k-1])/static_cast<float>(GLOBAL_MAX_INT);
    float estimated_card = (k/normalized_kth_min)-1;
    int rounded_up_card = static_cast<int>(std::ceil(estimated_card));
    return rounded_up_card;
}

float exact_jaccard(std::vector<int> raw1, std::vector<int> raw2) {
    // code
    std::set<int> set1(raw1.begin(), raw1.end());
    std::set<int> set2(raw2.begin(), raw2.end());
    std::vector<int> intersection;
    std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), std::back_inserter(intersection));
    std::vector<int> union_set;
    std::set_union(set1.begin(), set1.end(), set2.begin(), set2.end(), std::back_inserter(union_set));
    return union_set.empty() ? 0.0f : static_cast<float>(intersection.size())/static_cast<float>(union_set.size());
}

int exact_cardinality(std::vector<int> raw) {
    // code
    std::set<int> unique_elements(raw.begin(), raw.end());
    return unique_elements.size();
}

////////////////// PART 2 ////////////////////
MM::MM(const cs225::PNG& input, unsigned numTiles, unsigned k, hashFunction h): input_(input), numTiles_(numTiles), k_(k), h_(h) {
    
    unsigned int tileWidth = std::ceil(static_cast<float>(input.width()) / numTiles);
    unsigned int tileHeight = std::ceil(static_cast<float>(input.height()) / numTiles);

    min_hash_sketches_.resize(numTiles, std::vector<std::vector<uint64_t>>(numTiles));

    for (unsigned int tileRow = 0; tileRow < numTiles; ++tileRow) {
        processRowTiles(tileRow, tileWidth, tileHeight);
    }
}

void MM::processRowTiles(unsigned int tileRow, unsigned int tileWidth, unsigned int tileHeight) {
    for (unsigned int tileCol = 0; tileCol < numTiles_; ++tileCol) {
        min_hash_sketches_[tileRow][tileCol] = computeTileMinHash(tileRow, tileCol, tileWidth, tileHeight);
    }
}

std::vector<uint64_t> MM::computeTileMinHash(unsigned int tileRow, unsigned int tileCol, unsigned int tileWidth, unsigned int tileHeight) {
    std::vector<int> lvals;
    for (unsigned int x = tileRow * tileWidth; x < std::min((tileRow + 1) * tileWidth, input_.width()); ++x) {
        for (unsigned int y = tileCol * tileHeight; y < std::min((tileCol + 1) * tileHeight, input_.height()); ++y) {
            lvals.push_back(static_cast<int>(input_.getPixel(x, y).l * 255));
        }
    }
    return kminhash(lvals, k_, h_);
}


std::vector<uint64_t> MM::getMinHash(unsigned width, unsigned height) const {
    return min_hash_sketches_.at(width).at(height);
}

int MM::countMatchTiles(const MM & other, float threshold) const {
    int matchingTilesCount = 0;
    for (unsigned i = 0; i < numTiles_; ++i) {
        for (unsigned j = 0; j < numTiles_; ++j) {
            std::vector<uint64_t> thisTileMinHash = min_hash_sketches_[i][j];
            std::vector<uint64_t> otherTileMinHash = other.min_hash_sketches_[i][j];
            float jaccardSimilarity = minhash_jaccard(thisTileMinHash, otherTileMinHash);
            if (jaccardSimilarity >= threshold) {
                matchingTilesCount++;
            }
        }
    }
    return matchingTilesCount;
}

std::vector<std::tuple<int, int, int>> build_minhash_graph(std::vector<std::string> flist, unsigned numTiles, unsigned k, hashFunction h, float threshold) {
    std::vector<MM> minHashMatrices;
    for (const auto& file: flist) {
        cs225::PNG image;
        image.readFromFile(file);
        minHashMatrices.emplace_back(image, numTiles, k, h);
    }
    
    std::vector<std::tuple<int, int, int>> graph;
    for (size_t i = 0; i < minHashMatrices.size(); ++i) {
        for (size_t j = i + 1; j < minHashMatrices.size(); ++j) {
            int matchingTiles = minHashMatrices[i].countMatchTiles(minHashMatrices[j], threshold);
            graph.emplace_back(i, j, matchingTiles);
        }
    }
    return graph;
}

