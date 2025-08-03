/**
* @file mpht.hpp
 *
 * @brief This file defines the Minimal Perfect Hash Table (MPHT) class.
 *
 * @author Yi Chen
 * @contact: yi.chen.01@outlook.com
 * @feat: Added support for vertical minimizer sections and hashtable construction.
 *
 */

#pragma once

#include <vector>

#include "pthash.hpp"
#include "logfault.h"


typedef pthash::single_phf<pthash::murmurhash2_64,   // base hasher
        pthash::dictionary_dictionary,               // encoder type
        true                                         // minimal
>       pthash_type;

template<typename K, typename V>
class MPHT {
private:
    pthash::build_configuration config;
    pthash_type mphf;

public:
    std::vector<V> hashtable;
    MPHT();
    ~MPHT();
    uint64_t size();
    void build(const std::vector<K>& keys, const std::vector<V>& values);
    V find(K key);
    V operator[](K key);
    void save(const std::string& filename);
    void load(const std::string& filename);
};


/* Implementation */


template<typename K, typename V>
MPHT<K, V>::MPHT() {
    config.c = 6.0;
    config.alpha = 0.94;
    config.minimal_output = true;  // mphf
    config.verbose_output = false;
}

template<typename K, typename V>
MPHT<K, V>::~MPHT() = default;

template<typename K, typename V>
uint64_t MPHT<K, V>::size() {
    return hashtable.size();
}

template<typename K, typename V>
void MPHT<K, V>::build(const std::vector<K>& keys, const std::vector<V>& values) {
    if (keys.size() != values.size()) {
        throw std::invalid_argument("Keys and values must have the same size.");
    }
    // LFLOG_DEBUG << "Building MPHT with " << keys.size() << " keys and values.";
    // Build the minimal perfect hash function
    mphf.build_in_internal_memory(keys.begin(), keys.size(), config);

    // Map the values to the keys
    hashtable.resize(keys.size());
    for (uint64_t i = 0; i < size(); i++) {
        hashtable[mphf(keys[i])] = values[i];
    }
}

template<typename K, typename V>
V MPHT<K, V>::find(K key) {
    return hashtable[mphf(key)];
}

template<typename K, typename V>
V MPHT<K, V>::operator[](K key) {
    return find(key);
}

/**
 * Save the MPHT to a file in big-endian format
 * Schema: [mphf_size: 8 bytes]
 *         [mphf: mphf_size bytes]
 *         [len(size of hashtable): 8 bytes]
 *         [hashtable: len * sizeof(V) bytes]
 * @tparam K
 * @tparam V
 * @param filename
 */
template<typename K, typename V>
void MPHT<K, V>::save(const std::string& filename) {
    // Save the temporary mphf file
    essentials::save(mphf, filename.c_str());
}

/**
 * Load the MPHT from a file in big-endian format
 * Schema: [mphf_size: 8 bytes]
 *         [mphf: mphf_size bytes]
 *         [len(size of hashtable): 8 bytes]
 *         [hashtable: len * sizeof(V) bytes]
 * @tparam K
 * @tparam V
 * @param filename
 */
template<typename K, typename V>
void MPHT<K, V>::load(const std::string& filename) {
    essentials::load(mphf, filename.c_str());
}