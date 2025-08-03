/**
* @file util.cpp
 *
 * @brief This file defines utility functions for handling k-mer and minimizer data.
 *
 * This file is based on Kmer File Format (KFF) by Yoann Dufresne and has been extended by Yi Chen.
 *
 * @author Yi Chen
 * @contact: yi.chen.01@outlook.com
 * @feat: Added support for vertical minimizer sections and hashtable construction.
 *
 * @author Yoann Dufresne (Original Author)
 * @contact: yoann.dufresne0@gmail.com
 *
 */

#include "kero-cpp-api/detail/util.hpp"

uint64_t kero::get_mini_mask(uint64_t m) {
    // Note: left shift a 64-bit integer by 64 is undefined behavior
    // So the case m == 32 is handled separately
    // a literal 1 is treated as an int, so we use 1ull to make it a 64-bit integer
    return m >= 32 ? 0xFFFFFFFFFFFFFFFF : (1ull << (2 * m)) - 1;
}

uint64_t kero::mask_mini(uint64_t minimizer, uint64_t m) {
    return minimizer & get_mini_mask(m);
}

uint64_t kero::mask_mini(const uint8_t* mini_arr, uint64_t m) {
    uint64_t minimizer = 0;
    int nb_bytes_mini = static_cast<int>(((2 * m - 1) / 8) + 1);
    for (int i = 0; i < nb_bytes_mini; i++) {
        minimizer = (minimizer << 8) + mini_arr[i];
    }
    return mask_mini(minimizer, m);
}