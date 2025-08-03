/**
* @file util.hpp
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

#pragma once

// Utils

#include <cstdint>
#include <cstddef>

namespace kero {

    template<typename T>
    void store_big_endian(uint8_t *buff, size_t size, const T &data) {
        for (int b = size - 1; b >= 0; --b) {
            *buff++ = data >> (8 * b);
        }
    }

    template<typename T>
    void load_big_endian(const uint8_t *buff, size_t size, T &data) {
        data = 0;
        for (int b = 0; b < size; b++) {
            data <<= 8;
            data |= buff[b];
        }
    }

    uint64_t get_mini_mask(uint64_t m);

    uint64_t mask_mini(uint64_t minimizer, uint64_t m);

    uint64_t mask_mini(const uint8_t* mini_arr, uint64_t m);
}