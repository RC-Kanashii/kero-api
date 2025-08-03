/**
* @file kero_mmap.hpp
 *
 * @brief This file defines the Kero_Mmap_Accessor class for memory-mapped file access.
 *
 * @author Yi Chen
 * @contact: yi.chen.01@outlook.com
 * @feat: Added support for vertical minimizer sections and hashtable construction.
 *
 */

#ifndef KERO_MMAP_HPP
#define KERO_MMAP_HPP

#include <string>
#include <cstdint>
#include <stdexcept>

namespace kero {

    class Kero_Mmap_Accessor {
    private:
        int fd;                 // File descriptor
        uint8_t* file_ptr;      // Pointer to the mapped memory
        size_t file_size;       // Total size of the mapped file

    public:
        /**
         * @brief Construct a new Kero Mmap Accessor object and map the file into memory.
         * @param filename The path to the file to map.
         */
        Kero_Mmap_Accessor(const std::string& filename);

        /**
         * @brief Destroy the Kero Mmap Accessor object, unmapping the memory and closing the file.
         */
        ~Kero_Mmap_Accessor();

        // Disable copy and assign
        Kero_Mmap_Accessor(const Kero_Mmap_Accessor&) = delete;
        Kero_Mmap_Accessor& operator=(const Kero_Mmap_Accessor&) = delete;

        /**
         * @brief Get a pointer to the beginning of the mapped file.
         * @return const uint8_t* A read-only pointer to the mapped memory.
         */
        const uint8_t* get_ptr() const {
            return file_ptr;
        }

        /**
         * @brief Get the total size of the mapped file.
         * @return size_t The file size in bytes.
         */
        size_t get_size() const {
            return file_size;
        }
    };

} // namespace kero

#endif //KERO_MMAP_HPP