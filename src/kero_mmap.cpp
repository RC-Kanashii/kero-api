/**
* @file kero_mmap.cpp
 *
 * @brief This file defines the Kero_Mmap_Accessor class for memory-mapped file access.
 *
 * @author Yi Chen
 * @contact: yi.chen.01@outlook.com
 * @feat: Added support for vertical minimizer sections and hashtable construction.
 *
 */

#include "kero-api/kero_mmap.hpp"

// Required for mmap, fstat, etc.
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

namespace kero {

    Kero_Mmap_Accessor::Kero_Mmap_Accessor(const std::string& filename)
        : fd(-1), file_ptr(nullptr), file_size(0) {
        fd = open(filename.c_str(), O_RDONLY);
        if (fd == -1) {
            throw std::runtime_error("Mmap_Accessor: Failed to open file: " + filename);
        }

        struct stat sb{};
        if (fstat(fd, &sb) == -1) {
            close(fd);
            throw std::runtime_error("Mmap_Accessor: Failed to get file size.");
        }
        file_size = sb.st_size;

        // MAP_PRIVATE ensures that writes to the mapping are not propagated to the file.
        // It's good practice for read-only access.
        file_ptr = static_cast<uint8_t*>(mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (file_ptr == MAP_FAILED) {
            close(fd);
            throw std::runtime_error("Mmap_Accessor: Failed to map file to memory.");
        }
    }

    Kero_Mmap_Accessor::~Kero_Mmap_Accessor() {
        if (file_ptr != nullptr) {
            munmap(file_ptr, file_size);
        }
        if (fd != -1) {
            close(fd);
        }
    }

} // namespace kero