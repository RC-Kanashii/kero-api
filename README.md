# C++ Library for Kero Files

This is a C++ library for reading and writing kero files. It provides a high-level API for simple k-mer iteration and a low-level API for fine-grained control over the file structure.

## High-Level API: Easy K-mer Enumeration

The high-level `Kero_reader` API allows you to iterate through all k-mers in a file without managing the underlying section-based format.

### Iterating Through K-mers

Use a simple while loop with `next_kmer` to process each k-mer and its associated data.

```cpp
#include "kero-cpp-api/kero_io.hpp"

Kero_reader reader("my_file.kero");
uint8_t* kmer;
uint8_t* data;

// next_kmer returns true if a k-mer was successfully read
while (reader.next_kmer(kmer, data)) {
    // Use the kmer and its data
}
```

### Accessing File Properties

Key properties from the file's value sections ('v') are easily accessible.

```cpp
// Get variables like k-mer or minimizer size
uint64_t k = reader.get_var("k");
uint64_t m = reader.get_var("m");
uint64_t data_size = reader.get_var("data_size");

// Get the nucleotide encoding table (in order A, C, G, T)
uint8_t* encoding = reader.get_encoding();
```

## Low-Level API: Full Control

The low-level API gives you direct access to read and write individual sections of a kero file using the `Kero_file` object.

### File Handling

Open files in read ("r") or write ("w") mode.

```cpp
// Open a file for reading
Kero_file infile("path/to/file.kero", "r");
infile.close();

// Open a file for writing
// Note: The directory path must exist beforehand.
Kero_file outfile("path/to/new_file.kero", "w");
outfile.close(); // Automatically writes the footer and index
```

### Writing a Kero File

#### 1. Write the Header

Start by defining the nucleotide encoding and optional metadata.

```cpp
// write_encoding must be called before other data.
// Example: A=0, C=1, G=2, T=3
outfile.write_encoding(0, 1, 2, 3);

// Optionally, write metadata
std::string meta = "My file metadata";
outfile.write_metadata(meta.length(), (uint8_t*)meta.c_str());
```

#### 2. Write Data Sections

Create and write one or more data sections.

##### Value Section ('v')

Used for storing key-value metadata like k.

```cpp
Section_GV sgv(outfile);
sgv.write_var("k", 31);
sgv.write_var("m", 15);
sgv.write_var("data_size", 1);
sgv.close();
```

##### Raw Sequence Section ('r')

For writing blocks of overlapping k-mers.

```cpp
Section_Raw sr(outfile);

// Assumes 'byte_seq' is an encoded sequence and 'counts' is the data array
// Example: sequence "GATTACA" contains one 7-mer
sr.write_compacted_sequence(byte_seq, 7, counts);

sr.close();
```

##### Vertical Minimizer Section ('M')

This specialized section stores all super-k-mers for a single minimizer in a compressed, columnar format. Data is buffered, and only written to disk when the section is closed.

```cpp
Section_Minimizer sm(outfile);

// 1. Define the minimizer for this entire section
uint8_t* minimizer_bytes = encode("ACA"); // Assumes an encode function
sm.write_minimizer(minimizer_bytes);

// 2. Add a super-k-mer's data. This consists of the sequence parts
//    *before* and *after* the minimizer.
//    (e.g., for "GATT-ACA-G", the sequence to write is "GATTG")
std::string seq_without_mini = "GATTG";
uint8_t* seq_bytes = encode(seq_without_mini);
uint64_t minimizer_position_in_superkmer = 4;

sm.write_compacted_sequence_without_mini(
    seq_bytes,
    seq_without_mini.length(),
    minimizer_position_in_superkmer,
    counts // Data for the k-mers in this super-k-mer
);

// 3. Close the section to write all buffered data to disk
sm.close();
```

#### 3. Automatic Footer and Indexing

When the `Kero_file` is closed, the library automatically generates a footer containing:

- **A Hashtable Section ('h')**: Contains a Minimal Perfect Hash Function (MPHF) that provides fast lookups for 'M' sections based on their minimizer.

- **An Index Section ('i')**: Maps the locations of all sections in the file.

This makes the output file efficiently seekable.

### Reading a Kero File

Process a file by reading its sections sequentially.

#### 1. Detect Section Type

Check the type character to determine how to read the next section.

```cpp
char section_type = infile.read_section_type();
```

#### 2. Read Sections

Based on the type, construct the appropriate section reader to process its contents.

##### Value Section ('v')

```cpp
Section_GV sgv(&infile);
uint64_t k = sgv.vars["k"];
sgv.close();
// The values are also loaded into infile.global_vars
```

##### Raw Sequence Section ('r')

```cpp
Section_Raw sr(&infile);

uint8_t* seq_buf = new uint8_t[...];
uint8_t* data_buf = new uint8_t[...];

while(sr.remaining_blocks > 0) {
    uint64_t nb_kmers = sr.read_compacted_sequence(seq_buf, data_buf);
    // Use the sequence block and data
}
sr.close();
```

##### Vertical Minimizer Section ('M')

Reading a minimizer section populates the minimizer and allows iteration through its associated super-k-mers.

```cpp
Section_Minimizer sm(&infile);

// The minimizer for this section is now in sm.minimizer
uint8_t* minimizer = sm.minimizer;
uint64_t minimizer_pos;

while(sm.remaining_blocks > 0) {
    // Reads one super-k-mer's data (sequence parts and k-mer data)
    sm.read_compacted_sequence_without_mini(seq_buf, data_buf, minimizer_pos);
    // You can reconstruct the full super-k-mer using the minimizer and its position
}
sm.close();
```

### Index and Hashtable Handling

When a file is opened in read mode, its index ('i') and hashtable ('h') sections are automatically discovered and loaded into memory to enable fast navigation. This process is transparent to the user.