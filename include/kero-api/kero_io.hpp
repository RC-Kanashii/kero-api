/**
 * @file kero_io.hpp
 *
 * @brief This file defines the classes needed to create and read kero files.
 * This file contains a low level and a high level APIs.
 * The Kero_file class is the base class for the low level API.
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

#include <fstream>
#include <unordered_map>
#include <map>
#include <vector>

#include "kero-api/detail/mpht.hpp"
#include "ic.h"

#ifdef _WIN32
#include <iso646.h>
using uint = unsigned long;
#endif


#ifndef KERO_IO
#define KERO_IO



class Section_GV;
class Section_Raw;
class Section_Minimizer_Old;
class Section_Minimizer;
class Section_Index;
class Section_Hashtable;
class Kero_reader;

/**
 * This class is the central class for the low level kero file API.
 *
 * Kero_file directly contains the functions needed to create/read the kero file headers.
 * The class also contains functions to open section objects (one per section type).
 *
 */
class Kero_file {
public:
	std::fstream fs;
	unsigned long current_position;

	bool writing_started;
	uint8_t * file_buffer;
	unsigned long file_size;
	unsigned long buffer_size;
	unsigned long max_buffer_size;
	unsigned long next_free;
	bool delete_on_destruction;

	bool tmp_closed;	

	/**
	 * Read encoding from file and save it to the public argument "encoding".
	 */
	void read_encoding();

	/**
	 * Read the metadata size from the file and save the value in metadate_size.
	 */
	void read_size_metadata();

	void write_footer();
	void footer_discovery();
	void index_discovery();
	void read_index(long position);

public:
	std::string filename;
	uint8_t major_version;
	uint8_t minor_version;
	bool uniqueness;
	bool canonicity;
	
	// True when the header have been read/write and the fs pointer is after.
	bool header_over;
	bool footer_discovery_ended;
	bool index_discovery_ended;

	unsigned long end_position;
	
	bool is_writer;
	bool is_reader;

	Section_GV * footer;
	std::vector<Section_Index *> index;

	bool indexed;
	std::map<long, char> section_positions;

	// encoding:        A:0  C:1 G:3 T:2
	uint8_t encoding[4] = {0, 1, 3, 2};
	
	uint32_t metadata_size = 0;

	std::unordered_map<std::string, uint64_t> global_vars;

    // For minimizer section registration and hashtable construction
    std::vector<uint64_t> mini_list;
    std::vector<uint64_t> mini_pos;

	// --- Filesystem functions ---
	/** Open the file filename with the mode.
	 * mode must be chosen in the set of values {r: read, w: write}
	 *
	 * @param filename The path to the file to construct/read.
	 * @param mode Opening mode of the file. w for writing, r for reading.
   *
	 */
	Kero_file(const std::string filename, const std::string mode);
	/** Reopen from the beginning the KERO file defined at the construction.
	 * The opening mode can differ from the original one.
	 * 
	 * @param mode r for reading, w for writing (overwrite any previous file).
	 */
	void open(const std::string mode);
	/**
	 * Destroy properly a Kero object, closing the file behind it.
	 */
	~Kero_file();
	/**
	 * Close the file.
	 * 
	 * @param write_buffer Write the buffer in writing mode. If set to false, the buffer is never
	 * saved on disk and the file will be deleted on object destruction.
	 */
	void close(bool write_buffer=true);

	/** Read a slice of the file and copy it into bytes
	 * @param bytes Bytes array already allocated by the user (minimum size bytes)
	 * @param Size in bytes to copy
	 * 
	 * @return true if everything is ok.
	 */
	void read(uint8_t * bytes, unsigned long size);
	/** Writes size bytes into the file buffer. Written on disk if size > 1MB or on file closing.
	 * @param bytes Bytes to write
	 * @param Size in bytes to copy
	 */
	void write(const uint8_t * bytes, unsigned long size);
	/** Writes size bytes at position and returns at the same position than initially.
	 * @param bytes Bytes to write.
	 * @param size Number of bytes to write
	 * @param position Where to write in the file (will overwrite the previous data)
	 */
	void write_at(const uint8_t * bytes, unsigned long size, unsigned long position);
	/** Give the current position in the file
	 * 
	 * @return position from the beginning of the file.
	 */
	unsigned long tellp();
	/** Relative jump. Can be negative.
	 * @param size Number of bytes to jump
	 **/
	void jump(long size);
	/** Jump to a specific position in the file.
	 * 
	 * @param position Position to reach in the file.
	 * @param from_end If true, position from the end of the file.
	 */
	void jump_to(unsigned long position, bool from_end=false);

	/**
	 * In writing mode, the KERO files are indexed by default.
	 * This function allows the deactivation of the footer index generation.
	 */
	void tmp_close();
	/**
	 * Reopen the file and point to the end of it.
	 * This function is automatically called on write events if the file has been temporarily closed.
	 */
	void reopen();
	

	// --- Index related ---

	void set_indexation(bool indexed);
	/**
	 * Register a section into index
	 */
	void register_position(char section_type);
	/**
	 * Release the file pointer by temporarily close the file stream.
	 * The usage of this function increase the execution time.
	 * It should be use only when the user need a large amount of simulataneous kero files.
	 */

	// --- header functions ---
	/**
	 * Set the encoding used to compact the nucleotides into 2-bits values.
	 * Only the two lower bits of each uint8_t will be used.
	 * The 4 2-bits values must be diferent to each other.
	 *
	 * @param a encoding for A or a letters.
	 * @param c encoding for C or c letters.
	 * @param g encoding for G or g letters.
	 * @param t encoding for T or t letters.
	 */
	void write_encoding(uint8_t a, uint8_t c, uint8_t g, uint8_t t);
	/**
	 * Set the encoding used to compact the nucleotides into 2-bits values.
	 * Only the two lower bits of each uint8_t will be used.
	 * The 4 2-bits values must be diferent to each other.
	**/
	void write_encoding(uint8_t * encoding);
	/** Set the uniqueness of the kmers in the file.
	  * true means that no kmer will be present 2 times in the file.
	  **/
	void set_uniqueness(bool uniqueness);
	/** Set the cononicity of the kmers in the file.
	  * true means that if a kmer is present in the file, its reverse complement is absent.
	  **/
	void set_canonicity(bool canonicity);
	/**
	 * Write the metadata entered by the user.
	 * Also write all the needed values like metadata section size into the file.
	 *
	 * @param size The size of data array.
	 * @param data An array of data. Can either be plain text or binary data.
	 *
	 */
	void write_metadata(uint32_t size, const uint8_t * data);
	/**
	 * Read the next size Bytes in the file and return them into the data array.
	 * This function is used to read the metadata field.
	 * The size of the metadata is saved in the metadata_size field.
	 *
	 * @param data Array filled with the file content (up to size).
	 *
	 */
	void read_metadata(uint8_t * data);
	/**
	 * Function called to complete a header reading or writing before the first section io
	 */
	void complete_header();


	// --- general section ---
	/**
	 * Read the next Byte and return it.
	 * If the file pointer is correctly aligned with the beginning of a section, this Byte is a char that represents the section type.
	 *
	 * @return A char corresponding to the type of the following section.
	 *
	 */
	char read_section_type();
	/**
	 * Jump over the next section if it is section containing kmers.
	 *
	 * @return True if a section is skipped, false if the section can't be skipped, eof is reached or the file is not read ready.
	 */
	bool jump_next_section();


    // --- Section Hashtable related ---
    /**
     * Register a minimizer section into the hashtable
     */
    void register_minimizer_section(uint64_t minimizer);
};


class Section {
	protected:
		uint8_t * buffer;
	public:
		Kero_file * file;
		long beginning;

		bool active_buffer;
		long max_buffer_size;

		Section(Kero_file * file);
		virtual ~Section() {};
		
		/** Copy the current section in the file pointed by the function.
		 * 
		 * @param file The file where to copy the section
		 **/
		virtual void copy(Kero_file * file) {};

    virtual void close();
};


class SectionBuilder {
public:
	/** Build the next section in the file.
	 * The file pointer must be aligned on the first Byte of a section.
	 * 
	 * @param file Pointer to a file to read.
	 * 
	 * @return The correct section from the file.
	 **/
	static Section * build(Kero_file * file);
};


/**
 * File manipulator for Global Variable sections.
 *
 */
class Section_GV : public Section {
private:
	void read_section();
	void read_var();

	// friend class Kero_file;

public:
	uint64_t nb_vars;
	/**
	 * A map containing all the declared variables for this section.
	 * If the file is opened in r mode, all the variables of this section have been loaded during the object construction.
	 */
	std::map<std::string, uint64_t> vars;

	Section_GV(Kero_file * file);
	/**
	 * Write a variable into the file.
	 *
	 * @param var_name String name of the variable.
	 * @param value The 64 bits value of the constant.
	 *
	 */
	void write_var(const std::string & var_name, uint64_t value);
	/** Copy the current section in the file pointed by the function.
	 * 
	 * @param file The file where to copy the section
	 **/
	void copy(Kero_file * file);
	/**
	 * Closes the section.
	 * If w mode, go back to the beginning of the section to write the correct number of variables.
	 *
	 */
	void close();
};

class Block_section_reader {
protected:
	uint8_t nb_kmers_bytes;

	friend class Kero_reader;

public:
	/**
	 * The number of blocks in the section.
	 * mode r: filled during the construction of the object.
	 * mode w: Increased each time a block is added. Used on close to write the number at the beginning of the section.
	 */
	uint64_t nb_blocks;
	uint64_t remaining_blocks;
	uint64_t k;
	uint64_t max;
	uint64_t data_size;

	virtual ~Block_section_reader(){};

	static Block_section_reader * construct_section(Kero_file * file);
	/**
	 * Read the next block of the section.
	 * The sequence of the block is pushed in the seq array and the data in the data array.
	 * These arrays must be already allocated.
	 * 
	 * @param seq Array filled with the compacted sequence (2 bit / nucl).
	 * @param data Array filled with data linked to the kmers of the sequence.
	 *
	 * @return The number of kmers in the sequence.
	 */
	virtual uint64_t read_compacted_sequence(uint8_t* seq, uint8_t* data) {return 0;}
	/**
	 * Read the next block of the section.
	 * The sequence of the block is pushed in the seq_data array, immediatly followed by data.
	 * The array must be already allocated.
	 * 
	 * @param seq_data Array filled with the compacted sequence (2 bit / nucl) immediatly followed
	 * with data linked to the kmers of the sequence.
	 *
	 * @return The number of kmers in the sequence.
	 */
	virtual uint64_t read_compacted_sequence(uint8_t* seq_data) {return 0;}
	/**
	 * Jumb over the next block of the section.
	 */
	virtual void jump_sequence() {}
	/**
	 * Jump over the full section
	 */
	void jump_section() {
		while (this->remaining_blocks > 0)
			this->jump_sequence();
	}
};

/**
 * File manipulator for Index sections.
 *
 */
class Section_Index : public Section {
private:

	friend class Kero_file;
	using Section::copy;

public:
	std::map<int64_t, char> index;
	int64_t next_index;

	Section_Index(Kero_file * file);
	~Section_Index(){};
	void register_section(char section_type, int64_t index);
	void set_next_index(int64_t index);
	void close();
};

/**
 * File manipulator for Raw sections.
 *
 */
class Section_Raw: public Section, public Block_section_reader {
private:
	friend class Kero_file;
	friend class Block_section_reader;

	uint32_t read_section_header();

public:
	Section_Raw(Kero_file * file);
	~Section_Raw(){};

	/**
	 * Read the next block of the section.
	 * The sequence of the block is pushed in the seq array and the data in the data array.
	 * These arrays must be already allocated.
	 * 
	 * @param seq Array filled with the compacted sequence (2 bit / nucl).
	 * @param data Array filled with data linked to the kmers of the sequence.
	 *
	 * @return The number of kmers in the sequence.
	 *
	 */
	uint64_t read_compacted_sequence(uint8_t* seq, uint8_t* data);
	/**
	 * Read the next block of the section.
	 * The sequence of the block is pushed in the seq_data array, immediatly followed by data.
	 * The array must be already allocated.
	 * 
	 * @param seq_data Array filled with the compacted sequence (2 bit / nucl) immediatly followed
	 * with data linked to the kmers of the sequence.
	 *
	 * @return The number of kmers in the sequence.
	 */
	uint64_t read_compacted_sequence(uint8_t* seq_data);
	/**
	 * Write a block containing the sequence seq with kmer data associated.
	 * 
	 * @param seq A compacted sequence (2 bit / nucl).
	 * @param seq_size Size of the sequence (in nucleotides).
	 * @param data Data array of the kmers in the sequence.
	 *
	 */
	void write_compacted_sequence(uint8_t* seq, uint64_t seq_size, uint8_t * data_array);
	/**
	 * Jumb over the next block of the section.
	 */
	void jump_sequence();
	/** Copy the current section in the file pointed by the function.
	 * 
	 * @param file The file where to copy the section
	 **/
	void copy(Kero_file * file);
	/**
	 * Close the section.
	 * If w mode, go back to the beginning of the section to write the correct number of blocks.
	 *
	 */
	void close();
};


/**
 * File manipulator for Minimizer_Vertical sections.
 *
 */
class Section_Minimizer : public Section, public Block_section_reader {
private:
	// Buffers
    std::vector<uint64_t> n_value_buffer;  // the number of k-mers in each block
    std::vector<uint64_t> m_idx_buffer;    // the minimizer index for each block
    std::vector<uint8_t> seq_buffer;       // the sequence buffer for each block
	std::vector<uint8_t> data_buffer;      // the data buffer for each block

	// For sequence reading and writing
	uint64_t cur_skmer_idx;                // current super k-mer index to read or write
	uint64_t last_n_pos;				   // last n position to read or write
	uint64_t last_m_idx_pos;			   // last m_idx position to read or write
	uint64_t last_seq_pos;				   // last seq position to read or write
	uint64_t last_data_pos;    			   // last data position to read or write

    void read_section_header();
    void write_section_header();
    void write_columns();
    void backfill_column_offsets();

public:
    Section_Minimizer(Kero_file* file);
	Section_Minimizer& operator= (Section_Minimizer && smv);
    ~Section_Minimizer();

    uint64_t k;                           // k value of kmer
    uint64_t max;                         // max block size
    uint64_t data_size;                   // data size
	uint64_t m;                           // m value of minimizer
	uint8_t* minimizer;                   // minimizer

	// Useful variables
    uint8_t nb_bytes_mini;                 // the number of bytes used to store the minimizer
    uint8_t mini_pos_bytes;                // the number of bytes used to store the minimizer position

	// Relative offsets of the columns to the start of the section
    uint64_t n_col_offset;
    uint64_t m_idx_col_offset;
    uint64_t seq_col_offset;
    uint64_t data_col_offset;
	uint64_t start_pos;

	// Public methods
    void write_minimizer(uint8_t* minimizer);
    void write_compacted_sequence_without_mini(uint8_t* seq, uint64_t seq_size, uint64_t mini_pos, uint8_t* data_array);
    void write_compacted_sequence(uint8_t* seq, uint64_t seq_size, uint64_t mini_pos, uint8_t* data_array);
    void add_minimizer(uint64_t nb_kmer, uint8_t* seq, uint64_t mini_pos);
    uint64_t read_compacted_sequence(uint8_t* seq, uint8_t* data);
    uint64_t read_compacted_sequence(uint8_t* seq_data);
    uint64_t read_compacted_sequence_without_mini(uint8_t *seq, uint8_t *data, uint64_t &mini_pos);
	void copy(Kero_file * file);
    void jump_sequence();
    void close();

	/**
	 * @brief Reads and decompresses all column data (n, m_idx, data) from a memory-mapped file.
	 * This method is designed to be called once to pre-cache data for parallel access.
	 *
	 * @param mmap_ptr Pointer to the start of the memory-mapped file.
	 */
	void precache_columns_from_mmap(const uint8_t* mmap_ptr);
};


/**
 * File manipulator for Hashtable sections.
 *
 * Schema:
 * ascii(h): 1B
 * nb_mphf: 8B
 * mphf: nB
 * nb_hashtable: 8B
 * hashtable: 8*n B
 *
 */
class Section_Hashtable : public Section {
private:
    friend class Kero_file;
    using Section::copy;

    static constexpr uint64_t BUFF_CHUNK_SIZE = 1 * 1024 * 1024;  // 1MB

public:
    MPHT<uint64_t, uint64_t> mpht;
    std::vector<uint64_t> minimizers;
    std::vector<uint64_t> positions;
	uint64_t nb_mphf{};  // size of the mphf in bytes

    explicit Section_Hashtable(Kero_file * file);
    ~Section_Hashtable() override;
    void reg_sm(uint64_t minimizer, uint64_t index);
    void close();
};


class Kero_reader {
private:
	// Space alocated for copying the current kmer given to the user.
	uint8_t * current_kmer;
	// Current sequence
	uint8_t * current_seq_data;
	// uint8_t * current_sequence;
	// Current sequence shifted to match the 4 different alignements
	uint8_t ** current_shifts;
	// Size in nucleotides of current sequence
	uint64_t current_seq_nucleotides;
	// Size in bytes of current sequence
	uint64_t current_seq_bytes;
	// Number of kmers in the current sequence
	uint64_t current_seq_kmers;
	// Number of kmer remaining in the current block
	uint64_t remaining_kmers;
	// Data array for the current block
	// uint8_t * current_data;
	// Section currently read.
	Block_section_reader * current_section;
	// Remaining blocks before end of the section
	uint64_t remaining_blocks;


	void read_until_first_section_block();
	void read_next_block();

public:
	uint64_t k;
	uint64_t data_size;
	uint64_t max;
	
	Kero_file * file;
	
	Kero_reader(std::string filename);
	~Kero_reader();

	bool has_next();
	uint64_t next_block(uint8_t* & sequence, uint8_t* & data);
	bool next_kmer(uint8_t* & sequence, uint8_t* & data);

	uint64_t get_var(std::string name);
	uint8_t * get_encoding();
};

uint64_t bytes_from_bit_array(uint64_t bits_per_elem, uint64_t nb_elem);

#endif
