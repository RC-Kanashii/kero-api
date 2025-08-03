/**
* @file kero_io.cpp
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

#include <cstdint>
#include <iostream>
#include <cassert>
#include <cstring>
#include <sstream>
#include <cmath>

#include <map>
#include <vector>

#include "kero-api/kero_io.hpp"
#include "kero-api/detail/util.hpp"
#include "ic.h"

using namespace std;
using namespace kero;

#define KERO_VERSION_MAJOR 0
#define KERO_VERSION_MINOR 1

uint64_t bytes_from_bit_array(uint64_t bits_per_elem, uint64_t nb_elem) {
	if (bits_per_elem == 0 or nb_elem == 0)
		return 0;
	else
		return ((bits_per_elem * nb_elem - 1) / 8) + 1;
}

static void leftshift8(uint8_t * bitarray, size_t length, size_t bitshift);
static void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift);
static uint8_t fusion8(uint8_t left_bits, uint8_t right_bits, size_t merge_index);
static inline size_t p4nenc_bound(size_t n, size_t size);
static inline size_t round_up(size_t n, size_t a);


// ----- Open / Close functions -----

Kero_file::Kero_file(const string filename, const string mode) {
	// Variable init
	this->filename = filename;
	
	this->is_writer = false;
	this->is_reader = false;
	
	this->writing_started = false;
	this->next_free = 0;
	this->buffer_size = 1 << 10; // 1 KB
	// this->buffer_size = 1 << 4;
	this->max_buffer_size = 1 << 20; // 1 MB
	// this->max_buffer_size = 1 << 6;
	this->file_buffer = new uint8_t[this->buffer_size];
	this->file_size = 0;
	this->delete_on_destruction = false;

	this->open(mode);
}

void Kero_file::open(string mode) {
	this->writing_started = false;
	this->current_position = 0;

	// Determine the mode and open the file
	if (mode[0] == 'w') {
		this->is_writer = true;
		this->file_size = 0;
		this->next_free = 0;
	} else if (mode[0] == 'r') {
		this->is_reader = true;
		// If no info on the file
		if (this->file_size == 0 and this->next_free == 0) {
			// Open the fp
			this->fs.open(this->filename, fstream::binary | fstream::in);
			if (this->fs.fail()) {
				std::string msg = "Cannot open file " + this->filename;
				throw std::runtime_error(msg);
			}
			// Compute the file length
			long position = this->fs.tellp();
			this->fs.seekg(0, this->fs.end);
			this->file_size = (long)(this->fs.tellp()) - position;
			// Go back to the beginning
			this->fs.seekg(0, this->fs.beg);
		}
	} else {
		cerr << "Unsupported mode " << mode << endl;
		exit(1);
	}

	this->tmp_closed = false;
	this->header_over = false;
	this->indexed = false;
	this->footer = nullptr;
	this->footer_discovery_ended = true;

	// Write the signature and the version at the beginning of the file
	if (this->is_writer) {
		uint8_t default_encoding = 0b00011110;
		// Signature
		uint8_t buff[] = {	'K', 'E', 'R', 'O',
							KERO_VERSION_MAJOR, KERO_VERSION_MINOR,
							default_encoding,
							0 /*uniqueness*/, 0 /*canonicity*/
						};

		this->write(buff, 9);

		this->indexed = true;
		this->end_position = 0;
	}
	// Read the header
	else if (this->is_reader) {
		// Header integrity marker
		uint8_t buff[4];
		this->read(buff, 4);
		if (buff[0] != 'K' or buff[1] != 'E' or buff[2] != 'R' or buff[3] != 'O') {
			cerr << "Absent KERO signature at the beginning of the file." << endl;
			cerr << "Please check that the file is not corrupted" << endl;
			throw "Absent signature at the beginning";
		}

		// Version reading
		this->read(&this->major_version, 1);
		this->read(&this->minor_version, 1);
		if (KERO_VERSION_MAJOR < this->major_version or (KERO_VERSION_MAJOR == this->major_version and KERO_VERSION_MINOR < this->minor_version)) {
			cerr << "The software version " << (uint)KERO_VERSION_MAJOR << "." << (uint)KERO_VERSION_MINOR << " can't read files writen in version " << (uint)this->major_version << "." << (uint)this->minor_version << endl;
			throw "Unexpected version number";
		}
		// Encoding load
		this->read_encoding();
		// Read global flags
		uint8_t flag;
		this->read(&flag, 1);
		this->uniqueness = flag != 0;

		this->read(&flag, 1);
		this->canonicity = flag != 0;
		// Read metadata size
		this->read(buff, 4);
		load_big_endian(buff, 4, this->metadata_size);


		// Footer integrity marker
		unsigned long saved_position = this->tellp();
		this->jump_to(4, true);
		this->end_position = this->tellp();
		this->read(buff, 4);
		this->jump_to(saved_position);
		if (buff[0] != 'K' or buff[1] != 'E' or buff[2] != 'R' or buff[3] != 'O') {
			cerr << "Absent KERO signature at the end of the file." << endl;
			cerr << "Please check that the file is not corrupted" << endl;
			throw "Absent signature at the end";
		}

		// Back to the start
		this->footer_discovery_ended = false;
		// Discover footer
		this->footer_discovery();
		this->index_discovery();
	}
}

void Kero_file::close(bool write_buffer) {
	if (this->is_writer) {
		// Write the index
		if (this->indexed)
			this->write_footer();
		// Write the signature
		char signature[] = {'K', 'E', 'R', 'O'};
		this->write((uint8_t *)signature, 4);
		
		// Write the end of the file
		if (write_buffer) {
			// The file was never opened
			if (not this->writing_started) {
				this->writing_started = true;
				this->fs.open(this->filename, fstream::binary | fstream::out);
			} else if (this->tmp_closed) {
				this->reopen();
			}
			// Write the buffer
			this->fs.write((char *)this->file_buffer, this->next_free);
			if (this->fs.fail()) {
				cerr << "Filesystem problem during buffer disk saving" << endl;
				exit(1);
			}
			this->file_size += this->next_free;
			this->next_free = 0;
		} else {
			this->delete_on_destruction = true;
		}

		// cout << "delete_on_destruction " << delete_on_destruction << endl;
		// cout << this->filename << endl;
		if (this->fs.is_open())
			this->fs.close();
	}
	else if (this->is_reader) {
		if (fs.is_open()) {
			fs.close();
		}
	}

	this->tmp_closed = false;
	this->is_writer = false;
	this->is_reader = false;
}


Kero_file::~Kero_file() {
	this->close();

	delete[] this->file_buffer;
	if (this->delete_on_destruction and this->file_size > 0) {
		remove(this->filename.c_str());
	}

	if (this->footer != nullptr)
		delete this->footer;

	for (Section_Index * si : this->index)
		delete si;
}


void Kero_file::set_indexation(bool indexed) {
	if (this->is_writer)
		this->indexed = indexed;
}


void Kero_file::register_position(char section_type) {
	if (this->is_writer and this->indexed) {
		this->section_positions[this->tellp()] = section_type;
	}
}


void Kero_file::complete_header() {
	if (this->header_over)
		return;

	// If the metadata has not been read, jump over
	if (this->is_reader) {
		this->jump(this->metadata_size);
	}

	// If metadata has not been write, write a 0 byte one.
	else if (this->is_writer) {
		this->write_metadata(0, nullptr);
	}

	this->header_over = true;
}


void Kero_file::footer_discovery() {
	long current_pos = this->tellp();

	// Look at the footer
	this->jump_to(23, true);
	// Try to extract the footer size
	stringstream ss;
	char c = 'o';
	for (uint i=0 ; i<11 ; i++) {
		this->read((uint8_t *)&c, 1);
		ss << c;
	}
	if (ss.str().compare("footer_size") != 0) {
		return;
	}
	this->jump(1); // remove the '\0'

	uint64_t size = 0;
	uint8_t buff[8];
	this->read(buff, 8);
	load_big_endian(buff, 8, size);
	// Jump to value section start
	this->jump_to(size+3, true);

	this->footer = new Section_GV(this);
	this->footer->close();
	this->footer_discovery_ended = true;

	this->jump_to(current_pos);
}


void Kero_file::index_discovery() {
	long current_pos = this->tellp();
	bool header_over = this->header_over;
	this->complete_header();

	// Search in footer
	if (this->footer != nullptr and this->footer->vars.find("first_index") != this->footer->vars.end()) {
		this->indexed = true;
		this->read_index((long)this->footer->vars["first_index"]);
	}

	// Search first section
	if (not this->indexed) {
		char type = this->fs.peek();
		if (type == 'i') {
			this->indexed = true;
			this->read_index(this->tellp());
		}

	}

	this->header_over = header_over;
	this->index_discovery_ended = true;

	this->jump_to(current_pos);
}


void Kero_file::read_index(long position) {
	long init_pos = this->tellp();

	while (position != 0) {
		// Move to the beginning
		this->jump_to(position);
		// read the local index content
		Section_Index * si = new Section_Index(this);
		this->index.push_back(si);
		si->close();
		// Update index position to the next index section
		if (si->next_index == 0)
			position = 0;
		else {
			position = this->tellp() + si->next_index;
		}
	}

	this->jump_to(init_pos);
}


void Kero_file::read(uint8_t * bytes, unsigned long size) {
	if (not this->is_reader) {
		cerr << "Cannot read a file in writing mode." << endl;
		exit(1);
	}

	// Read in the file
	if (this->current_position < this->file_size) {
		// Read the end of the file and the beginning of the buffer
		if (this->current_position + size > this->file_size) {
			uint64_t fs_read_size = this->file_size - this->current_position;
			this->read(bytes, fs_read_size);
			this->read(bytes + fs_read_size, size - fs_read_size);
			return;
		}
		// Read inside the file
		else {
			// File not opened
			if (not this->fs.is_open())
				this->fs.open(this->filename, fstream::binary | fstream::in);

			// long tp = this->fs.tellp();
			this->fs.read((char *)bytes, size);
			if (this->fs.fail()) {
				// cout << tp << endl;
				cerr << "Impossible to read the file " << this->filename << " on disk." << endl;
				exit(1);
			}
		}
	}
	// Read in the buffer
	else {
		// Compute the buffer positions to read
		uint64_t buffer_position = this->current_position - this->file_size;
		if (buffer_position + size > this->next_free) {
			string error = string("Read out of the file, Byte ") + to_string(this->file_size + this->next_free);
			throw out_of_range(error);
			exit(1);
		}

		memcpy(bytes, this->file_buffer + buffer_position, size);
	}
	
	this->current_position += size;
}

void Kero_file::write(const uint8_t * bytes, unsigned long size) {
	if (not this->is_writer) {
		if (this->is_reader)
			cerr << "Cannot write a file in reading mode." << endl;
		else
			cerr << "Cannot write a closed file" << endl;
		exit(1);
	}

	unsigned long buff_space = this->buffer_size - this->next_free;

	// Resize buffer
	while (buff_space < size and this->buffer_size < this->max_buffer_size) {
		// Enlarge the buffer
		this->buffer_size *= 2;
		uint8_t * next_buffer = new uint8_t[this->buffer_size];
		// Copy the previous values
		memcpy(next_buffer, this->file_buffer, this->next_free);
		buff_space = this->buffer_size - this->next_free;
		// Fill the empty part with 0s
		memset(next_buffer + this->next_free, 0, buff_space);
		// remove the previous space
		delete[] this->file_buffer;
		this->file_buffer = next_buffer;
	}

	// fill the buffer
	if (buff_space >= size) {
		memcpy(this->file_buffer + this->next_free, bytes, size);
		this->next_free += size;
	}
	// Not enought space, write the file
	else {
		// Open the file if needed
		if (not this->writing_started) {
			this->fs.open(this->filename, fstream::binary | fstream::out);
			this->writing_started = true;
		} else if (this->tmp_closed) {
			this->reopen();
		}

		this->fs.write((char*)this->file_buffer, this->next_free);
		this->fs.write((char*)bytes, size);
		this->file_size += this->next_free + size;
		this->next_free = 0;

		if (this->fs.fail()) {
			cerr << "File system error while writing " << this->filename << endl;
			exit(1);
		}
	}

	this->current_position += size;
}

void Kero_file::write_at(const uint8_t * bytes, unsigned long size, unsigned long position) {
	if (not this->is_writer) {
		if (this->is_reader)
			cerr << "Cannot write a file in reading mode." << endl;
		else
			cerr << "Cannot write a closed file" << endl;
		exit(1);
	}

	if (position > this->file_size + this->next_free) {
		cerr << "Cannot write after the last byte of the file." << endl;
		exit(1);
	}

	// Write the file on disk
	if (position < this->file_size) {
		// Only in file
		if (position + size <= this->file_size) {
			if (this->tmp_closed) {
				this->reopen();
			}
			this->fs.seekp(position);
			this->fs.write((char*)bytes, size);
			if (this->fs.fail()) {
				cerr << "File system error while writing " << this->filename << " at position " << position << endl;
				exit(1);
			}
			this->fs.seekp(this->file_size);
		}
		// On both file and buffer
		else {
			unsigned long in_file_size = this->file_size - position;
			// Write the file part
			this->write_at(bytes, in_file_size, position);
			// Write the buffer part
			this->write_at(bytes + in_file_size, size - in_file_size, position + in_file_size);
		}
	}
	// Write the buffer in RAM
	else {
		unsigned long corrected_position = position - this->file_size;
		
		// Write in the current buffer space
		if (corrected_position + size <= this->next_free) {
			memcpy(this->file_buffer + corrected_position, bytes, size);
		}
		// Spillover the buffer
		else {
			this->next_free = corrected_position;
			this->write(bytes, size);
		}
	}
}

unsigned long Kero_file::tellp() {
	return this->current_position;
}


void Kero_file::jump(long size) {
	// cout << "Jump " << this->current_position << " " << size << " / " << this->file_size << " " << this->next_free << endl;
	this->jump_to(this->current_position + size);
}

void Kero_file::jump_to(unsigned long position, bool from_end) {
	if (this->file_size + this->next_free < position) {
		cerr << "Jump out of the file." << endl;
		exit(1);
	}

	// Determine absolute position
	if (from_end) {
		position = this->file_size + this->next_free - position;
	}
	// cout << "position " << position << endl;

	// Jump into the written file
	if (position < this->file_size) {
		this->fs.seekp(position);
	}
	// Jump into the buffer
	else /*if (this->current_position < this->file_size)*/ {
		this->fs.seekg(0, this->fs.end);
	}
	this->current_position = position;
}


void Kero_file::tmp_close() {
	if (this->is_writer and this->fs.is_open()) {
		this->fs.close();
		this->fs.clear();
		this->tmp_closed = true;
	}
}


void Kero_file::reopen() {
	if (this->tmp_closed) {
		auto streammode = fstream::binary | fstream::out | fstream::in | fstream::ate;

		// Open the file
		this->fs.open(this->filename, streammode);
		this->tmp_closed = false;
	}
}


void Kero_file::write_footer() {
    // Write the hashtable
    assert(mini_list.size() == mini_pos.size());
    Section_Hashtable sh(this);
    auto len = mini_list.size();
    for (size_t i = 0; i < len; i++) {
        sh.reg_sm(mini_list[i], mini_pos[i]);
    }
    sh.close();

    // Write the index section
    Section_Index si(this);

    // Compute end position
    long position = si.beginning + 17 + 9 * this->section_positions.size();

	// Add the values
	for (auto & it : this->section_positions) {
		si.register_section(it.second, it.first - position);
	}
	si.close();

	// Write a value section to register everything
	Section_GV sgv(this);
	sgv.write_var("first_index", si.beginning);
	sgv.write_var("footer_size", 9 + 2 * (12 + 8));
	sgv.close();
}


// ----- Header functions -----

void Kero_file::write_encoding(uint8_t a, uint8_t c, uint8_t g, uint8_t t) {
	// Value masking
	a &= 0b11; c &= 0b11; g &= 0b11; t &= 0b11;

	// Verification of the differences.
	assert(a != c); assert(a != g); assert(a != t);
	assert(c != g); assert(g != t);
	assert(g != t);

	// set values
	this->encoding[0] = a;
	this->encoding[1] = c;
	this->encoding[2] = g;
	this->encoding[3] = t;

	// Write to file
	uint8_t code = (a << 6) | (c << 4) | (g << 2) | t;
	this->write_at(&code, 1, 5);
}

void Kero_file::set_uniqueness(bool uniqueness) {
	uint8_t bit_uniq = uniqueness ? 1 : 0;
	this->write_at(&bit_uniq, 1, 6);
}
void Kero_file::set_canonicity(bool canonicity) {
	uint8_t bit_canon = canonicity ? 1 : 0;
	this->write_at(&bit_canon, 1, 7);
}

void Kero_file::write_encoding(uint8_t * encoding) {
	this->write_encoding(encoding[0], encoding[1], encoding[2], encoding[3]);
}

void Kero_file::read_encoding() {
	uint8_t code, a, c, g, t;
	// Get code values
	this->read(&code, 1);

	// Split each nucleotide encoding
	this->encoding[0] = a = (code >> 6) & 0b11;
	this->encoding[1] = c = (code >> 4) & 0b11;
	this->encoding[2] = g = (code >> 2) & 0b11;
	this->encoding[3] = t = code & 0b11;

	// Verification of the encoding
	if (a == c or a == g or a == t or c == g or c == t or g == t) {
		throw "Wrong encoding. The 4 2-bits values must be different.";
	}
}

void Kero_file::write_metadata(uint32_t size, const uint8_t * data) {
	if (this->header_over) {
		cerr << "The metadata have to be written prior to other content." << endl;
		exit(1);
	}

	uint8_t buff[4];
	store_big_endian(buff, 4, size);
	this->write(buff, 4);
	this->write(data, size);

	this->header_over = true;
}


void Kero_file::read_metadata(uint8_t * data) {
	this->read(data, this->metadata_size);
	this->header_over = true;
}

bool Kero_file::jump_next_section() {
	if (not is_reader)
		return false;
	char section_type = read_section_type();
	if (this->current_position == this->file_size + this->next_free)
		return false;
	if (section_type == 'r' or section_type == 'm' or section_type == 'M') {
		Block_section_reader * section = Block_section_reader::construct_section(this);
		section->jump_section();
		delete section;
		return true;
	}
	return false;
}


// ----- Sections -----

char Kero_file::read_section_type() {
	// Verify that header has been read.
	if (not this->header_over) {
		this->complete_header();
	}

	if (this->current_position < this->file_size) {
		return this->fs.peek();
	}
	else {
		return (char)this->file_buffer[this->current_position - this->file_size];
	}
}

void Kero_file::register_minimizer_section(uint64_t minimizer) {
    if (this->is_writer and this->indexed) {
        this->mini_list.push_back(minimizer);
        this->mini_pos.push_back(this->tellp());
    }
}


Section::Section(Kero_file * file) {
	this->file = file;

	if (not file->header_over and file->footer_discovery_ended) {
		file->complete_header();
	}

	this->beginning = file->tellp();
}

void Section::close() {
	this->file = nullptr;
}


Section * SectionBuilder::build(Kero_file * file) {
	char type = file->read_section_type();
	switch (type) {
		case 'i':
			return new Section_Index(file);
		case 'v':
			return new Section_GV(file);
		case 'r':
			return new Section_Raw(file);
		case 'm':
        case 'M':
			return new Section_Minimizer(file);
        case 'h':
            return new Section_Hashtable(file);
		default:
			cerr << "Unknown section " << type << "(" << (uint)type << ")" << endl;
			throw std::runtime_error("Unknown section type " + std::string(1, type));
	}
}


// ----- Global variables sections -----

Section_GV::Section_GV(Kero_file * file) : Section(file) {
	this->nb_vars = 0;
	this->file->global_vars.clear();

	if (this->file->is_reader) {
		this->read_section();
	}

	if (file->is_writer) {
		if (file->indexed)
			file->register_position('v');
		char type = 'v';
		this->file->write((uint8_t *)&type, 1);
	}
}

void Section_GV::write_var(const string & var_name, uint64_t value) {
	this->nb_vars += 1;
	this->vars[var_name] = value;
	this->file->global_vars[var_name] = value;
}

void Section_GV::read_section() {
	char type = '\0';
	this->file->read((uint8_t *)&type, 1);
	if (type != 'v')
		throw "The section do not start with the 'v' char, you can't open a Global Variable section.";

	uint8_t buff[8];
	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->nb_vars);
	for (uint64_t i=0 ; i<nb_vars ; i++) {
		this->read_var();
	}
}

void Section_GV::read_var() {
	if (file->tellp() >= file->end_position)
		throw "eof reached before the end of the variable section";

	// Name reading
	stringstream ss;
	char c = 'o';
	this->file->read((uint8_t *)&c, 1);
	while (c != '\0') {
		ss << c;
		this->file->read((uint8_t *)&c, 1);
	}

	// Value reading
	uint64_t value = 0;
	uint8_t buff[8];
	this->file->read(buff, 8);
	load_big_endian(buff, 8, value);

	// Saving
	string name = ss.str();
	this->vars[name] = value;
	this->file->global_vars[name] = value;
}

void Section_GV::copy(Kero_file * file) {
	// Remove empty variable sections
	if (this->vars.size() == 0)
		return;

	// Open the copy
	Section_GV sgv(file);
	// Copy all the variables
	for (const auto & it : this->vars) {
		sgv.write_var(it.first, it.second);
	}
	// Clos the copy
	sgv.close();
}

void Section_GV::close() {
	if (file->is_writer) {
		uint8_t buff[8];
		// write the number of block values
		store_big_endian(buff, 8, this->nb_vars);
		this->file->write(buff, 8);
		// Write the variables
		for (std::map<std::string,uint64_t>::iterator var_tuple=this->vars.begin() ; var_tuple != this->vars.end() ; var_tuple++) {
			const string & name = var_tuple->first;
			this->file->write((const uint8_t *)name.c_str(), name.length()+1);
			store_big_endian(buff, 8, var_tuple->second);
			this->file->write(buff, 8);
		}
	}

	Section::close();
}



Section_Index::Section_Index(Kero_file * file) : Section(file) {
	char type;
	uint8_t buff[8];

	this->next_index = 0;

	if (this->file->is_reader) {
		this->file->read((uint8_t *)&type, 1);
		if (type != 'i')
			throw "The section do not start with the 'i' char, you can not open an Index section.";

		uint64_t nb_vars;
		this->file->read(buff, 8);
		load_big_endian(buff, 8, nb_vars);
		for (uint64_t i=0 ; i<nb_vars ; i++) {
			int64_t idx = 0;
			this->file->read((uint8_t *)&type, 1);
			this->file->read(buff, 8);
			load_big_endian(buff, 8, idx);
			this->index[idx] = type;
		}

		if (nb_vars != this->index.size())
			throw "index collision in i section";

		this->file->read(buff, 8);
		load_big_endian(buff, 8, this->next_index);
	}
}

void Section_Index::register_section(char section_type, int64_t pos) {
	this->index[pos] = section_type;
}

void Section_Index::set_next_index(int64_t index) {
	this->next_index = index;
}

// void Section_Index::copy(Kero_file * file) {
// 	cerr << "You are trying to copy an index from a file to another. ";
// 	cerr << "As the positions can be different between the files, this operation is not allowed." << endl;
// 	exit(2);
// }

void Section_Index::close() {
	if (this->file->is_writer) {
		uint8_t buff[8];
		// Section header
		char type = 'i';
		this->file->write((uint8_t *)&type, 1);
		store_big_endian(buff, 8, this->index.size());
		this->file->write(buff, 8);
		// Write index
		for (std::map<int64_t, char>::iterator it=this->index.begin(); it!=this->index.end(); ++it) {
			// Section type
			type = it->second;
			this->file->write((uint8_t *)&type, 1);
			// Section index
			store_big_endian(buff, 8, it->first);
			this->file->write(buff, 8);
		}
		store_big_endian(buff, 8, this->next_index);
		this->file->write(buff, 8);
	}

	Section::close();
}



Block_section_reader * Block_section_reader::construct_section(Kero_file * file) {
	// Very and complete if needed the header
	file->complete_header();

	char type = file->read_section_type();
	if (type == 'r') {
		return new Section_Raw(file);
	} else if (type == 'M') {
        return new Section_Minimizer(file);
    } else
		return nullptr;
}


// ----- Raw sequence section -----

Section_Raw::Section_Raw(Kero_file * file) : Section(file){
	if (file->global_vars.find("k") == file->global_vars.end())
		throw "Impossible to read the raw section due to missing k variable";
	if(file->global_vars.find("max") == file->global_vars.end())
		throw "Impossible to read the raw section due to missing max variable";
	if(file->global_vars.find("data_size") == file->global_vars.end())
		throw "Impossible to read the raw section due to missing data_size variable";
	
	uint64_t k = file->global_vars["k"];
	uint64_t max = file->global_vars["max"];
	uint64_t data_size = file->global_vars["data_size"];

	this->nb_blocks = 0;

	this->k = k;
	this->max = max;
	this->data_size = data_size;

	// Computes the number of bytes needed to store the number of kmers in each block
	uint64_t nb_bits = static_cast<uint64_t>(ceil(log2(max)));
	this->nb_kmers_bytes = static_cast<uint8_t>(bytes_from_bit_array(nb_bits, 1));

	if (file->is_reader) {
		this->read_section_header();
	}

	if (file->is_writer) {
		if (file->indexed)
			file->register_position('r');
		char type = 'r';
		this->file->write((uint8_t *)&type, 1);
		this->file->write((uint8_t *)&this->nb_blocks, 8);
	}
}

uint32_t Section_Raw::read_section_header() {
	char type;
	this->file->read((uint8_t *)&type, 1);
	if (type != 'r')
		throw "The section do not start with the 'r' char, you can't open a Raw sequence section.";

	uint8_t buff[8];
	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->nb_blocks);
	this->remaining_blocks = this->nb_blocks;

	return this->nb_blocks;
}

void Section_Raw::write_compacted_sequence(uint8_t* seq, uint64_t seq_size, uint8_t * data_array) {
	uint8_t buff[8];
	// 1 - Write nb kmers
	uint64_t nb_kmers = seq_size - k + 1;
	store_big_endian(buff, this->nb_kmers_bytes, nb_kmers);
	this->file->write(buff, this->nb_kmers_bytes);
	// 2 - Write sequence
	uint64_t seq_bytes_needed = (seq_size + 3) / 4;
	this->file->write(seq, seq_bytes_needed);
	// 3 - Write data
	uint64_t data_bytes_needed = data_size * nb_kmers;
	this->file->write(data_array, data_bytes_needed);

	this->nb_blocks += 1;
}

uint64_t Section_Raw::read_compacted_sequence(uint8_t* seq, uint8_t* data) {
	uint8_t buff[8];
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 2 - Read the sequence
	size_t seq_size = nb_kmers_in_block + k - 1;
	size_t seq_bytes_needed = (seq_size + 3) / 4;
	this->file->read(seq, seq_bytes_needed);
	// 3 - Read the data.
	uint64_t data_bytes_used = data_size * nb_kmers_in_block;
	this->file->read(data, data_bytes_used);

	this->remaining_blocks -= 1;

	return nb_kmers_in_block;
}


uint64_t Section_Raw::read_compacted_sequence(uint8_t* seq_data) {
	uint8_t buff[8];
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 2 - Read the sequence
	size_t seq_size = nb_kmers_in_block + k - 1;
	size_t seq_bytes_needed = (seq_size + 3) / 4;
	uint64_t data_bytes_used = data_size * nb_kmers_in_block;
	this->file->read(seq_data, seq_bytes_needed + data_bytes_used);

	this->remaining_blocks -= 1;

	return nb_kmers_in_block;
}


void Section_Raw::copy(Kero_file * file) {
	uint max_nucl = this->k + this->max - 1;
	uint8_t * seq_buffer = new uint8_t[(max_nucl + 3) / 4];
	uint8_t * data_buffer = new uint8_t[this->max * this->data_size];

	// Open the copy
	Section_Raw sr(file);
	// Copy all the variables
	for (uint i=0 ; i<this->nb_blocks ; i++) {
		// Read
		uint64_t size = this->read_compacted_sequence(seq_buffer, data_buffer);
		// Rewrite
		sr.write_compacted_sequence(seq_buffer, this->k + size - 1, data_buffer);
	}
	// Clos the copy
	sr.close();
	delete[] seq_buffer;
	delete[] data_buffer;
}


void Section_Raw::jump_sequence() {
	uint8_t buff[8];
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 2 - Determine the sequence size
	size_t seq_size = nb_kmers_in_block + k - 1;
	size_t seq_bytes_needed = (seq_size + 3) / 4;
	// 3 - Determine the data size
	size_t data_bytes_used = data_size * nb_kmers_in_block;
	// 4 - Jumb over the 
	file->jump(seq_bytes_needed + data_bytes_used);
	this->remaining_blocks -= 1;
}


void Section_Raw::close() {
	if (this->file->is_writer) {
		uint8_t buff[8];
		store_big_endian(buff, 8, this->nb_blocks);
		this->file->write_at(buff, 8, this->beginning + 1);
	}

	if (file->is_reader) {
		// Jump over remaining sequences of the section
		while (this->remaining_blocks > 0)
			this->jump_sequence();
	}

	Section::close();
}


/* Bitshift to the left all the bits in the array with a maximum of 7 bits.
 * Overflow on the left will be set into the previous cell.
 */
static void leftshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	if (length > 0) {
		for (uint64_t i=0 ; i<length-1 ; i++) {
			bitarray[i] = (bitarray[i] << bitshift) | (bitarray[i+1] >> (8-bitshift));
		}
		bitarray[length-1] <<= bitshift;
	}
}

/* Similar to the previous function but on the right */
static void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	if (length > 0) {
		for (uint64_t i=length-1 ; i>0 ; i--) {
			bitarray[i] = (bitarray[i-1] << (8-bitshift)) | (bitarray[i] >> bitshift);
		}
		bitarray[0] >>= bitshift;
	}
}

/* Fusion to bytes into one.
 * The merge_index higher bits are from left_bits the others from right_bits
 */
static uint8_t fusion8(uint8_t left_bits, uint8_t right_bits, size_t merge_index) {
	uint8_t mask = 0xFF << (8-merge_index);
	return (left_bits & mask) | (right_bits & ~mask);
}

/* Computes the number of bytes needed to store a bit array of size bits.
 * The result is rounded up to the next byte.
 */
static size_t p4nenc_bound(size_t n, size_t size) {
	// return ((n + 127) / 128 + (n + 32) * sizeof(uint64_t));
	return (n + 127) / 128 + (n + 32) * size;
}

// ----- Vertical Minimizer Section -----

/* Vertical Minimizer Section is a section that contains the minimizers of a sequence.
 * The minimizers are stored in a compacted form, with the following columns:
 * - n value: the number of super k-mers in the column
 * - m index: the index of the minimizer in the column
 * - data: the data associated with the minimizer
 * - seq: the sequence associated with the minimizer
 */
void Section_Minimizer::read_section_header() {
	// 1. Verify section type
	char type;
	this->file->read((uint8_t*)&type, 1);
	if (type != 'M')
		throw std::runtime_error("The section do not start with the 'M' char, you can't open a vertical minimizer sequence section.");

	// 2. Read the minimizer
	this->file->read(this->minimizer, this->nb_bytes_mini);

	// 3. Read the number of following super k-mers
	uint8_t buff[8];
	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->nb_blocks);
	this->remaining_blocks = this->nb_blocks;

	// 4. Read offsets of columns
	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->n_col_offset);
	this->n_col_offset += this->start_pos;

	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->m_idx_col_offset);
	this->m_idx_col_offset += this->start_pos;

	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->data_col_offset);
	this->data_col_offset += this->start_pos;

	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->seq_col_offset);
	this->seq_col_offset += this->start_pos;
}


/* Write the section header for a vertical minimizer section.
 * The header contains:
 * - Section type
 * - Minimizer
 * - Number of super k-mers
 * - Column offsets (placeholders for now)
 */
void Section_Minimizer::write_section_header() {
	// 1. Write Section type
	char type = 'M';
	this->file->write(reinterpret_cast<uint8_t *>(&type), 1);

	// 2. Write minimizer
	this->file->write(this->minimizer, this->nb_bytes_mini);

	// 3. Write super k-mer count
	uint8_t buff[8];
	store_big_endian(buff, 8, this->nb_blocks);
	this->file->write(buff, 8);

	// 4. Write column offset placeholder
	memset(buff, 0, 8);
	this->n_col_offset = this->file->tellp();
	this->file->write(buff, 8); // n_col_offset placeholder
	this->m_idx_col_offset = this->file->tellp();
	this->file->write(buff, 8); // m_idx_col_offset placeholder
	this->data_col_offset = this->file->tellp();
	this->file->write(buff, 8); // data_col_offset placeholder
	this->seq_col_offset = this->file->tellp();
	this->file->write(buff, 8); // seq_col_offset placeholder

	// Note: at this time, these four offsets point to the locations of themselves,
	// not the start of the actual columns. They will be updated in the close() function.
}


/* Write the columns of the vertical minimizer section.
 * The columns are written in the following order:
 * 1. n value column
 * 2. m index column
 * 3. data column
 * 4. seq column
 */
void Section_Minimizer::write_columns() {
	uint8_t buff[8];

	// Pre-allocate buffers for compression
	size_t compressed_buf_size = std::max(p4nenc_bound(n_value_buffer.size(), sizeof(uint64_t)),
		std::max(p4nenc_bound(m_idx_buffer.size(), sizeof(uint64_t)), p4nenc_bound(data_buffer.size(), sizeof(uint8_t))));
	auto* compressed_buf = new uint8_t[compressed_buf_size];

	// 1. Write n value column
	this->n_col_offset = this->file->tellp();
	{
		// Compress n_value_buffer
		uint64_t compressed_n_size = p4nenc64(n_value_buffer.data(), n_value_buffer.size(), compressed_buf);
		// Write the size of the compressed data
		store_big_endian(buff, 8, compressed_n_size);
		this->file->write(buff, 8);
		// Write the compressed data
		this->file->write(compressed_buf, compressed_n_size);
	}

	// 2. Write m_idx column
	{
		this->m_idx_col_offset = this->file->tellp();
		// Compress m_idx_buffer
		uint64_t compressed_m_idx_size = p4nenc64(m_idx_buffer.data(), m_idx_buffer.size(), compressed_buf);
		// Write the size of the compressed data
		store_big_endian(buff, 8, compressed_m_idx_size);
		this->file->write(buff, 8);
		// Write the compressed data
		this->file->write(compressed_buf, compressed_m_idx_size);
	}

	// 3. Write data column
	{
		this->data_col_offset = this->file->tellp();
		// Write the size of the data
		store_big_endian(buff, 8, this->data_buffer.size());
		this->file->write(buff, 8);
		// Compress data_buffer
		std::vector<unsigned char> compressed_data_buf(p4nenc_bound(data_buffer.size(), sizeof(uint8_t)));
		uint64_t compressed_data_size = p4nenc8(data_buffer.data(), data_buffer.size(), compressed_data_buf.data());
		// Write the size of the compressed data
		store_big_endian(buff, 8, compressed_data_size);
		this->file->write(buff, 8);
		// Write the compressed data
		this->file->write(compressed_data_buf.data(), compressed_data_size);
	}

	// 4. Write seq column
	this->seq_col_offset = this->file->tellp();
	this->file->write(this->seq_buffer.data(), this->seq_buffer.size());
}


/* Backfill the column offsets in the section header.
 * This function updates the offsets of the columns to point to their actual data locations.
 * It is called at the end of the section writing process.
 */
void Section_Minimizer::backfill_column_offsets() {
	// Save the original position
	uint64_t original_pos = this->file->tellp();

	// Backfill the offsets
	// Save the relative position to the start of this section
	uint8_t buff[8];
	uint64_t n_col_offset_idx = this->start_pos + 1 + this->nb_bytes_mini + 8;
	store_big_endian(buff, 8, n_col_offset - start_pos);
	this->file->write_at(buff, 8, n_col_offset_idx);

	store_big_endian(buff, 8, m_idx_col_offset - start_pos);
	this->file->write_at(buff, 8, n_col_offset_idx + 8);

	store_big_endian(buff, 8, data_col_offset - start_pos);
	this->file->write_at(buff, 8, n_col_offset_idx + 16);

	store_big_endian(buff, 8, seq_col_offset - start_pos);
	this->file->write_at(buff, 8, n_col_offset_idx + 24);

	// Return to the original position
	this->file->jump_to(original_pos);
}


/* Section_Minimizer constructor
 * Initializes the section with the file and reads the header if necessary.
 * Throws an exception if required global variables are missing.
 */
Section_Minimizer::Section_Minimizer(Kero_file* file) : Section(file) {
	this->start_pos = file->tellp();

	if (file->global_vars.find("k") == file->global_vars.end())
		throw "Impossible to read the minimizer section due to missing k variable";
	if (file->global_vars.find("m") == file->global_vars.end())
		throw "Impossible to read the minimizer section due to missing m variable";
	if(file->global_vars.find("max") == file->global_vars.end())
		throw "Impossible to read the minimizer section due to missing max variable";
	if(file->global_vars.find("data_size") == file->global_vars.end())
		throw "Impossible to read the minimizer section due to missing data_size variable";

	uint64_t k = file->global_vars["k"];
	uint64_t m = file->global_vars["m"];
	uint64_t max = file->global_vars["max"];
	uint64_t data_size = file->global_vars["data_size"];

	this->nb_blocks = 0;
	this->remaining_blocks = 0;

	this->k = k;
	this->m = m;
	this->max = max;
	this->data_size = data_size;

	this->cur_skmer_idx = 0;
	this->last_n_pos = 0;
	this->last_m_idx_pos = 0;
	this->last_seq_pos = 0;
	this->last_data_pos = 0;

	this->nb_bytes_mini = 0;
	this->mini_pos_bytes = 0;
	this->nb_blocks = 0;

	this->n_col_offset = 0;
	this->m_idx_col_offset = 0;
	this->seq_col_offset = 0;
	this->data_col_offset = 0;

	// Computes the number of bytes needed to store the number of kmers in each block
	auto nb_bits = static_cast<uint64_t>(ceil(log2(max)));
	this->nb_kmers_bytes = static_cast<uint8_t>(bytes_from_bit_array(nb_bits, 1));
	this->nb_bytes_mini = static_cast<uint8_t>(bytes_from_bit_array(2, m));
	this->minimizer = new uint8_t[nb_bytes_mini];
	memset(this->minimizer, 0, nb_bytes_mini);
	uint64_t mini_pos_bits = static_cast<uint8_t>(ceil(log2(k+max-1)));
	this->mini_pos_bytes = bytes_from_bit_array(mini_pos_bits, 1);

	if (file->is_reader) {
		this->read_section_header();
	}
}


/* Move constructor for Section_Minimizer
 * Transfers ownership of the file and other members from the source object.
 */
Section_Minimizer& Section_Minimizer::operator= (Section_Minimizer && smv) {
	file = smv.file;
	smv.file = nullptr;
	beginning = smv.beginning;
	nb_blocks = smv.nb_blocks;
    nb_blocks = smv.nb_blocks;

	m = smv.m;
	k = smv.k;
	max = smv.max;
	data_size = smv.data_size;

	this->remaining_blocks = smv.remaining_blocks;
	this->nb_kmers_bytes = nb_kmers_bytes;

	nb_bytes_mini = smv.nb_bytes_mini;

	n_col_offset = smv.n_col_offset;
	m_idx_col_offset = smv.m_idx_col_offset;
	seq_col_offset = smv.seq_col_offset;
	data_col_offset = smv.data_col_offset;

	std::swap(minimizer, smv.minimizer);

	return *this;
}


/* * Destructor for Section_Minimizer
 * Cleans up the minimizer array if it was allocated.
 */
Section_Minimizer::~Section_Minimizer() {
	if (minimizer != nullptr) {  // 添加保护
		delete[] minimizer;
		minimizer = nullptr;
	}
}


/* Write the minimizer to the internal minimizer variable.
 * This function does not write to the file directly, but stores the minimizer
 * in the internal variable for later writing in the close() function.
 */
void Section_Minimizer::write_minimizer(uint8_t* minimizer) {
	// Don't directly write the minimizer. Copy it to the internal minimizer variable instead
	// Writing will be done in the close() function.
	memcpy(this->minimizer, minimizer, this->nb_bytes_mini);
}


/* Write a compacted sequence without the minimizer.
 * This function writes the number of k-mers, the minimizer index, the data array,
 * and the compacted sequence into their respective buffers.
 * It is used when the minimizer is not included in the sequence.
 */
void Section_Minimizer::write_compacted_sequence_without_mini(
	uint8_t* seq, uint64_t seq_size, uint64_t mini_pos, uint8_t* data_array) {
	// 1. Calculate the number of k-mers in the current super k-mer
	uint64_t nb_kmers = seq_size + this->m - this->k + 1;

	// 2. Write the number of k-mers and the minimizer index into memory buffers column-wise
	n_value_buffer.push_back(nb_kmers);
	m_idx_buffer.push_back(mini_pos);

	// 3. Process the k-mer counts
	size_t data_bytes = this->data_size * nb_kmers;
	data_buffer.insert(data_buffer.end(), data_array, data_array + data_bytes);

	// 4. Process the compacted sequence without minimizer
	size_t seq_bytes = bytes_from_bit_array(2, seq_size);
	seq_buffer.insert(seq_buffer.end(), seq, seq + seq_bytes);

	// 5. Update the number of super k-mers
	this->nb_blocks++;
}


/* Write a compacted sequence with the minimizer included.
 * This function compacts the sequence by moving the suffix to the position of the minimizer,
 * shifting it to align with the minimizer, and writing the compacted sequence into the file.
 */
void Section_Minimizer::write_compacted_sequence(
	uint8_t* seq, uint64_t seq_size, uint64_t mini_pos, uint8_t* data_array) {
	// Compute all the bit and byte quantities needed.
	uint64_t seq_bytes = bytes_from_bit_array(2, seq_size);
	uint left_offset_nucl = (4 - seq_size % 4) % 4;

	// 1 - Prepare space for sequence manipulation
	uint8_t * seq_copy = new uint8_t[seq_bytes];
	memcpy(seq_copy, seq, seq_bytes);

	// 2 - Move the suffix to the bytes where the minimiser started
	uint mini_start_byte = (mini_pos + left_offset_nucl) / 4;
	uint suff_start_byte = (mini_pos + m + left_offset_nucl) / 4;
	uint suff_bytes = seq_bytes-suff_start_byte;
	for (uint i=0 ; i<suff_bytes ; i++) {
		seq_copy[mini_start_byte + i] = seq_copy[suff_start_byte + i];
	}

	// 3 - shift the suffix to align with the minimizer position
	uint mini_offset = (mini_pos + left_offset_nucl) % 4;
	uint suff_offset = (mini_pos + m + left_offset_nucl) % 4;
	if (mini_offset < suff_offset)
		leftshift8(seq_copy + mini_start_byte, seq_bytes - mini_start_byte, (suff_offset - mini_offset) * 2);
	else
		rightshift8(seq_copy + mini_start_byte, seq_bytes - mini_start_byte, (mini_offset - suff_offset) * 2);

	// 4 - fusion the common byte
	seq_copy[mini_start_byte] = fusion8(seq[mini_start_byte], seq_copy[mini_start_byte], mini_offset * 2);

	// 5 - put all the offset bits on the left
	leftshift8(seq_copy, seq_bytes, left_offset_nucl * 2);
	rightshift8(seq_copy, seq_bytes, ((4 - ((seq_size - m) % 4)) % 4) * 2);

	// Write the compacted sequence into file.
	this->write_compacted_sequence_without_mini(seq_copy, seq_size-m, mini_pos, data_array);

	delete[] seq_copy;
}


/* Add a minimizer to the sequence.
 * This function takes a sequence, the number of k-mers, and the position of the minimizer,
 * and adds the minimizer to the sequence at the specified position.
 * It shifts the suffix to align with the minimizer and merges them into the sequence.
 */
void Section_Minimizer::add_minimizer(uint64_t nb_kmer, uint8_t* seq, uint64_t mini_pos) {
	uint64_t seq_size = nb_kmer + k - 1;
	uint64_t seq_bytes = bytes_from_bit_array(2, seq_size);
	uint64_t seq_left_offset = (4 - (seq_size % 4)) % 4;
	uint64_t no_mini_size = seq_size - m;
	uint64_t no_mini_bytes = bytes_from_bit_array(2, no_mini_size);
	uint64_t no_mini_left_offset = (4 - (no_mini_size % 4)) % 4;
	// Shift the whole sequence to the left to have no padding on byte 0.
	leftshift8(seq, no_mini_bytes, no_mini_left_offset*2);


	// Prepare the suffix
	uint8_t * suffix = new uint8_t[seq_bytes];
	memset(suffix, 0, seq_bytes);
	uint suff_nucl = seq_size - m - mini_pos;
	// Values inside seq before any change
	uint no_mini_suff_start_nucl = mini_pos;
	uint no_mini_suff_start_byte = no_mini_suff_start_nucl / 4;
	uint no_mini_suff_bytes = no_mini_bytes - no_mini_suff_start_byte;
	memcpy(suffix, seq + no_mini_suff_start_byte, no_mini_suff_bytes);
	// Shift to the left
	uint no_mini_suff_offset = no_mini_suff_start_nucl % 4;
	leftshift8(suffix, no_mini_suff_bytes, no_mini_suff_offset * 2);


	// Prepare the minimizer
	uint8_t * mini = new uint8_t[seq_bytes];
	memset(mini, 0, seq_bytes);
	memcpy(mini, this->minimizer, nb_bytes_mini);
	// Shift to the left
	uint mini_offset = (4 - (m % 4)) % 4;
	leftshift8(mini, nb_bytes_mini, mini_offset * 2);


	// Align the minimizer
	uint final_mini_start_nucl = mini_pos;
	uint final_mini_start_byte = mini_pos / 4;
	uint final_mini_offset = final_mini_start_nucl % 4;
	uint final_mini_byte_size = (m + final_mini_offset + 3) / 4;
	rightshift8(mini, seq_bytes, final_mini_offset * 2);

	// Merge minimizer
	seq[final_mini_start_byte] = fusion8(
		seq[final_mini_start_byte],
		mini[0],
		final_mini_offset * 2
	);
	for (uint idx=1 ; idx<final_mini_byte_size ; idx++) {
		seq[final_mini_start_byte+idx] = mini[idx];
	}


	// Align the suffix with the end of the minimizer
	uint final_suff_start_nucl = final_mini_start_nucl + m;
	uint final_suff_start_byte = final_suff_start_nucl / 4;
	uint final_suff_offset = final_suff_start_nucl % 4;
	uint final_suff_byte_size = (suff_nucl + final_suff_offset + 3) / 4;
	if (final_suff_byte_size > 0)
	{
		rightshift8(suffix, seq_bytes, final_suff_offset * 2);

		// Merge the suffix
		seq[final_suff_start_byte] = fusion8(
			seq[final_suff_start_byte],
			suffix[0],
			final_suff_offset * 2
		);
		for (uint64_t idx=1 ; idx<final_suff_byte_size ; idx++) {
			seq[final_suff_start_byte+idx] = suffix[idx];
		}
	}

	// Align everything to the right
	rightshift8(seq, seq_bytes, seq_left_offset * 2);

	delete[] suffix;
	delete[] mini;
}


/* Read a compacted sequence without the minimizer.
 * This function reads the sequence and data from the file, and returns the number of k-mers in the sequence.
 * It also updates the position of the minimizer in the sequence.
 */
uint64_t Section_Minimizer::read_compacted_sequence(uint8_t* seq, uint8_t* data) {
	// Read the block
	uint64_t mini_pos;
	uint64_t nb_kmers_in_skmer = this->read_compacted_sequence_without_mini(seq, data, mini_pos);
	this->add_minimizer(nb_kmers_in_skmer, seq, mini_pos);

	return nb_kmers_in_skmer;
}


/* Read a compacted sequence with the minimizer.
 * This function reads the sequence and data from the file, and returns the number of k-mers in the sequence.
 * It also updates the position of the minimizer in the sequence.
 * The sequence is returned in a single buffer, with the sequence followed by the data.
 */
uint64_t Section_Minimizer::read_compacted_sequence(uint8_t* seq_data) {
	// Read the block
	uint64_t mini_pos;
	uint8_t* seq = new uint8_t[bytes_from_bit_array(2, this->k + this->max - 1)];
	uint8_t* data = new uint8_t[this->max * this->data_size];
	uint64_t nb_kmers_in_skmer = this->read_compacted_sequence_without_mini(seq, data, mini_pos);

	// Concatenate the sequence and the data
	uint64_t seq_size = bytes_from_bit_array(2, nb_kmers_in_skmer + this->k - this->m - 1);
	memcpy(seq_data, seq, seq_size);
	memcpy(seq_data + seq_size, data, this->data_size * nb_kmers_in_skmer);

	delete[] seq;
	delete[] data;

	// Determine the number of new bytes needed for minimizer insertion
	uint64_t seq_size_nucls = k - m + nb_kmers_in_skmer - 1;
	uint64_t free_nucls = (4 - seq_size_nucls) % 4;
	uint64_t bytes_needed = (m - free_nucls + 3) / 4;

	if (bytes_needed > 0)
	{
		// Move the data to the right to have the hole needed for minimizer insertion
		uint64_t seq_size_bytes = (seq_size_nucls + 3) / 4;
		uint64_t data_size_bytes = data_size * nb_kmers_in_skmer;
		for (size_t i=0 ; i<data_size_bytes ; i++)
		{
			// Shift the ith byte to "byte_needed" bytes on the right
			size_t byte_idx = seq_size_bytes + data_size_bytes - 1 - i;
			seq_data[byte_idx + bytes_needed] = seq_data[byte_idx];
			seq_data[byte_idx] = 0;
		}
	}

	// Insert minimizer
	this->add_minimizer(nb_kmers_in_skmer, seq_data, mini_pos);

	return nb_kmers_in_skmer;
}


/* Read a compacted sequence without the minimizer.
 * This function reads the sequence and data from the file, and returns the number of k-mers in the sequence.
 * It also updates the position of the minimizer in the sequence.
 * The sequence is returned in a single buffer, with the sequence followed by the data.
 */
uint64_t Section_Minimizer::read_compacted_sequence_without_mini(
	uint8_t *seq, uint8_t *data, uint64_t &mini_pos) {
	if (this->cur_skmer_idx >= this->nb_blocks) return 0;

	uint8_t buff[8];

	// Initialize the related positions
	if (this->cur_skmer_idx == 0) {
		this->last_n_pos = 0;
		this->last_m_idx_pos = 0;
		if (this->data_size > 0) {
			// this->last_data_pos = data_col_offset;
			this->last_data_pos = 0;
		}

		// last_seq_pos is special because seq is read from the file, not from the buffer in memory
		this->last_seq_pos = seq_col_offset;

		// Uncompress the n_value column
		this->file->jump_to(this->n_col_offset);
		this->file->read(buff, 8);
		uint64_t compressed_n_size;
		load_big_endian(buff, 8, compressed_n_size);
		{
			std::vector<uint8_t> compressed_n_buf(compressed_n_size);
			this->file->read(compressed_n_buf.data(), compressed_n_size);
			// Align the buffer to 8 bytes for p4ndec64
			while (compressed_n_buf.size() < 8) {
				compressed_n_buf.push_back(0);
			}
			this->n_value_buffer.resize(this->nb_blocks);
			p4ndec64(compressed_n_buf.data(), this->nb_blocks, this->n_value_buffer.data());
		}

		// Uncompress the m_idx column
		this->file->jump_to(this->m_idx_col_offset);
		this->file->read(buff, 8);
		uint64_t compressed_m_idx_size;
		load_big_endian(buff, 8, compressed_m_idx_size);
		{
			std::vector<uint8_t> compressed_m_idx_buf(compressed_m_idx_size);
			this->file->read(compressed_m_idx_buf.data(), compressed_m_idx_size);
			// Align the buffer to 8 bytes for p4ndec64
			while (compressed_m_idx_buf.size() < 8) {
				compressed_m_idx_buf.push_back(0);
			}
			this->m_idx_buffer.resize(this->nb_blocks);
			p4ndec64(compressed_m_idx_buf.data(), this->nb_blocks, this->m_idx_buffer.data());
		}

		// Uncompress the data column
		if (this->data_size > 0) {
			this->file->jump_to(this->data_col_offset);
			// Read the size of the data
			this->file->read(buff, 8);
			uint64_t nb_data_buf;
			load_big_endian(buff, 8, nb_data_buf);
			// Read the size of the compressed data
			this->file->read(buff, 8);
			uint64_t compressed_data_size;
			load_big_endian(buff, 8, compressed_data_size);
			{
				std::vector<uint8_t> compressed_data_buf(compressed_data_size);
				this->file->read(compressed_data_buf.data(), compressed_data_size);
				// Align the buffer to 8 bytes for p4ndec64
				while (compressed_data_buf.size() < 8) {
					compressed_data_buf.push_back(0);
				}
				this->data_buffer.resize(nb_data_buf);
				p4ndec8(compressed_data_buf.data(), nb_data_buf, this->data_buffer.data());
			}
		}
	}

	// Read n
	uint64_t n = this->n_value_buffer[this->last_n_pos++];

	// Read m_idx
	mini_pos = this->m_idx_buffer[this->last_m_idx_pos++];

	// Read data
	if (data != nullptr && this->data_size > 0) {
		uint64_t nb_data_bytes = this->data_size * n;
		for (int i = 0; i < nb_data_bytes; i++) {
			data[i] = this->data_buffer[this->last_data_pos++];
		}
	}

	// Read seq
	uint64_t nb_seq_bytes = bytes_from_bit_array(2, n + this->k - this->m - 1);
	this->file->jump_to(this->last_seq_pos);
	this->file->read(seq, nb_seq_bytes);
	this->last_seq_pos += nb_seq_bytes;

	this->cur_skmer_idx++;
	this->remaining_blocks--;
	return n;
}


/* Copy the current Section_Minimizer to another Kero_file.
 * This function creates a new Section_Minimizer in the provided file and copies the minimizer and sequences.
 * It does not write the minimizer directly, but stores it in the new section.
 */
void Section_Minimizer::copy(Kero_file * file) {
	uint max_nucl = this->k + this->max - 1;
	uint8_t * tmp_seq_buffer = new uint8_t[(max_nucl + 3) / 4];
	uint8_t * tmp_data_buffer = new uint8_t[this->max * this->data_size];
	uint64_t mini_pos = 0;

	// Open the copy
	// Note: file is in write mode
	Section_Minimizer smv(file);

	memcpy(smv.minimizer, this->minimizer, this->nb_bytes_mini);

	for (uint64_t i = 0; i < this->nb_blocks; i++) {
        uint64_t nb_kmers = this->read_compacted_sequence_without_mini(tmp_seq_buffer, tmp_data_buffer, mini_pos);
        smv.write_compacted_sequence_without_mini(tmp_seq_buffer, nb_kmers, mini_pos, tmp_data_buffer);
    }

	// Close
	smv.close();
	delete[] tmp_seq_buffer;
	delete[] tmp_data_buffer;
}


/* Jump to the next sequence in the minimizer section.
 * This function is used when reading the section in a reader mode.
 * It skips the current sequence and prepares for the next one.
 * Note: This function does not actually read any data, it just jumps to the next sequence.
 */
void Section_Minimizer::jump_sequence() {
    // These variables are not used, just as placeholders
    uint64_t seq_size = this->k + this->max - 1;
    auto* seq = new uint8_t[bytes_from_bit_array(2, seq_size)];
    auto* data = new uint8_t[this->max * this->data_size];
    uint64_t mini_pos = 0;
    this->read_compacted_sequence_without_mini(seq, data, mini_pos);
    delete[] seq;
    delete[] data;
}


/* Close the Section_Minimizer.
 * This function writes the section header, writes the columns, and backfills the column offsets.
 * It also cleans up the minimizer array and closes the section.
 */
void Section_Minimizer::close() {
	if (this->file->is_writer) {
		// 1. Register the position in the hashtable section
		if (this->file->indexed) {
			this->file->register_minimizer_section(mask_mini(this->minimizer, this->m));
		}

		// 2. Write the section header
		this->write_section_header();

		// 3. Write the columns
		this->write_columns();

		// 4. Backfill the column offsets
		this->backfill_column_offsets();
	}

	if (this->file->is_reader) {
		while (this->remaining_blocks > 0)
			this->jump_sequence();
	}

	// Do cleanup
	delete[] this->minimizer;
	this->minimizer = nullptr;
	Section::close();
}

// Helper function to read data from a mmap pointer instead of a file stream
static void mmap_read(const uint8_t* mmap_ptr, uint64_t offset, uint8_t* dest, uint64_t size) {
    memcpy(dest, mmap_ptr + offset, size);
}

/* Precache columns from a memory-mapped file.
 * This function reads the compressed columns from the mmap pointer and decompresses them into internal buffers.
 * It is used to avoid repeated file I/O when accessing the columns multiple times.
 */
void Section_Minimizer::precache_columns_from_mmap(const uint8_t* mmap_ptr) {
    if (!n_value_buffer.empty()) return; // Already cached

    uint8_t buff[8];

    // Uncompress the n_value column
    uint64_t compressed_n_size;
    mmap_read(mmap_ptr, this->n_col_offset, buff, 8); // Read compressed size
    load_big_endian(buff, 8, compressed_n_size);
    {
        std::vector<uint8_t> compressed_n_buf(compressed_n_size);
        mmap_read(mmap_ptr, this->n_col_offset + 8, compressed_n_buf.data(), compressed_n_size);

        this->n_value_buffer.resize(this->nb_blocks);
        if (compressed_n_size > 0) {
           p4ndec64(compressed_n_buf.data(), this->nb_blocks, this->n_value_buffer.data());
        }
    }

    // Uncompress the m_idx column
    uint64_t compressed_m_idx_size;
    mmap_read(mmap_ptr, this->m_idx_col_offset, buff, 8); // Read compressed size
    load_big_endian(buff, 8, compressed_m_idx_size);
    {
        std::vector<uint8_t> compressed_m_idx_buf(compressed_m_idx_size);
        mmap_read(mmap_ptr, this->m_idx_col_offset + 8, compressed_m_idx_buf.data(), compressed_m_idx_size);

        this->m_idx_buffer.resize(this->nb_blocks);
        if (compressed_m_idx_size > 0) {
            p4ndec64(compressed_m_idx_buf.data(), this->nb_blocks, this->m_idx_buffer.data());
        }
    }

    // Uncompress the data column
    if (this->data_size > 0) {
        uint64_t nb_data_buf;
        mmap_read(mmap_ptr, this->data_col_offset, buff, 8); // Read total data buffer size
        load_big_endian(buff, 8, nb_data_buf);

        uint64_t compressed_data_size;
        mmap_read(mmap_ptr, this->data_col_offset + 8, buff, 8); // Read compressed data size
        load_big_endian(buff, 8, compressed_data_size);

        if (compressed_data_size > 0) {
            std::vector<uint8_t> compressed_data_buf(compressed_data_size);
            mmap_read(mmap_ptr, this->data_col_offset + 16, compressed_data_buf.data(), compressed_data_size);

            this->data_buffer.resize(nb_data_buf);
            p4ndec8(compressed_data_buf.data(), nb_data_buf, this->data_buffer.data());
        }
    }
}

// ----- Hash Table Section -----

/* Section_Hashtable constructor
 * Initializes the section with the file and reads the header if necessary.
 * It reads the section type, the length of the mphf, and the hashtable.
 * Throws an exception if the section type is not 'h'.
 */
Section_Hashtable::Section_Hashtable(Kero_file *file) : Section(file) {
    char type;
    uint8_t buff[8];

//    this->next_index = 0;

    if (this->file->is_reader) {
        // Read the section type
        this->file->read((uint8_t *)&type, 1);
        if (type != 'h')
            throw "The section do not start with the 'h' char, you can not open a Hashtable section.";

        // Read the length of the mphf
        this->file->read(buff, 8);
        load_big_endian(buff, 8, this->nb_mphf);

        // Read the mphf part and generate a temporary file
        std::ofstream temp_file("mphf.bin", std::ios::binary);
        if (!temp_file.is_open())
            throw "Impossible to open the temporary file mphf.bin.";
    	std::vector<uint8_t> buff_chunk(BUFF_CHUNK_SIZE);
        uint64_t nb_bytes_read = 0;
        while (nb_bytes_read < nb_mphf) {
            uint64_t nb_bytes_to_read = std::min(nb_mphf - nb_bytes_read, BUFF_CHUNK_SIZE);
            this->file->read(buff_chunk.data(), nb_bytes_to_read);
            temp_file.write(reinterpret_cast<char*>(buff_chunk.data()), nb_bytes_to_read);
            nb_bytes_read += nb_bytes_to_read;
        }
        temp_file.close();
        mpht.load("mphf.bin");
        std::remove("mphf.bin");

        // Read the length of the hashtable
        uint64_t nb_hashtable;
        this->file->read(buff, 8);
        load_big_endian(buff, 8, nb_hashtable);
    	this->mpht.hashtable.resize(nb_hashtable);

        // Read the hashtable
        for (uint64_t i = 0 ; i < nb_hashtable ; i++) {
            uint64_t key;
            this->file->read(buff, 8);
            load_big_endian(buff, 8, key);
            this->mpht.hashtable[i] = key;
        }
    }
}

Section_Hashtable::~Section_Hashtable() = default;

/* Register a minimizer and its index in the hashtable.
 * This function adds a minimizer and its corresponding index to the internal vectors.
 * It is used when writing the hashtable section to store the minimizers and their positions.
 */
void Section_Hashtable::reg_sm(uint64_t minimizer, uint64_t index) {
    minimizers.push_back(minimizer);
    positions.push_back(index);
}

/* Close the Section_Hashtable.
 * This function builds the minimal perfect hash function (MPHF) from the stored minimizers and positions,
 * writes the section type, the length of the MPHF, and the hashtable to the file.
 * It also cleans up temporary files used during the process.
 */
void Section_Hashtable::close() {
    if (file->is_writer && !minimizers.empty()) {
        uint8_t buff[8];

        // Build mpht
        assert(minimizers.size() == positions.size());
        this->mpht.build(minimizers, positions);

        // Write the section type
        char type = 'h';
        this->file->register_position('h');
        this->file->write((uint8_t *)&type, 1);

		// Save temporary  mphf.bin
    	this->mpht.save("mpht.bin");
        std::ifstream temp_file("mpht.bin", std::ios::binary);
        if (!temp_file.is_open())
            throw "Impossible to open the temporary file mphf.bin.";

    	// Write the length of the mphf
        temp_file.seekg(0, std::ios::end);
        uint64_t nb_mphf = temp_file.tellg();
        temp_file.seekg(0, std::ios::beg);
        store_big_endian(buff, 8, nb_mphf);
        this->file->write(buff, 8);

        // Write the mphf
        uint8_t buff_chunk[BUFF_CHUNK_SIZE];  // 1MB
        uint64_t nb_bytes_written = 0;
        while (nb_bytes_written < nb_mphf) {
            uint64_t nb_bytes_to_write = std::min(nb_mphf - nb_bytes_written, BUFF_CHUNK_SIZE);
            temp_file.read(reinterpret_cast<char*>(buff_chunk), nb_bytes_to_write);
            this->file->write(buff_chunk, nb_bytes_to_write);
            nb_bytes_written += nb_bytes_to_write;
        }
        temp_file.close();
        std::remove("mphf.bin");

        // Write the length of the hashtable
        store_big_endian(buff, 8, this->mpht.size());
        this->file->write(buff, 8);

        // Write the hashtable
        for (uint64_t i = 0 ; i < this->mpht.size() ; i++) {
            store_big_endian(buff, 8, this->mpht.hashtable[i]);
            this->file->write(buff, 8);
        }
    }

    Section::close();
}


// -------- Start of the high level API -----------

Kero_reader::Kero_reader(std::string filename) {
	// Open the file
	this->file = new Kero_file(filename, "r");

	this->current_seq_data = new uint8_t[1];
	this->current_seq_data[0] = 0;

	// Create fake small datastrucutes waiting for the right values.
	this->current_shifts = new uint8_t*[4];
	this->current_shifts[0] = this->current_seq_data;
	for (uint8_t i=1 ; i<4 ; i++) {
		this->current_shifts[i] = new uint8_t[1];
		this->current_shifts[i][0] = 0;
	}

	this->current_section = NULL;
	this->current_kmer = new uint8_t[1];
	this->remaining_kmers = 0;

	this->k = 0;
	this->max = 0;
	this->data_size = 0;

	this->has_next();
}

Kero_reader::~Kero_reader() {
	delete[] this->current_kmer;
	delete[] this->current_seq_data;
	if (this->current_section != NULL)
		delete this->current_section;

	for (uint i=1 ; i<4 ; i++)
		delete[] this->current_shifts[i];
	delete[] this->current_shifts;

	delete this->file;
}

void Kero_reader::read_until_first_section_block() {
	while (current_section == NULL or remaining_blocks == 0) {
		if (this->file->tellp() == this->file->end_position) {
			break;
		}

		// char section_type = this->file->read_section_type();
		char section_type = file->read_section_type();
		// --- Update data structure sizes ---
		if (section_type == 'v')
		{
			// Read the global variable block
			Section_GV gvs(file);
			// Update sequence size if k or max change
			if (gvs.vars.find("k") != gvs.vars.end()
				or gvs.vars.find("max") != gvs.vars.end())
			{
				// Compute the max size of a sequence
				this->k = this->file->global_vars["k"];
				this->max = this->file->global_vars["max"];
				uint64_t seq_max_size = bytes_from_bit_array(2, max + k - 1);
				uint64_t data_max_size = data_size * max;
				// sequence + data buffer
				delete[] this->current_seq_data;
				this->current_seq_data = new uint8_t[seq_max_size + data_max_size];
				
				// Shifts
				this->current_shifts[0] = this->current_seq_data;
				for (uint8_t i=1 ; i<4 ; i++)
				{
					delete[] this->current_shifts[i];
					this->current_shifts[i] = new uint8_t[seq_max_size];
					memset(this->current_shifts[i], 0, seq_max_size);
				}

				// Current kmer
				delete[] this->current_kmer;
				this->current_kmer = new uint8_t[k/4 + 1];
				memset(this->current_kmer, 0, (k/4+1));
			}

			// Update the data array size
			if (gvs.vars.find("data_size") != gvs.vars.end()
				or gvs.vars.find("max") != gvs.vars.end())
			{
				// Update global variables
				this->max = this->file->global_vars["max"];
				this->data_size = this->file->global_vars["data_size"];
				// Update buffer size
				uint64_t seq_max_size = bytes_from_bit_array(2, max + k - 1);
				uint64_t data_max_size = this->data_size * max;
				delete[] this->current_seq_data;
				this->current_seq_data = new uint8_t[seq_max_size + data_max_size];
				this->current_shifts[0] = this->current_seq_data;
				memset(this->current_seq_data, 0, seq_max_size + data_max_size);
			}
		}
		// Mount data from the files to the datastructures.
		else if (section_type == 'i') {
			Section_Index index(file);
			index.close();
		}
        else if (section_type == 'h') {
            Section_Hashtable hashtable(file);
            hashtable.close();
        }
        else {
			current_section = Block_section_reader::construct_section(file);
            if (current_section)
			    remaining_blocks = current_section->nb_blocks;
		}
	}
}


void Kero_reader::read_next_block() {
	// Read from the file
	current_seq_kmers = remaining_kmers = current_section->read_compacted_sequence(current_seq_data);
	current_seq_nucleotides = remaining_kmers + this->k - 1;
	current_seq_bytes = bytes_from_bit_array(2, current_seq_nucleotides);

	// Create the 4 possible shifts of the sequence for easy use.
	for (uint8_t i=1 ; i<min((uint64_t)4, remaining_kmers) ; i++) {
		// Copy
		memcpy(current_shifts[i], current_shifts[0], current_seq_bytes);
		// Shift
		rightshift8(current_shifts[i], current_seq_bytes, 2 * i);
	}
}

bool Kero_reader::has_next() {
	if (current_section == NULL and (file->end_position > file->tellp()))
		read_until_first_section_block();
	return file->end_position > file->tellp();
}

uint64_t Kero_reader::next_block(uint8_t* & sequence, uint8_t* & data) {
	// Verify the abylity to find another kmer in the file.
	if (!this->has_next()){
		sequence = NULL;
		data = NULL;
		return 0;
	}

	uint64_t nb_kmers = current_section->read_compacted_sequence(sequence, data);
	
	remaining_kmers = 0;
	remaining_blocks -= 1;
	if (remaining_blocks == 0) {
		delete current_section;
		current_section = NULL;
	}

	return nb_kmers;
}

bool Kero_reader::next_kmer(uint8_t* & kmer, uint8_t* & data) {
	// Verify the ability to find another kmer in the file.
	if (!this->has_next()){
		kmer = nullptr;
		data = nullptr;
		return false;
	}

	// Load the next block
	if (remaining_kmers == 0) {
		read_next_block();
	}

	uint64_t right_shift = (remaining_kmers - 1) % 4;
	uint64_t prefix_offset = (4 - (current_seq_nucleotides % 4)) % 4;
	uint64_t kmer_idx = current_seq_kmers - remaining_kmers;

	uint64_t start_nucl = prefix_offset + right_shift + kmer_idx;
	uint64_t start_byte = start_nucl / 4;
	uint64_t end_nucl = start_nucl + this->k - 1;
	uint64_t end_byte = end_nucl / 4;

	memcpy(current_kmer, current_shifts[right_shift]+start_byte, end_byte-start_byte+1);
	kmer = current_kmer;
	data = current_seq_data + current_seq_bytes + (current_seq_kmers - remaining_kmers) * this->data_size;
	
	// Read the next block if needed.
	remaining_kmers -= 1;
	if (remaining_kmers == 0) {
		remaining_blocks -= 1;
		if (remaining_blocks == 0) {
			delete current_section;
			current_section = nullptr;
		}
	}

	return true;
}


uint64_t Kero_reader::get_var(string name) {
	if (file->global_vars.find(name) != file->global_vars.end())
		return file->global_vars[name];

	cerr << "Variable " << name << " is absent from the file." << endl;
	exit(2);

	return 0;
}


uint8_t * Kero_reader::get_encoding() {
	return file->encoding;
}
