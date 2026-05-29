// BinarySnapshot.cpp — Implementation of compact binary snapshot format
// Apache-2.0 — No GPL dependencies.
// Copyright (c) 2025-2026 Louis-Philippe Morency, Universite de Montreal

#include "BinarySnapshot.h"
#include "flexaid.h"
#include "fileio.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstring>
#include <new>

namespace flexaids {

// ═══════════════════════════════════════════════════════════════════
//  SnapshotWriter
// ═══════════════════════════════════════════════════════════════════

SnapshotWriter::SnapshotWriter(const std::string& path, uint32_t n_atoms)
    : path_(path)
    , n_atoms_(n_atoms)
    , count_(0)
    , file_(nullptr)
{
    if (n_atoms == 0) {
        throw std::runtime_error("SnapshotWriter: n_atoms must be > 0");
    }
    file_ = std::fopen(path.c_str(), "wb");
    if (!file_) {
        throw std::runtime_error(
            std::string("SnapshotWriter: failed to open '") + path +
            "' for writing: " + std::strerror(errno));
    }
    // Write a placeholder header — will be overwritten on flush/close.
    write_header();
}

SnapshotWriter::~SnapshotWriter() {
    if (file_) {
        flush();
        std::fclose(file_);
        file_ = nullptr;
    }
}

SnapshotWriter::SnapshotWriter(SnapshotWriter&& other) noexcept
    : path_(std::move(other.path_))
    , n_atoms_(other.n_atoms_)
    , count_(other.count_)
    , file_(other.file_)
{
    other.file_ = nullptr;
    other.count_ = 0;
}

SnapshotWriter& SnapshotWriter::operator=(SnapshotWriter&& other) noexcept {
    if (this != &other) {
        if (file_) {
            flush();
            std::fclose(file_);
        }
        path_ = std::move(other.path_);
        n_atoms_ = other.n_atoms_;
        count_ = other.count_;
        file_ = other.file_;
        other.file_ = nullptr;
        other.count_ = 0;
    }
    return *this;
}

void SnapshotWriter::write_header() {
    BinaryFileHeader hdr{};
    std::memcpy(hdr.magic, "FASN", 4);
    hdr.version     = host_to_le_u16(1);
    hdr.flags       = host_to_le_u16(0);
    hdr.n_snapshots = host_to_le_u32(count_);
    hdr.n_atoms     = host_to_le_u32(n_atoms_);

    std::rewind(file_);
    size_t written = std::fwrite(&hdr, sizeof(hdr), 1, file_);
    if (written != 1) {
        // Header write failure is non-fatal during construction;
        // fatal during flush() — caller checks return.
    }
    std::fflush(file_);
}

void SnapshotWriter::write_snapshot(float score, uint32_t generation,
                                     const std::vector<float>& coords) {
    if (coords.size() != static_cast<size_t>(3) * n_atoms_) {
        throw std::runtime_error(
            "SnapshotWriter::write_snapshot: coords size " +
            std::to_string(coords.size()) + " != expected " +
            std::to_string(3u * n_atoms_));
    }
    write_snapshot(score, generation, coords.data(), coords.size());
}

void SnapshotWriter::write_snapshot(float score, uint32_t generation,
                                     const float* coords, size_t coords_len) {
    if (!coords) {
        throw std::runtime_error("SnapshotWriter::write_snapshot: null coords");
    }
    if (coords_len != static_cast<size_t>(3) * n_atoms_) {
        throw std::runtime_error(
            "SnapshotWriter::write_snapshot: coords length " +
            std::to_string(coords_len) + " != expected " +
            std::to_string(3u * n_atoms_));
    }

    std::lock_guard<std::mutex> lock(mutex_);

    // Write record header
    BinarySnapshotRecord rec;
    rec.score      = host_to_le_f32(score);
    rec.generation = host_to_le_u32(generation);

    size_t w = std::fwrite(&rec, sizeof(rec), 1, file_);
    if (w != 1) {
        throw std::runtime_error(
            "SnapshotWriter: failed to write record header: " +
            std::string(std::strerror(errno)));
    }

    // Write coordinates — on little-endian hosts, no conversion needed
    // for a float32 array.  On big-endian, we must swap each float.
    const size_t n_floats = 3 * static_cast<size_t>(n_atoms_);
    if (host_is_little_endian()) {
        w = std::fwrite(coords, sizeof(float), n_floats, file_);
    } else {
        // Big-endian: swap each float
        std::vector<float> swapped(n_floats);
        for (size_t i = 0; i < n_floats; ++i) {
            swapped[i] = host_to_le_f32(coords[i]);
        }
        w = std::fwrite(swapped.data(), sizeof(float), n_floats, file_);
    }

    if (w != n_floats) {
        throw std::runtime_error(
            "SnapshotWriter: failed to write coordinates: " +
            std::string(std::strerror(errno)));
    }

    ++count_;
}

void SnapshotWriter::flush() {
    if (!file_) return;
    // Seek to beginning, rewrite header with current count.
    write_header();
    std::fflush(file_);
}

// ═══════════════════════════════════════════════════════════════════
//  SnapshotReader
// ═══════════════════════════════════════════════════════════════════

SnapshotReader::SnapshotReader(const std::string& path)
    : path_(path)
    , file_size_(0)
#ifdef _WIN32
    , file_handle_(INVALID_HANDLE_VALUE)
    , mapping_handle_(NULL)
#else
    , fd_(-1)
#endif
    , mapped_data_(nullptr)
    , mapped_size_(0)
{
    // First, check magic bytes
    validate_header();

    // Try mmap for zero-copy access
    mmap_file();
}

SnapshotReader::~SnapshotReader() {
    munmap_file();
}

void SnapshotReader::validate_header() {
#ifdef _WIN32
    HANDLE h = CreateFileA(path_.c_str(), GENERIC_READ, FILE_SHARE_READ,
                           NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (h == INVALID_HANDLE_VALUE) {
        throw std::runtime_error(
            "SnapshotReader: cannot open '" + path_ + "'");
    }
    DWORD bytesRead;
    BinaryFileHeader raw;
    BOOL ok = ReadFile(h, &raw, sizeof(raw), &bytesRead, NULL);
    LARGE_INTEGER fileSize;
    GetFileSizeEx(h, &fileSize);
    CloseHandle(h);

    if (!ok || bytesRead != sizeof(raw)) {
        throw std::runtime_error(
            "SnapshotReader: failed to read header from '" + path_ + "'");
    }
    file_size_ = static_cast<size_t>(fileSize.QuadPart);
#else
    FILE* f = std::fopen(path_.c_str(), "rb");
    if (!f) {
        throw std::runtime_error(
            "SnapshotReader: cannot open '" + path_ +
            "': " + std::strerror(errno));
    }
    BinaryFileHeader raw;
    size_t nread = std::fread(&raw, sizeof(raw), 1, f);

    // Get file size
    std::fseek(f, 0, SEEK_END);
    file_size_ = static_cast<size_t>(std::ftell(f));
    std::fclose(f);

    if (nread != 1) {
        throw std::runtime_error(
            "SnapshotReader: failed to read header from '" + path_ +
            "': file too small");
    }
#endif

    // Validate magic
    if (std::memcmp(raw.magic, "FASN", 4) != 0) {
        throw std::runtime_error(
            "SnapshotReader: invalid magic in '" + path_ +
            "' (expected 'FASN')");
    }

    header_.magic[0] = 'F'; header_.magic[1] = 'A';
    header_.magic[2] = 'S'; header_.magic[3] = 'N';
    header_.version     = le_to_host_u16(raw.version);
    header_.flags       = le_to_host_u16(raw.flags);
    header_.n_snapshots = le_to_host_u32(raw.n_snapshots);
    header_.n_atoms     = le_to_host_u32(raw.n_atoms);

    // Validate version
    if (header_.version != 1) {
        throw std::runtime_error(
            "SnapshotReader: unsupported version " +
            std::to_string(header_.version) + " in '" + path_ + "'");
    }

    // Validate file size against declared content
    size_t expected = file_header_bytes() +
                      static_cast<size_t>(header_.n_snapshots) *
                      snapshot_record_bytes(header_.n_atoms);
    if (file_size_ < expected) {
        throw std::runtime_error(
            "SnapshotReader: file truncated in '" + path_ +
            "' (expected " + std::to_string(expected) + " bytes, got " +
            std::to_string(file_size_) + ")");
    }
}

void SnapshotReader::mmap_file() {
    if (file_size_ == 0) return;

#ifdef _WIN32
    file_handle_ = CreateFileA(path_.c_str(), GENERIC_READ, FILE_SHARE_READ,
                               NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (file_handle_ == INVALID_HANDLE_VALUE) return;

    mapping_handle_ = CreateFileMappingA(file_handle_, NULL, PAGE_READONLY,
                                          0, 0, NULL);
    if (!mapping_handle_) {
        CloseHandle(file_handle_);
        file_handle_ = INVALID_HANDLE_VALUE;
        return;
    }

    mapped_data_ = MapViewOfFile(mapping_handle_, FILE_MAP_READ, 0, 0, 0);
    if (mapped_data_) {
        mapped_size_ = file_size_;
    }
#else
    fd_ = ::open(path_.c_str(), O_RDONLY);
    if (fd_ < 0) return;

    mapped_data_ = ::mmap(nullptr, file_size_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (mapped_data_ == MAP_FAILED) {
        mapped_data_ = nullptr;
        ::close(fd_);
        fd_ = -1;
        return;
    }
    mapped_size_ = file_size_;
#endif
}

void SnapshotReader::munmap_file() {
#ifdef _WIN32
    if (mapped_data_) {
        UnmapViewOfFile(mapped_data_);
        mapped_data_ = nullptr;
    }
    if (mapping_handle_) {
        CloseHandle(mapping_handle_);
        mapping_handle_ = NULL;
    }
    if (file_handle_ != INVALID_HANDLE_VALUE) {
        CloseHandle(file_handle_);
        file_handle_ = INVALID_HANDLE_VALUE;
    }
#else
    if (mapped_data_) {
        ::munmap(mapped_data_, mapped_size_);
        mapped_data_ = nullptr;
    }
    if (fd_ >= 0) {
        ::close(fd_);
        fd_ = -1;
    }
#endif
}

std::vector<SnapshotReader::Snapshot> SnapshotReader::read_all() const {
    std::vector<Snapshot> result;
    result.reserve(header_.n_snapshots);

    const size_t rec_bytes = snapshot_record_bytes(header_.n_atoms);
    const size_t n_floats = 3 * static_cast<size_t>(header_.n_atoms);

    // Use mmap if available, else fall back to stdio.
    if (mapped_data_) {
        const char* base = static_cast<const char*>(mapped_data_) +
                           file_header_bytes();
        for (uint32_t i = 0; i < header_.n_snapshots; ++i) {
            const char* rec_ptr = base + i * rec_bytes;

            BinarySnapshotRecord rec;
            std::memcpy(&rec, rec_ptr, sizeof(rec));
            rec.score      = le_to_host_f32(rec.score);
            rec.generation = le_to_host_u32(rec.generation);

            const float* coords_raw = reinterpret_cast<const float*>(
                rec_ptr + sizeof(BinarySnapshotRecord));

            Snapshot snap;
            snap.score = rec.score;
            snap.generation = rec.generation;
            snap.coords.resize(n_floats);

            if (host_is_little_endian()) {
                std::memcpy(snap.coords.data(), coords_raw,
                            n_floats * sizeof(float));
            } else {
                for (size_t j = 0; j < n_floats; ++j) {
                    snap.coords[j] = le_to_host_f32(coords_raw[j]);
                }
            }
            result.push_back(std::move(snap));
        }
    } else {
        // Fallback: seek-and-read via stdio
        FILE* f = std::fopen(path_.c_str(), "rb");
        if (!f) {
            throw std::runtime_error(
                "SnapshotReader::read_all: cannot reopen '" + path_ + "'");
        }
        std::fseek(f, static_cast<long>(file_header_bytes()), SEEK_SET);

        std::vector<float> tmp_coords(n_floats);

        for (uint32_t i = 0; i < header_.n_snapshots; ++i) {
            BinarySnapshotRecord rec;
            if (std::fread(&rec, sizeof(rec), 1, f) != 1) break;
            rec.score      = le_to_host_f32(rec.score);
            rec.generation = le_to_host_u32(rec.generation);

            if (std::fread(tmp_coords.data(), sizeof(float), n_floats, f)
                != n_floats) break;

            Snapshot snap;
            snap.score = rec.score;
            snap.generation = rec.generation;
            snap.coords.resize(n_floats);

            if (host_is_little_endian()) {
                snap.coords = tmp_coords;
            } else {
                for (size_t j = 0; j < n_floats; ++j) {
                    snap.coords[j] = le_to_host_f32(tmp_coords[j]);
                }
            }
            result.push_back(std::move(snap));
        }
        std::fclose(f);
    }

    return result;
}

SnapshotReader::Snapshot SnapshotReader::read_snapshot(uint32_t index) const {
    if (index >= header_.n_snapshots) {
        throw std::out_of_range(
            "SnapshotReader::read_snapshot: index " +
            std::to_string(index) + " >= n_snapshots " +
            std::to_string(header_.n_snapshots));
    }

    const size_t rec_bytes = snapshot_record_bytes(header_.n_atoms);
    const size_t n_floats  = 3 * static_cast<size_t>(header_.n_atoms);
    const size_t offset    = file_header_bytes() + index * rec_bytes;

    Snapshot snap;

    if (mapped_data_ && offset + rec_bytes <= mapped_size_) {
        const char* rec_ptr = static_cast<const char*>(mapped_data_) + offset;

        BinarySnapshotRecord rec;
        std::memcpy(&rec, rec_ptr, sizeof(rec));
        snap.score      = le_to_host_f32(rec.score);
        snap.generation = le_to_host_u32(rec.generation);

        const float* coords_raw = reinterpret_cast<const float*>(
            rec_ptr + sizeof(BinarySnapshotRecord));

        snap.coords.resize(n_floats);
        if (host_is_little_endian()) {
            std::memcpy(snap.coords.data(), coords_raw,
                        n_floats * sizeof(float));
        } else {
            for (size_t j = 0; j < n_floats; ++j) {
                snap.coords[j] = le_to_host_f32(coords_raw[j]);
            }
        }
    } else {
        FILE* f = std::fopen(path_.c_str(), "rb");
        if (!f) {
            throw std::runtime_error(
                "SnapshotReader::read_snapshot: cannot open '" + path_ + "'");
        }
        std::fseek(f, static_cast<long>(offset), SEEK_SET);

        BinarySnapshotRecord rec;
        if (std::fread(&rec, sizeof(rec), 1, f) != 1) {
            std::fclose(f);
            throw std::runtime_error("SnapshotReader: failed to read record");
        }
        snap.score      = le_to_host_f32(rec.score);
        snap.generation = le_to_host_u32(rec.generation);

        snap.coords.resize(n_floats);
        if (std::fread(snap.coords.data(), sizeof(float), n_floats, f)
            != n_floats) {
            std::fclose(f);
            throw std::runtime_error(
                "SnapshotReader: failed to read coordinates");
        }
        std::fclose(f);

        if (!host_is_little_endian()) {
            for (size_t j = 0; j < n_floats; ++j) {
                snap.coords[j] = le_to_host_f32(snap.coords[j]);
            }
        }
    }

    return snap;
}

const float* SnapshotReader::mmap_coordinates(uint32_t snapshot_index) const {
    if (!mapped_data_) return nullptr;
    if (snapshot_index >= header_.n_snapshots) return nullptr;

    const size_t rec_bytes = snapshot_record_bytes(header_.n_atoms);
    const size_t offset = file_header_bytes() +
                          snapshot_index * rec_bytes +
                          sizeof(BinarySnapshotRecord);

    if (offset + 3 * sizeof(float) * header_.n_atoms > mapped_size_) {
        return nullptr;
    }

    return reinterpret_cast<const float*>(
        static_cast<const char*>(mapped_data_) + offset);
}

bool SnapshotReader::is_binary_snapshot(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return false;

    char magic[4];
    size_t nread = std::fread(magic, 1, 4, f);
    std::fclose(f);

    return (nread == 4 && std::memcmp(magic, "FASN", 4) == 0);
}

// ═══════════════════════════════════════════════════════════════════
//  snapshot_to_pdb — generic version
// ═══════════════════════════════════════════════════════════════════

void snapshot_to_pdb(const SnapshotReader::Snapshot& snapshot,
                     const std::vector<std::string>& atom_names,
                     const std::vector<std::string>& res_names,
                     const std::vector<int>& res_numbers,
                     const std::vector<char>& chains,
                     const std::vector<std::string>& elements,
                     const std::vector<int>& atom_numbers,
                     const std::string& remark,
                     FILE* out)
{
    if (!out) return;

    // REMARK header
    if (!remark.empty()) {
        std::fprintf(out, "REMARK %s\n", remark.c_str());
    }

    // ATOM/HETATM records
    const uint32_t n_atoms = static_cast<uint32_t>(snapshot.coords.size() / 3);
    const uint32_t n_meta = static_cast<uint32_t>(atom_names.size());

    for (uint32_t i = 0; i < n_atoms && i < n_meta; ++i) {
        const char* field = "ATOM  ";  // default; could be HETATM depending on caller
        float x = snapshot.coords[i * 3 + 0];
        float y = snapshot.coords[i * 3 + 1];
        float z = snapshot.coords[i * 3 + 2];
        int anum = (i < atom_numbers.size()) ? atom_numbers[i] : static_cast<int>(i + 1);
        const char* aname = atom_names[i].c_str();
        const char* rname = res_names[i].c_str();
        char chn = (i < chains.size()) ? chains[i] : ' ';
        int rnum = (i < res_numbers.size()) ? res_numbers[i] : 1;
        const char* elem = (i < elements.size()) ? elements[i].c_str() : " C";

        std::fprintf(out,
            "%s%5d %-4s %-4s%c%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  \n",
            field, anum, aname, rname, chn, rnum, x, y, z, elem);
    }

    std::fprintf(out, "END\n");
}

// ═══════════════════════════════════════════════════════════════════
//  snapshot_to_pdb_from_global — FA_Global/resid-aware version
// ═══════════════════════════════════════════════════════════════════

void snapshot_to_pdb_from_global(const SnapshotReader::Snapshot& snapshot,
                                 ::FA_Global_struct* FA,
                                 ::atom_struct* atoms,
                                 ::residue_struct* residue,
                                 const std::string& remark,
                                 FILE* out)
{
    if (!out || !FA || !atoms || !residue) return;

    // REMARK header
    if (!remark.empty()) {
        std::fprintf(out, "REMARK %s\n", remark.c_str());
    }

    // Build per-atom metadata from FA_Global/residue arrays,
    // then overlay snapshot coordinates on top.
    // We iterate residues -> atoms just like write_pdb() does.
    char field[7];

    // We need to apply the snapshot coordinates temporarily.
    // Save original coordinates, overlay snapshot coords, write, then restore.
    // This avoids copying the entire atoms array.
    const uint32_t n_atoms = FA->atm_cnt_real;

    // Save originals
    std::vector<std::array<float, 3>> saved_coords(n_atoms + 1);  // 1-based
    for (int i = 1; i <= static_cast<int>(n_atoms); ++i) {
        saved_coords[i][0] = atoms[i].coor[0];
        saved_coords[i][1] = atoms[i].coor[1];
        saved_coords[i][2] = atoms[i].coor[2];
    }

    // Overlay snapshot coordinates (flat [x0,y0,z0,x1,y1,z1,...] -> atoms)
    // The snapshot stores ALL atoms (protein + ligand) for full-structure output.
    for (uint32_t i = 0; i < snapshot.coords.size() / 3 && i < n_atoms; ++i) {
        atoms[i + 1].coor[0] = snapshot.coords[i * 3 + 0];
        atoms[i + 1].coor[1] = snapshot.coords[i * 3 + 1];
        atoms[i + 1].coor[2] = snapshot.coords[i * 3 + 2];
    }

    // Write using the same loop structure as write_pdb()
    for (int k = 1; k <= FA->res_cnt; ++k) {
        int rot = residue[k].rot;
        for (int i = residue[k].fatm[rot]; i <= residue[k].latm[rot]; ++i) {
            if (!FA->output_scored_only || atoms[i].optres != NULL) {
                if (residue[atoms[i].ofres].type == 0) std::strcpy(field, "ATOM  ");
                if (residue[atoms[i].ofres].type == 1) std::strcpy(field, "HETATM");

                std::fprintf(out,
                    "%s%5d %s %s %c%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  \n",
                    field,
                    atoms[i].number,
                    atoms[i].name,
                    residue[atoms[i].ofres].name,
                    residue[atoms[i].ofres].chn,
                    residue[atoms[i].ofres].number,
                    atoms[i].coor[0],
                    atoms[i].coor[1],
                    atoms[i].coor[2],
                    get_element(atoms[i].type));
            }
        }
    }

    // CONECT records for ligand
    for (int i = 1; i <= FA->num_het; ++i) {
        int rot = residue[FA->het_res[i]].rot;
        for (int j = residue[FA->het_res[i]].fatm[rot];
             j <= residue[FA->het_res[i]].latm[rot]; ++j) {
            std::fprintf(out, "CONECT%5d", atoms[j].number);
            for (int kk = 1; kk <= atoms[j].bond[0]; ++kk) {
                std::fprintf(out, "%5d", atoms[atoms[j].bond[kk]].number);
            }
            std::fprintf(out, "\n");
        }
    }

    std::fprintf(out, "END\n");

    // Restore original coordinates
    for (int i = 1; i <= static_cast<int>(n_atoms); ++i) {
        atoms[i].coor[0] = saved_coords[i][0];
        atoms[i].coor[1] = saved_coords[i][1];
        atoms[i].coor[2] = saved_coords[i][2];
    }
}

} // namespace flexaids
