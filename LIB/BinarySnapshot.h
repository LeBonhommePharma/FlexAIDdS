// BinarySnapshot.h — Compact binary snapshot format for FlexAIDdS GA
// Replaces intermediate text-PDB snapshots with ~10x faster binary I/O.
//
// Apache-2.0 — No GPL dependencies.
// Copyright (c) 2025-2026 Louis-Philippe Morency, Universite de Montreal

#ifndef BINARY_SNAPSHOT_H
#define BINARY_SNAPSHOT_H

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <string>
#include <vector>
#include <stdexcept>

// Platform-specific mmap
#ifdef _WIN32
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#  include <io.h>
#else
#  include <sys/mman.h>
#  include <sys/stat.h>
#  include <fcntl.h>
#  include <unistd.h>
#endif

namespace flexaids {

// ── File format constants ──────────────────────────────────────────
// Layout:
//   [FileHeader]             — 16 bytes
//   [SnapshotRecord] × N     — variable size (record_bytes() each)
//
// FileHeader: magic(4) + version(2) + flags(2) + n_snapshots(4) + n_atoms(4)
// SnapshotRecord: score(f32) + generation(u32) + coords(f32 × 3 × n_atoms)

#pragma pack(push, 1)

/// Binary file header — 16 bytes, little-endian on disk.
struct BinaryFileHeader {
    char     magic[4];       ///< "FASN"
    uint16_t version;        ///< format version (currently 1)
    uint16_t flags;          ///< reserved, must be 0
    uint32_t n_snapshots;    ///< total snapshots in file
    uint32_t n_atoms;        ///< atoms per snapshot (coordinate array length / 3)
};

/// Single snapshot record — header + coordinate payload.
/// On disk: score(4B) + generation(4B) + coords(12*n_atoms B)
struct BinarySnapshotRecord {
    float    score;          ///< energy / CF score (kcal/mol)
    uint32_t generation;     ///< GA generation index
    // Followed by float coords[3 * n_atoms] — NOT packed inline
    // because the array length is variable (depends on header.n_atoms).
};

#pragma pack(pop)

/// Size of a complete on-disk record for a given n_atoms.
inline constexpr size_t snapshot_record_bytes(uint32_t n_atoms) {
    return sizeof(BinarySnapshotRecord) + sizeof(float) * 3 * n_atoms;
}

/// Size of the file header.
inline constexpr size_t file_header_bytes() {
    return sizeof(BinaryFileHeader);
}

// ── Endianness helpers ─────────────────────────────────────────────
// All on-disk values are little-endian.  On big-endian hosts we byte-swap.
// (x86/x86-64/ARM LE are the overwhelmingly common case — hot path is a no-op.)

inline bool host_is_little_endian() {
    const uint16_t test = 0x0001;
    return reinterpret_cast<const uint8_t*>(&test)[0] == 0x01;
}

inline uint16_t le_to_host_u16(uint16_t v) {
    if (host_is_little_endian()) return v;
    return static_cast<uint16_t>((v >> 8) | (v << 8));
}
inline uint32_t le_to_host_u32(uint32_t v) {
    if (host_is_little_endian()) return v;
    return ((v >> 24) & 0x000000FFu) |
           ((v >>  8) & 0x0000FF00u) |
           ((v <<  8) & 0x00FF0000u) |
           ((v << 24) & 0xFF000000u);
}
inline float le_to_host_f32(float v) {
    if (host_is_little_endian()) return v;
    uint32_t tmp;
    std::memcpy(&tmp, &v, sizeof(tmp));
    tmp = le_to_host_u32(tmp);
    float result;
    std::memcpy(&result, &tmp, sizeof(result));
    return result;
}

inline uint16_t host_to_le_u16(uint16_t v) { return le_to_host_u16(v); }
inline uint32_t host_to_le_u32(uint32_t v) { return le_to_host_u32(v); }
inline float    host_to_le_f32(float v)     { return le_to_host_f32(v); }

// ── SnapshotWriter ─────────────────────────────────────────────────
/// Thread-safe writer for binary snapshot files.
/// Opens file on construction, appends records, updates header on flush/close.
class SnapshotWriter {
public:
    /// Open a binary snapshot file for writing.
    /// @param path     File path
    /// @param n_atoms  Number of atoms per snapshot
    /// @throws std::runtime_error on I/O failure
    SnapshotWriter(const std::string& path, uint32_t n_atoms);

    /// Flush and close the file, updating the header with final count.
    ~SnapshotWriter();

    // Non-copyable, movable
    SnapshotWriter(const SnapshotWriter&) = delete;
    SnapshotWriter& operator=(const SnapshotWriter&) = delete;
    SnapshotWriter(SnapshotWriter&&) noexcept;
    SnapshotWriter& operator=(SnapshotWriter&&) noexcept;

    /// Write a single snapshot. Thread-safe.
    /// @param score       Energy / CF score
    /// @param generation  GA generation index
    /// @param coords      Flat coordinate array [x0,y0,z0, x1,y1,z1, ...]
    ///                     Must have exactly 3 * n_atoms elements.
    /// @throws std::runtime_error if coords size mismatches or write fails
    void write_snapshot(float score, uint32_t generation,
                        const std::vector<float>& coords);

    /// Overload for raw pointer (avoids vector copy in hot path).
    void write_snapshot(float score, uint32_t generation,
                        const float* coords, size_t coords_len);

    /// Flush buffered data and update the file header with current count.
    void flush();

    /// Number of snapshots written so far.
    uint32_t count() const { return count_; }

    /// Number of atoms per snapshot.
    uint32_t n_atoms() const { return n_atoms_; }

private:
    std::string path_;
    uint32_t n_atoms_;
    uint32_t count_;
    FILE* file_;
    std::mutex mutex_;

    void write_header();
};

// ── SnapshotReader ─────────────────────────────────────────────────
/// Read-only reader for binary snapshot files.
/// Supports random access by index and optional mmap for zero-copy.
class SnapshotReader {
public:
    /// Result struct for a single snapshot.
    struct Snapshot {
        float    score;
        uint32_t generation;
        std::vector<float> coords;   ///< x0,y0,z0, x1,y1,z1, ...
    };

    /// Open and validate a binary snapshot file.
    /// @throws std::runtime_error if file not found, truncated, or magic mismatch
    explicit SnapshotReader(const std::string& path);

    /// Close file and release mmap if any.
    ~SnapshotReader();

    // Non-copyable, non-movable (owns mmap)
    SnapshotReader(const SnapshotReader&) = delete;
    SnapshotReader& operator=(const SnapshotReader&) = delete;

    /// Read all snapshots into memory.
    std::vector<Snapshot> read_all() const;

    /// Read a single snapshot by index.
    /// @throws std::out_of_range if index >= n_snapshots
    Snapshot read_snapshot(uint32_t index) const;

    /// Number of snapshots in the file.
    uint32_t n_snapshots() const { return header_.n_snapshots; }

    /// Number of atoms per snapshot.
    uint32_t n_atoms() const { return header_.n_atoms; }

    /// Access mmapped coordinate data for zero-copy.
    /// Returns nullptr if mmap failed or not available.
    /// The pointer is valid for the lifetime of this SnapshotReader.
    /// Layout: contiguous array of SnapshotRecord + coords data.
    const float* mmap_coordinates(uint32_t snapshot_index) const;

    /// Check if the file at path starts with the "FASN" magic bytes.
    /// Returns true if it is a valid binary snapshot file.
    static bool is_binary_snapshot(const std::string& path);

private:
    std::string path_;
    BinaryFileHeader header_;
    size_t file_size_;

    // mmap state
#ifdef _WIN32
    HANDLE file_handle_;
    HANDLE mapping_handle_;
#else
    int fd_;
#endif
    void* mapped_data_;
    size_t mapped_size_;

    void mmap_file();
    void munmap_file();
    void validate_header();
};

// ── Conversion to text PDB ────────────────────────────────────────
/// Convert a single binary snapshot to text PDB format.
/// @param snapshot    The snapshot data
/// @param atom_names  Per-atom PDB atom names (e.g. " CA ", " N  ")
/// @param res_names   Per-atom residue names (e.g. "ALA", "GLY")
/// @param res_numbers Per-atom residue numbers
/// @param chains      Per-atom chain IDs
/// @param elements    Per-atom element strings (e.g. " C", " N")
/// @param atom_numbers Per-atom PDB atom serial numbers
/// @param remark      REMARK line to include (may be empty)
/// @param out         FILE* to write to (must be opened by caller)
void snapshot_to_pdb(const SnapshotReader::Snapshot& snapshot,
                     const std::vector<std::string>& atom_names,
                     const std::vector<std::string>& res_names,
                     const std::vector<int>& res_numbers,
                     const std::vector<char>& chains,
                     const std::vector<std::string>& elements,
                     const std::vector<int>& atom_numbers,
                     const std::string& remark,
                     FILE* out);

/// Convenience overload: extract atom metadata from FA_Global/resid arrays
/// and write PDB.  Forward-declared here; implemented in BinarySnapshot.cpp
/// after including flexaid.h.
struct FA_Global_struct;
struct atom_struct;
struct residue_struct;
void snapshot_to_pdb_from_global(const SnapshotReader::Snapshot& snapshot,
                                 FA_Global_struct* FA,
                                 atom_struct* atoms,
                                 residue_struct* residue,
                                 const std::string& remark,
                                 FILE* out);

} // namespace flexaids

#endif // BINARY_SNAPSHOT_H
