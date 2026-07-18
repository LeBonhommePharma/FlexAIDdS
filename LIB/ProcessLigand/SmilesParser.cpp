// SmilesParser.cpp — OpenSMILES-compliant parser implementation
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "SmilesParser.h"

#include <cctype>
#include <cassert>
#include <unordered_map>
#include <stdexcept>
#include <limits>
#include <cmath>

namespace bonmol {

// ---------------------------------------------------------------------------
// Organic-subset: these elements may appear unbracketed in SMILES.
// Lowercase versions signal aromaticity.
// Expected valences for implicit-H computation.
// ---------------------------------------------------------------------------
struct OrganicAtomInfo {
    Element  elem;
    int      normal_valence[4]; // ordered list of normal valences, terminated by -1
};

static const std::unordered_map<char, OrganicAtomInfo> kOrganicSubset = {
    {'B',  {Element::B,  {3,-1,-1,-1}}},
    {'C',  {Element::C,  {4,-1,-1,-1}}},
    {'N',  {Element::N,  {3, 5,-1,-1}}},
    {'O',  {Element::O,  {2,-1,-1,-1}}},
    {'P',  {Element::P,  {3, 5,-1,-1}}},
    {'S',  {Element::S,  {2, 4, 6,-1}}},
    {'F',  {Element::F,  {1,-1,-1,-1}}},
    {'I',  {Element::I,  {1,-1,-1,-1}}},
};

// Aromatic organic subset (lowercase), map to same elements
static const std::unordered_map<char, OrganicAtomInfo> kAromaticSubset = {
    {'b',  {Element::B,  {3,-1,-1,-1}}},
    {'c',  {Element::C,  {4,-1,-1,-1}}},
    {'n',  {Element::N,  {3, 5,-1,-1}}},
    {'o',  {Element::O,  {2,-1,-1,-1}}},
    {'p',  {Element::P,  {3, 5,-1,-1}}},
    {'s',  {Element::S,  {2, 4, 6,-1}}},
};

// ---------------------------------------------------------------------------
// SmilesParser implementation
// ---------------------------------------------------------------------------

void SmilesParser::reset() {
    smiles_.clear();
    pos_             = 0;
    mol_             = BonMol{};
    warnings_.clear();
    ring_open_.clear();
    branch_stack_.clear();
    prev_atom_        = -1;
    pending_order_    = std::nullopt;
    pending_aromatic_ = false;
}

char SmilesParser::peek(int offset) const noexcept {
    int idx = pos_ + offset;
    if (idx < 0 || idx >= static_cast<int>(smiles_.size())) return '\0';
    return smiles_[idx];
}

char SmilesParser::consume() noexcept {
    if (pos_ >= static_cast<int>(smiles_.size())) return '\0';
    return smiles_[pos_++];
}

bool SmilesParser::at_end() const noexcept {
    return pos_ >= static_cast<int>(smiles_.size());
}

// ---------------------------------------------------------------------------
// Main parse entry point
// ---------------------------------------------------------------------------

SmilesParseResult SmilesParser::parse(const std::string& smiles) {
    reset();
    // Fail closed on empty / whitespace-only input (tests + pipeline).
    {
        bool any = false;
        for (char c : smiles) {
            if (!std::isspace(static_cast<unsigned char>(c))) { any = true; break; }
        }
        if (!any)
            throw SmilesParseError("empty SMILES string", 0);
    }
    smiles_ = smiles;
    mol_.smiles = smiles;

    while (!at_end()) {
        if (!parse_token()) break;
    }

    // Verify all ring closures were matched
    for (const auto& [rnum, ro] : ring_open_) {
        throw SmilesParseError("unclosed ring closure %" + std::to_string(rnum), pos_);
    }

    // Verify branch stack is empty
    if (!branch_stack_.empty()) {
        throw SmilesParseError("unclosed branch '('", pos_);
    }

    if (mol_.num_atoms() == 0) {
        throw SmilesParseError("no atoms parsed from SMILES", pos_);
    }

    // Compute implicit H for all atoms that don't have explicit brackets
    // (bracket atoms set implicit_h_count during parse_bracket_atom)
    for (int i = 0; i < mol_.num_atoms(); ++i) {
        if (mol_.atoms[i].implicit_h_count < 0) {
            // was set explicitly in bracket — keep as-is (flip sign back)
            mol_.atoms[i].implicit_h_count = -mol_.atoms[i].implicit_h_count - 1;
        } else {
            // organic subset: compute from valence
            mol_.atoms[i].implicit_h_count = compute_implicit_h(i);
        }
    }

    mol_.finalize();

    SmilesParseResult result;
    result.mol      = std::move(mol_);
    result.warnings = std::move(warnings_);
    return result;
}

// ---------------------------------------------------------------------------
// Token dispatch
// ---------------------------------------------------------------------------

bool SmilesParser::parse_token() {
    char c = peek();

    // ---- Bond tokens ----
    if (c == '-') {
        consume();
        pending_order_    = BondOrder::SINGLE;
        pending_aromatic_ = false;
        return true;
    }
    if (c == '=') {
        consume();
        pending_order_    = BondOrder::DOUBLE;
        pending_aromatic_ = false;
        return true;
    }
    if (c == '#') {
        consume();
        pending_order_    = BondOrder::TRIPLE;
        pending_aromatic_ = false;
        return true;
    }
    if (c == ':') {
        consume();
        pending_order_    = BondOrder::AROMATIC;
        pending_aromatic_ = true;
        return true;
    }
    if (c == '.') {
        // Disconnected component — next atom has no prev
        consume();
        prev_atom_        = -1;
        pending_order_    = std::nullopt;
        pending_aromatic_ = false;
        return true;
    }

    // ---- Branch tokens ----
    if (c == '(') {
        consume();
        branch_stack_.push_back(prev_atom_);
        return true;
    }
    if (c == ')') {
        consume();
        if (branch_stack_.empty())
            throw SmilesParseError("unexpected ')'", pos_);
        prev_atom_ = branch_stack_.back();
        branch_stack_.pop_back();
        pending_order_    = std::nullopt;
        pending_aromatic_ = false;
        return true;
    }

    // ---- Ring closure: digit or %NN ----
    if (std::isdigit(c)) {
        consume();
        int ring_num = c - '0';
        BondOrder ring_bond = pending_order_.value_or(BondOrder::SINGLE);
        bool ring_arom      = pending_aromatic_;
        pending_order_      = std::nullopt;
        pending_aromatic_   = false;
        handle_ring_closure(ring_num, ring_bond, ring_arom);
        return true;
    }
    if (c == '%') {
        consume();
        if (!std::isdigit(peek()) || !std::isdigit(peek(1)))
            throw SmilesParseError("expected two digits after '%'", pos_);
        int d1 = consume() - '0';
        int d2 = consume() - '0';
        int ring_num = d1 * 10 + d2;
        BondOrder ring_bond = pending_order_.value_or(BondOrder::SINGLE);
        bool ring_arom      = pending_aromatic_;
        pending_order_      = std::nullopt;
        pending_aromatic_   = false;
        handle_ring_closure(ring_num, ring_bond, ring_arom);
        return true;
    }

    // ---- Bracket atom ----
    if (c == '[') {
        int new_atom = parse_bracket_atom();
        connect(new_atom);
        prev_atom_ = new_atom;
        return true;
    }

    // ---- Organic-subset or aromatic atom ----
    // Check two-character symbols first: Cl, Br
    if ((c == 'C' && peek(1) == 'l') || (c == 'B' && peek(1) == 'r')) {
        int new_atom = parse_organic_atom();
        connect(new_atom);
        prev_atom_ = new_atom;
        return true;
    }
    if (kOrganicSubset.count(c) || kAromaticSubset.count(c)) {
        int new_atom = parse_organic_atom();
        connect(new_atom);
        prev_atom_ = new_atom;
        return true;
    }

    // Whitespace is ignorable; any other unknown token is a hard parse error
    // (do not silently accept garbage like "not_a_smiles_!!!").
    if (std::isspace(static_cast<unsigned char>(c))) {
        consume();
        return true;
    }
    throw SmilesParseError(std::string("unknown SMILES token '") + c + "'", pos_);
}

// ---------------------------------------------------------------------------
// Organic-subset atom parser
// ---------------------------------------------------------------------------

int SmilesParser::parse_organic_atom() {
    char c = peek();

    // Two-character organic atoms
    if (c == 'C' && peek(1) == 'l') { consume(); consume();
        int idx = mol_.add_atom(Element::Cl); return idx; }
    if (c == 'B' && peek(1) == 'r') { consume(); consume();
        int idx = mol_.add_atom(Element::Br); return idx; }

    // Aromatic subset (lowercase)
    if (kAromaticSubset.count(c)) {
        consume();
        const auto& info = kAromaticSubset.at(c);
        int idx = mol_.add_atom(info.elem, 0, 0, true /*aromatic*/);
        return idx;
    }

    // Normal organic subset
    if (kOrganicSubset.count(c)) {
        consume();
        const auto& info = kOrganicSubset.at(c);
        int idx = mol_.add_atom(info.elem);
        return idx;
    }

    throw SmilesParseError(std::string("unexpected character '") + c + "'", pos_);
}

// ---------------------------------------------------------------------------
// Bracket atom parser: [isotope? symbol @? @? H count? charge? :mapnum?]
// ---------------------------------------------------------------------------

int SmilesParser::parse_bracket_atom() {
    assert(peek() == '[');
    consume(); // consume '['

    int bracket_pos = pos_;

    // Optional isotope
    int isotope = 0;
    if (std::isdigit(peek())) {
        isotope = parse_int();
    }

    // Element symbol
    bool is_aromatic = false;
    Element elem = parse_element(is_aromatic);
    if (elem == Element::Unknown && !is_aromatic) {
        // Special case: 'H' alone is valid
        // (already handled by parse_element)
        throw SmilesParseError("unknown element in bracket atom", bracket_pos);
    }

    // Optional chirality: @ or @@
    if (peek() == '@') {
        consume();
        if (peek() == '@') consume();
        // Chirality parsed, not enforced
    }

    // Optional hydrogen count
    int h_count = 0;
    if (peek() == 'H') {
        consume();
        if (std::isdigit(peek()))
            h_count = parse_int();
        else
            h_count = 1;
    }

    // Optional formal charge
    int formal_charge = 0;
    if (peek() == '+' || peek() == '-') {
        char sign = consume();
        if (std::isdigit(peek())) {
            int mag = parse_int();
            formal_charge = (sign == '+') ? mag : -mag;
        } else {
            // Count repeated +/-
            int mag = 1;
            while (peek() == sign) { consume(); ++mag; }
            formal_charge = (sign == '+') ? mag : -mag;
        }
    }

    // Optional atom-map number
    int atom_map = 0;
    if (peek() == ':') {
        consume();
        atom_map = parse_int();
    }

    if (peek() != ']')
        throw SmilesParseError("expected ']' closing bracket atom", pos_);
    consume(); // consume ']'

    int idx = mol_.add_atom(elem, formal_charge, isotope, is_aromatic);
    // Store explicit H count using the negative-flag convention understood by parse()
    // Negative means "explicit": -(h_count + 1)
    mol_.atoms[idx].implicit_h_count = -(h_count + 1);
    mol_.atoms[idx].atom_map_num     = atom_map;
    mol_.atoms[idx].isotope          = isotope;

    return idx;
}

// ---------------------------------------------------------------------------
// Element symbol parser (inside or outside brackets)
// ---------------------------------------------------------------------------

Element SmilesParser::parse_element(bool& is_aromatic_out) {
    is_aromatic_out = false;
    char c = peek();

    // Aromatic atom symbols (lowercase)
    if (std::islower(c)) {
        is_aromatic_out = true;
        // Map lowercase → element
        char upper = static_cast<char>(std::toupper(c));
        consume();
        // Check two-char first: 'se' 'te'
        if (upper == 'S' && peek() == 'e') { consume(); return Element::Se; }
        switch (upper) {
            case 'B': return Element::B;
            case 'C': return Element::C;
            case 'N': return Element::N;
            case 'O': return Element::O;
            case 'P': return Element::P;
            case 'S': return Element::S;
            default:  return Element::Unknown;
        }
    }

    // Uppercase: try two-char first
    if (std::isupper(c)) {
        char next = peek(1);
        // Two-character uppercase symbols
        if (next != '\0' && std::islower(next)) {
            char sym[3] = { c, next, '\0' };
            Element e = element_from_symbol(sym);
            if (e != Element::Unknown) {
                consume(); consume();
                return e;
            }
        }
        // Single character
        consume();
        return element_from_symbol(std::string(1, c));
    }

    // Wildcard '*' — treat as unknown
    if (c == '*') {
        consume();
        return Element::Unknown;
    }

    return Element::Unknown;
}

// ---------------------------------------------------------------------------
// Integer parser
// ---------------------------------------------------------------------------

int SmilesParser::parse_int() {
    int value = 0;
    while (std::isdigit(peek())) {
        value = value * 10 + (consume() - '0');
    }
    return value;
}

// ---------------------------------------------------------------------------
// Bond connection
// ---------------------------------------------------------------------------

void SmilesParser::connect(int new_atom) {
    if (prev_atom_ < 0) {
        // No previous atom (start of SMILES or after '.')
        pending_order_    = std::nullopt;
        pending_aromatic_ = false;
        return;
    }

    BondOrder order;
    bool      aromatic;
    if (pending_order_.has_value()) {
        order   = *pending_order_;
        aromatic = pending_aromatic_;
    } else {
        // Implicit bond: aromatic if both atoms are marked aromatic
        order   = default_bond_order(prev_atom_, new_atom);
        aromatic = (order == BondOrder::AROMATIC);
    }

    mol_.add_bond(prev_atom_, new_atom, order, aromatic);

    pending_order_    = std::nullopt;
    pending_aromatic_ = false;
}

// ---------------------------------------------------------------------------
// Ring closure
// ---------------------------------------------------------------------------

void SmilesParser::handle_ring_closure(int ring_num, BondOrder explicit_order,
                                        bool aromatic) {
    if (prev_atom_ < 0)
        throw SmilesParseError("ring closure without preceding atom", pos_);

    auto it = ring_open_.find(ring_num);
    if (it == ring_open_.end()) {
        // Opening a new ring
        ring_open_[ring_num] = { prev_atom_, explicit_order, aromatic };
    } else {
        // Closing an existing ring
        RingOpen& ro = it->second;
        if (ro.atom_idx == prev_atom_)
            throw SmilesParseError("ring closure to same atom", pos_);

        // Resolve bond order: both ends must agree (or one is implicit)
        BondOrder order;
        bool is_arom;
        if (explicit_order != BondOrder::SINGLE || aromatic) {
            order   = explicit_order;
            is_arom = aromatic;
        } else if (ro.order != BondOrder::SINGLE || ro.aromatic) {
            order   = ro.order;
            is_arom = ro.aromatic;
        } else {
            order   = default_bond_order(ro.atom_idx, prev_atom_);
            is_arom = (order == BondOrder::AROMATIC);
        }

        mol_.add_bond(ro.atom_idx, prev_atom_, order, is_arom);
        ring_open_.erase(it);
    }
}

// ---------------------------------------------------------------------------
// Default bond order (used for implicit bonds)
// ---------------------------------------------------------------------------

BondOrder SmilesParser::default_bond_order(int atom_i, int atom_j) const noexcept {
    bool i_arom = mol_.atoms[atom_i].is_aromatic;
    bool j_arom = mol_.atoms[atom_j].is_aromatic;
    if (i_arom && j_arom) return BondOrder::AROMATIC;
    return BondOrder::SINGLE;
}

// ---------------------------------------------------------------------------
// Implicit hydrogen count computation
// Based on OpenSMILES specification: H_implicit = normal_valence - (bond_order_sum + |formal_charge|)
// The smallest normal valence >= bond_order_sum is used.
// ---------------------------------------------------------------------------

int SmilesParser::compute_implicit_h(int atom_idx) const noexcept {
    const Atom& a = mol_.atoms[atom_idx];

    // Atoms that don't get implicit H: unknowns, metals
    if (a.element == Element::Unknown) return 0;
    if (a.element == Element::Fe || a.element == Element::Zn ||
        a.element == Element::Ca || a.element == Element::Mg) return 0;

    // Compute current bond-order sum
    [[maybe_unused]] float bos = mol_.bond_order_sum(atom_idx);
    // Aromatic bonds: each contributes 1.5 but for H computation we use 1
    // (OpenSMILES uses integer arithmetic: aromatic bonds count as 1 for valence)
    // Count aromatic bonds and re-sum
    int   arom_count  = 0;
    float exact_bos   = 0.0f;
    for (int bidx : mol_.bond_adj[atom_idx]) {
        const Bond& b = mol_.bonds[bidx];
        if (b.is_aromatic || b.order == BondOrder::AROMATIC) {
            exact_bos += 1.5f;
            ++arom_count;
        } else {
            exact_bos += static_cast<float>(static_cast<int>(b.order));
        }
    }
    // For SMILES H counting, aromatic bonds count as 1
    float bos_for_h = exact_bos - 0.5f * arom_count;

    // Expected valences per element
    int normal_valences[4] = {-1,-1,-1,-1};
    switch (a.element) {
        case Element::B:  normal_valences[0]=3; break;
        case Element::C:  normal_valences[0]=4; break;
        case Element::N:  normal_valences[0]=3; normal_valences[1]=5; break;
        case Element::O:  normal_valences[0]=2; break;
        case Element::P:  normal_valences[0]=3; normal_valences[1]=5; break;
        case Element::S:  normal_valences[0]=2; normal_valences[1]=4; normal_valences[2]=6; break;
        case Element::F:  normal_valences[0]=1; break;
        case Element::Cl: normal_valences[0]=1; normal_valences[1]=3;
                          normal_valences[2]=5; normal_valences[3]=7; break;
        case Element::Br: normal_valences[0]=1; normal_valences[1]=3;
                          normal_valences[2]=5; break;
        case Element::I:  normal_valences[0]=1; normal_valences[1]=3;
                          normal_valences[2]=5; normal_valences[3]=7; break;
        default: return 0; // no implicit H for other elements
    }

    int target = -1;
    int bos_int = static_cast<int>(std::ceil(bos_for_h));
    for (int nv : normal_valences) {
        if (nv < 0) break;
        if (nv >= bos_int) { target = nv; break; }
    }
    if (target < 0) return 0; // over-valenced, no implicit H

    int h = target - bos_int;
    if (h < 0) h = 0;
    return h;
}

// ---------------------------------------------------------------------------
// Factory function — BonMol::from_smiles delegating to SmilesParser
// ---------------------------------------------------------------------------

BonMol from_smiles(const std::string& smiles) {
    SmilesParser parser;
    SmilesParseResult result = parser.parse(smiles);
    return std::move(result.mol);
}

} // namespace bonmol
