/* ============================================================================
   FlexAID∆S — terminal colour, matching the design system
   Le Bonhomme Pharma · Montréal

   The palette here is the SAME seven key colors the website and the Design
   System use, in the same bindings. A key color is never reassigned to a
   different quantity — that binding IS the system:

       mint       #45E0A8   ΔH · enthalpy · brand primary
       violet     #8B5CF6   ΔS · configurational entropy
       tangerine  #FF9300   ΔG · free energy · results
       firetruck  #F5232B   T  · temperature · and failure
       aqua       #00A2FF   ΔS_vib · vibrational · tENCoM
       strawberry #FF2F92   receptor · pocket · warnings
       magnesium  #DCDCE4   apo baseline · reference

   SERIES ramp order is ENERGY along the binding coordinate, not hue:
       apo → unbound ΔS → pocket contact → rigidification ΔS_vib
           → contacts ΔH → converged ΔG
   so a staged run reads as a reaction path. Firetruck stays OUT of the
   stage ramp: it is a scalar and a failure signal, never a data class.

   ── WHY NOTHING IS EVER WRITTEN TO A FILE ─────────────────────────────────
   Every accessor returns "" unless stdout is an interactive terminal. So
   `flexaidds ... > run.log`, a pipe, a CI job, or any redirect produces
   byte-identical output to before this header existed. Colour is a property
   of the terminal, never of the data. Result files (MOL2, surface, grand-
   canonical summaries) are written through their own FILE* handles and are
   not touched by anything here.

   Honours, in precedence order:
     FLEXAIDDS_NO_COLOR / NO_COLOR set  -> off  (https://no-color.org)
     TERM=dumb or unset                 -> off
     FORCE_COLOR / CLICOLOR_FORCE set   -> on, even when piped
     stdout is not a TTY                -> off
   `--no-color` on the command line calls tui::set_enabled(false).

   Depth degrades: 24-bit when COLORTERM says truecolor, else the nearest
   xterm-256 index, so the hues stay recognisable on older terminals.

   Header-only. No dependencies beyond libc.
   ============================================================================ */
#ifndef FLEXAIDDS_TUI_COLOR_H
#define FLEXAIDDS_TUI_COLOR_H

#include <cstdio>
#include <cstdlib>
#include <cstring>

#if defined(_WIN32)
#  include <io.h>
#  include <windows.h>
#  define FLEXAIDDS_ISATTY(fd) _isatty(fd)
#  define FLEXAIDDS_FILENO(f)  _fileno(f)
#else
#  include <unistd.h>
#  define FLEXAIDDS_ISATTY(fd) isatty(fd)
#  define FLEXAIDDS_FILENO(f)  fileno(f)
#endif

namespace tui {

/* ---- capability detection (computed once, then cached) ------------------- */

enum class Depth { None, Ansi256, TrueColor };

inline Depth& depth_ref() { static Depth d = Depth::None; return d; }
inline bool&  ready_ref() { static bool r = false;         return r; }

inline bool env_set(const char* k) {
    const char* v = std::getenv(k);
    return v != nullptr && v[0] != '\0';
}

inline Depth detect_depth() {
    /* NO_COLOR is a hard veto and takes precedence over everything. */
    if (env_set("FLEXAIDDS_NO_COLOR") || std::getenv("NO_COLOR") != nullptr)
        return Depth::None;

    const char* term = std::getenv("TERM");
    if (term == nullptr || std::strcmp(term, "dumb") == 0)
        return Depth::None;

    /* A forced request wins over the TTY test, so CI can still capture colour
       deliberately. Without it, a non-terminal stdout is always plain. */
    const bool forced = env_set("FORCE_COLOR") || env_set("CLICOLOR_FORCE");
    if (!forced && !FLEXAIDDS_ISATTY(FLEXAIDDS_FILENO(stdout)))
        return Depth::None;

#if defined(_WIN32)
    /* Windows terminals need VT sequences switched on explicitly. */
    HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode = 0;
    if (h != INVALID_HANDLE_VALUE && GetConsoleMode(h, &mode))
        SetConsoleMode(h, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
#endif

    const char* ct = std::getenv("COLORTERM");
    if (ct && (std::strstr(ct, "truecolor") || std::strstr(ct, "24bit")))
        return Depth::TrueColor;
    if (std::strstr(term, "256color") || std::strstr(term, "kitty") ||
        std::strstr(term, "alacritty") || std::strstr(term, "wezterm"))
        return Depth::Ansi256;
    /* Modern terminals overwhelmingly do truecolor even when COLORTERM is
       unset; 256 is the safe floor for anything that got this far. */
    return Depth::Ansi256;
}

inline Depth depth() {
    if (!ready_ref()) { depth_ref() = detect_depth(); ready_ref() = true; }
    return depth_ref();
}

inline bool enabled() { return depth() != Depth::None; }

/** Force colour off (for `--no-color`) or on. Call before any output. */
inline void set_enabled(bool on) {
    ready_ref() = true;
    depth_ref() = on ? (detect_depth() == Depth::TrueColor ? Depth::TrueColor
                                                           : Depth::Ansi256)
                     : Depth::None;
}

/* ---- the palette -------------------------------------------------------- */
/* Each entry carries its 24-bit hex and the nearest xterm-256 index. */

struct Hue { unsigned char r, g, b; unsigned char x256; };

inline constexpr Hue MINT       {0x45, 0xE0, 0xA8, 79};   /* ΔH        */
inline constexpr Hue VIOLET     {0x8B, 0x5C, 0xF6, 99};   /* ΔS        */
inline constexpr Hue TANGERINE  {0xFF, 0x93, 0x00, 208};  /* ΔG        */
inline constexpr Hue FIRETRUCK  {0xF5, 0x23, 0x2B, 196};  /* T · fail  */
inline constexpr Hue AQUA       {0x00, 0xA2, 0xFF, 39};   /* ΔS_vib    */
inline constexpr Hue STRAWBERRY {0xFF, 0x2F, 0x92, 198};  /* receptor  */
inline constexpr Hue MAGNESIUM  {0xDC, 0xDC, 0xE4, 253};  /* baseline  */
inline constexpr Hue FG         {0xE4, 0xE3, 0xF5, 189};  /* body text */
inline constexpr Hue MUTED      {0x8D, 0x8C, 0xB0, 103};  /* secondary */
inline constexpr Hue FAILTEXT   {0xFF, 0x6B, 0x6B, 203};  /* small fail labels */

/* Returns a pointer to a per-hue static buffer, or "" when colour is off.
   One buffer per hue, so several colours can be live in a single printf. */
template <const Hue* H>
inline const char* seq() {
    static char buf[24];
    static Depth built = Depth::None;
    const Depth d = depth();
    if (d == Depth::None) return "";
    if (built != d) {
        if (d == Depth::TrueColor)
            std::snprintf(buf, sizeof buf, "\033[38;2;%u;%u;%um",
                          (unsigned)H->r, (unsigned)H->g, (unsigned)H->b);
        else
            std::snprintf(buf, sizeof buf, "\033[38;5;%um", (unsigned)H->x256);
        built = d;
    }
    return buf;
}

inline const char* reset()      { return enabled() ? "\033[0m" : ""; }
inline const char* bold()       { return enabled() ? "\033[1m" : ""; }
inline const char* dim()        { return enabled() ? "\033[2m" : ""; }

inline const char* mint()       { return seq<&MINT>(); }
inline const char* violet()     { return seq<&VIOLET>(); }
inline const char* tangerine()  { return seq<&TANGERINE>(); }
inline const char* firetruck()  { return seq<&FIRETRUCK>(); }
inline const char* aqua()       { return seq<&AQUA>(); }
inline const char* strawberry() { return seq<&STRAWBERRY>(); }
inline const char* magnesium()  { return seq<&MAGNESIUM>(); }
inline const char* fg()         { return seq<&FG>(); }
inline const char* muted()      { return seq<&MUTED>(); }
inline const char* failtext()   { return seq<&FAILTEXT>(); }

/* Quantity aliases — prefer these at call sites, so the binding of colour to
   thermodynamic quantity stays legible in the source itself. */
inline const char* dH()    { return mint();      }
inline const char* dS()    { return violet();    }
inline const char* dG()    { return tangerine(); }
inline const char* T()     { return firetruck(); }
inline const char* dSvib() { return aqua();      }
inline const char* pocket(){ return strawberry();}
inline const char* apo()   { return magnesium(); }

/* ---- the series ramp ---------------------------------------------------- */
/* Energy order along the binding coordinate. Index 1..6; anything outside
   falls back to muted rather than wrapping into a wrong quantity. */
inline const char* series(int step) {
    switch (step) {
        case 1: return magnesium();   /* apo baseline            */
        case 2: return violet();      /* unbound · ΔS dominates  */
        case 3: return strawberry();  /* first pocket contact    */
        case 4: return aqua();        /* rigidification · ΔS_vib */
        case 5: return mint();        /* contacts formed · ΔH    */
        case 6: return tangerine();   /* converged · ΔG          */
        default: return muted();
    }
}

/* ---- composed elements, mirroring the website TUI ------------------------ */

/** Display width in terminal columns, not bytes. Labels carry Greek (Δ) and
    `·`, so printf's %-Ns would pad by BYTE count and shear the detail column
    to the left by however many continuation bytes the label happened to hold.
    Counts UTF-8 lead bytes; combining marks (U+0300..U+036F) take no column. */
inline int display_width(const char* s) {
    int w = 0;
    for (const unsigned char* p = (const unsigned char*)s; *p; ++p) {
        if ((*p & 0xC0) == 0x80) continue;              /* continuation byte */
        if ((p[0] & 0xFE) == 0xCC) { ++p; continue; }   /* U+0300..U+036F    */
        ++w;
    }
    return w;
}

inline void pad_to(int width, int used, std::FILE* out) {
    for (int i = used; i < width; ++i) std::fputc(' ', out);
}

/** `FlexAID∆S` wordmark: mint mark, tangerine ∆S — as everywhere else. */
inline void brand(std::FILE* out = stdout) {
    std::fprintf(out, "%s%sFlexAID%s%s∆S%s",
                 bold(), mint(), reset(), tangerine(), reset());
}

/** A violet section rule, the terminal form of the site's section divider. */
inline void rule(const char* title, std::FILE* out = stdout) {
    std::fprintf(out, "%s──%s %s%s%s %s", violet(), reset(),
                 bold(), title ? title : "", reset(), violet());
    int pad = 58 - display_width(title ? title : "");
    for (int i = 0; i < pad; ++i) std::fputc('-', out);
    std::fprintf(out, "%s\n", reset());
}

/** `[n/total] label            detail` — one step of the binding coordinate,
    coloured by its place on the series ramp. */
inline void stage(int n, int total, const char* label, const char* detail,
                  std::FILE* out = stdout) {
    const char* lbl = label ? label : "";
    std::fprintf(out, "  %s[%d/%d]%s %s%s%s",
                 muted(), n, total, reset(), series(n), lbl, reset());
    pad_to(28, display_width(lbl), out);
    std::fprintf(out, " %s%s%s\n", muted(), detail ? detail : "", reset());
}

/** The identity, read term by term in colour. */
inline void equation(unsigned temperature_K = 0, std::FILE* out = stdout) {
    std::fprintf(out, "  %sΔG%s = %sΔH%s − %sT%s%sΔS%s − %sT%s%sΔS_vib%s",
                 dG(), reset(), dH(), reset(), T(), reset(), dS(), reset(),
                 T(), reset(), dSvib(), reset());
    if (temperature_K) std::fprintf(out, "        %sT = %u K%s", T(), temperature_K, reset());
    std::fputc('\n', out);
}

/* Status lines. Severity reads by brightness: mint pass, strawberry caution,
   firetruck stop — never yellow. */
inline void ok(const char* msg, std::FILE* out = stdout) {
    std::fprintf(out, "  %s●%s %s%s%s\n", mint(), reset(), fg(), msg, reset());
}
inline void warn(const char* msg, std::FILE* out = stdout) {
    std::fprintf(out, "  %s● WARN%s %s%s%s\n", strawberry(), reset(), fg(), msg, reset());
}
inline void fail(const char* msg, std::FILE* out = stderr) {
    std::fprintf(out, "  %s● FAIL%s %s%s%s\n", failtext(), reset(), fg(), msg, reset());
}
/** Converged — the run reached ΔG. Tangerine, the end of the ramp. */
inline void converged(const char* msg, std::FILE* out = stdout) {
    std::fprintf(out, "  %s● CONVERGED%s %s%s%s\n",
                 tangerine(), reset(), muted(), msg ? msg : "", reset());
}

/** A `[TAG]` prefix in a chosen hue, matching the site's eyebrow treatment. */
inline void tag(const char* name, const char* hue, std::FILE* out = stdout) {
    std::fprintf(out, "%s[%s]%s ", hue, name, reset());
}

} // namespace tui

#endif /* FLEXAIDDS_TUI_COLOR_H */
