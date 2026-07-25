// shell_exec.h — Path-safe shell quoting + argv-based process helpers
//
// Prefer fork/execvp (or SubprocessGuard::fork_exec_argv) with argv arrays.
// When a shell is unavoidable (env KEY=VAL prefixes, stdout/stderr
// redirection), interpolate paths only through shell_quote() after
// validate_exec_path() succeeds.
//
// Security focus:
//   - Command injection via FLEXAIDDS_DATA_DIR or malicious manifest paths
//   - Provenance integrity (receipts must record real resolved paths, not
//     shell-mutated strings)
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
#pragma once

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#ifndef _WIN32
#include <fcntl.h>
#include <sys/select.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace flexaids {
namespace shell_exec {

/// Why a path was refused for shell interpolation or exec argv.
enum class PathReject : int {
    Ok = 0,
    Empty,
    ContainsNul,
    ContainsNewline,   // '\n' or '\r' — breaks sh -c line structure
    ContainsControl,   // other C0 controls (except TAB, which is allowed)
};

[[nodiscard]] inline const char* path_reject_message(PathReject r) {
    switch (r) {
        case PathReject::Ok:               return "ok";
        case PathReject::Empty:            return "path is empty";
        case PathReject::ContainsNul:      return "path contains NUL byte";
        case PathReject::ContainsNewline:  return "path contains newline/CR";
        case PathReject::ContainsControl:  return "path contains control character";
    }
    return "unknown path reject";
}

/// Inspect *path* for characters that must never reach a shell string or that
/// would corrupt path provenance. TAB is allowed (rare but valid in paths on
/// some filesystems); NUL, CR, LF, and other C0 controls are refused.
[[nodiscard]] inline PathReject validate_exec_path(std::string_view path) {
    if (path.empty()) return PathReject::Empty;
    for (unsigned char c : path) {
        if (c == '\0') return PathReject::ContainsNul;
        if (c == '\n' || c == '\r') return PathReject::ContainsNewline;
        // C0 controls except TAB (0x09)
        if (c < 0x20 && c != '\t') return PathReject::ContainsControl;
        if (c == 0x7f) return PathReject::ContainsControl;
    }
    return PathReject::Ok;
}

[[nodiscard]] inline bool is_safe_exec_path(std::string_view path) {
    return validate_exec_path(path) == PathReject::Ok;
}

/// POSIX-shell single-quote escape: `foo'bar` → `'foo'\''bar'`.
/// Does **not** validate; use shell_quote() for validated quoting.
[[nodiscard]] inline std::string shell_quote_raw(std::string_view s) {
    std::string o;
    o.reserve(s.size() + 8);
    o.push_back('\'');
    for (char c : s) {
        if (c == '\'') {
            // end quote, escaped literal quote, reopen quote
            o += "'\\''";
        } else {
            o.push_back(c);
        }
    }
    o.push_back('\'');
    return o;
}

/// Validate then single-quote. Returns nullopt if the path is unsafe.
[[nodiscard]] inline std::optional<std::string>
shell_quote_checked(std::string_view path) {
    if (validate_exec_path(path) != PathReject::Ok) return std::nullopt;
    return shell_quote_raw(path);
}

/// Validate then single-quote. Throws std::invalid_argument on unsafe paths.
[[nodiscard]] inline std::string shell_quote(std::string_view path) {
    PathReject r = validate_exec_path(path);
    if (r != PathReject::Ok) {
        throw std::invalid_argument(
            std::string("unsafe path for shell: ") + path_reject_message(r));
    }
    return shell_quote_raw(path);
}

/// Result of capturing a child process stdout/stderr (no shell).
struct CapturedOutput {
    int         exit_code = -1;
    std::string stdout_text;
    std::string stderr_text;  // filled when discard_stderr=false
    bool        ok = false;  // true if process started and exited normally
};

/// Run argv via fork+execvp with no shell. Captures stdout always.
/// If discard_stderr=true, child stderr → /dev/null.
/// If discard_stderr=false, child stderr is captured into stderr_text.
/// Returns empty-ok=false on fork/exec failure. On Windows, returns ok=false.
[[nodiscard]] inline CapturedOutput
run_argv_capture(const std::vector<std::string>& argv,
                 bool discard_stderr = true) {
    CapturedOutput out;
    if (argv.empty() || argv[0].empty()) return out;
    for (const auto& a : argv) {
        if (validate_exec_path(a) != PathReject::Ok) return out;
    }
#ifndef _WIN32
    int out_pipe[2];
    int err_pipe[2] = {-1, -1};
    if (::pipe(out_pipe) != 0) return out;
    if (!discard_stderr) {
        if (::pipe(err_pipe) != 0) {
            ::close(out_pipe[0]);
            ::close(out_pipe[1]);
            return out;
        }
    }

    pid_t pid = ::fork();
    if (pid < 0) {
        ::close(out_pipe[0]);
        ::close(out_pipe[1]);
        if (!discard_stderr) {
            ::close(err_pipe[0]);
            ::close(err_pipe[1]);
        }
        return out;
    }
    if (pid == 0) {
        ::close(out_pipe[0]);
        ::dup2(out_pipe[1], STDOUT_FILENO);
        ::close(out_pipe[1]);
        if (discard_stderr) {
            int devnull = ::open("/dev/null", O_WRONLY);
            if (devnull >= 0) {
                ::dup2(devnull, STDERR_FILENO);
                ::close(devnull);
            }
        } else {
            ::close(err_pipe[0]);
            ::dup2(err_pipe[1], STDERR_FILENO);
            ::close(err_pipe[1]);
        }
        ::close(STDIN_FILENO);
        int devnull_in = ::open("/dev/null", O_RDONLY);
        if (devnull_in >= 0 && devnull_in != STDIN_FILENO) {
            ::dup2(devnull_in, STDIN_FILENO);
            ::close(devnull_in);
        }

        std::vector<char*> c_argv;
        c_argv.reserve(argv.size() + 1);
        for (const auto& s : argv) {
            c_argv.push_back(const_cast<char*>(s.c_str()));
        }
        c_argv.push_back(nullptr);
        ::execvp(c_argv[0], c_argv.data());
        ::_exit(127);
    }

    ::close(out_pipe[1]);
    if (!discard_stderr) ::close(err_pipe[1]);

    // Drain both pipes to avoid child blocking on a full stderr pipe.
    char buf[4096];
    bool out_open = true;
    bool err_open = !discard_stderr;
    while (out_open || err_open) {
        fd_set rfds;
        FD_ZERO(&rfds);
        int maxfd = -1;
        if (out_open) {
            FD_SET(out_pipe[0], &rfds);
            maxfd = out_pipe[0];
        }
        if (err_open) {
            FD_SET(err_pipe[0], &rfds);
            if (err_pipe[0] > maxfd) maxfd = err_pipe[0];
        }
        int sel = ::select(maxfd + 1, &rfds, nullptr, nullptr, nullptr);
        if (sel < 0) {
            if (errno == EINTR) continue;
            break;
        }
        if (out_open && FD_ISSET(out_pipe[0], &rfds)) {
            ssize_t n = ::read(out_pipe[0], buf, sizeof(buf));
            if (n > 0) {
                out.stdout_text.append(buf, static_cast<size_t>(n));
            } else {
                out_open = false;
            }
        }
        if (err_open && FD_ISSET(err_pipe[0], &rfds)) {
            ssize_t n = ::read(err_pipe[0], buf, sizeof(buf));
            if (n > 0) {
                out.stderr_text.append(buf, static_cast<size_t>(n));
            } else {
                err_open = false;
            }
        }
    }
    ::close(out_pipe[0]);
    if (!discard_stderr) ::close(err_pipe[0]);

    int status = 0;
    while (::waitpid(pid, &status, 0) == -1) {
        if (errno != EINTR) {
            out.ok = false;
            return out;
        }
    }
    if (WIFEXITED(status)) {
        out.exit_code = WEXITSTATUS(status);
        out.ok = true;
    } else {
        out.exit_code = -1;
        out.ok = false;
    }
    return out;
#else
    (void)discard_stderr;
    return out;
#endif
}

/// First whitespace-delimited token of stdout from run_argv_capture, or "".
[[nodiscard]] inline std::string
run_argv_first_token(const std::vector<std::string>& argv) {
    auto cap = run_argv_capture(argv, /*discard_stderr=*/true);
    if (!cap.ok || cap.exit_code != 0) return {};
    std::string tok;
    for (char c : cap.stdout_text) {
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
            if (!tok.empty()) break;
        } else {
            tok.push_back(c);
        }
    }
    return tok;
}

/// Run argv and wait for exit code (no shell, no capture). -1 on failure.
[[nodiscard]] inline int run_argv(const std::vector<std::string>& argv) {
    if (argv.empty() || argv[0].empty()) return -1;
    for (const auto& a : argv) {
        if (validate_exec_path(a) != PathReject::Ok) return -1;
    }
#ifndef _WIN32
    pid_t pid = ::fork();
    if (pid < 0) return -1;
    if (pid == 0) {
        ::close(STDIN_FILENO);
        int devnull = ::open("/dev/null", O_RDONLY);
        if (devnull >= 0 && devnull != STDIN_FILENO) {
            ::dup2(devnull, STDIN_FILENO);
            ::close(devnull);
        }
        std::vector<char*> c_argv;
        c_argv.reserve(argv.size() + 1);
        for (const auto& s : argv) {
            c_argv.push_back(const_cast<char*>(s.c_str()));
        }
        c_argv.push_back(nullptr);
        ::execvp(c_argv[0], c_argv.data());
        ::_exit(127);
    }
    int status = 0;
    while (::waitpid(pid, &status, 0) == -1) {
        if (errno != EINTR) return -1;
    }
    return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
#else
    return -1;
#endif
}

}  // namespace shell_exec
}  // namespace flexaids
