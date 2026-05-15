#pragma once

#include <cctype>
#include <cstddef>
#include <cstdio>
#include <cstring>

inline void flexaids_copy_config_value(
    const char* line,
    char* dest,
    std::size_t dest_size)
{
    if (!dest || dest_size == 0) return;
    dest[0] = '\0';
    if (!line) return;

    const unsigned char* p = reinterpret_cast<const unsigned char*>(line);
    while (*p && !std::isspace(*p)) ++p;
    while (*p && std::isspace(*p)) ++p;

    std::snprintf(dest, dest_size, "%s", reinterpret_cast<const char*>(p));

    std::size_t len = std::strlen(dest);
    while (len > 0 && (dest[len - 1] == '\n' || dest[len - 1] == '\r')) {
        dest[--len] = '\0';
    }
}
