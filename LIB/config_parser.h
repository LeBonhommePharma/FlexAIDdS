// config_parser.h — JSON config file loader for FlexAIDdS
//
// Loads a user JSON config, merges with defaults, and applies
// the result to FA_Global / GB_Global structs.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "ProtocolConfig.h"
#include "json_value.h"
#include <string>

// Forward declarations
struct FA_Global_struct;
typedef struct FA_Global_struct FA_Global;
struct GB_Global_struct;
typedef struct GB_Global_struct GB_Global;

// Load a JSON config file and merge on top of defaults.
// If config_path is empty, returns defaults only.
json::Value load_config(const std::string& config_path);

// Apply merged JSON config to FA_Global and GB_Global structs.
// Optional `protocol` is a pre-snapshotted ProtocolConfig (no mid-apply getenv).
// When null, apply_config takes a single ProtocolConfig::from_env() snapshot at
// entry so env overrides stay dual-protocol-safe for this apply path.
void apply_config(const json::Value& config, FA_Global* FA, GB_Global* GB,
                  const flexaids::ProtocolConfig* protocol = nullptr);
