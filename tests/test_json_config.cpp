// tests/test_json_config.cpp — Unit tests for json_value.h (lightweight JSON parser)
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/json_value.h"
#include "../LIB/config_defaults.h"
#include <cmath>
#include <fstream>

static constexpr double EPSILON = 1e-9;

// ===========================================================================
// Basic value types
// ===========================================================================

TEST(JsonValue, NullDefault) {
    json::Value v;
    EXPECT_TRUE(v.is_null());
    EXPECT_FALSE(v.is_bool());
    EXPECT_FALSE(v.is_int());
}

TEST(JsonValue, BoolValues) {
    json::Value t(true), f(false);
    EXPECT_TRUE(t.is_bool());
    EXPECT_TRUE(t.as_bool());
    EXPECT_FALSE(f.as_bool());
}

TEST(JsonValue, IntValues) {
    json::Value v(42);
    EXPECT_TRUE(v.is_int());
    EXPECT_TRUE(v.is_number());
    EXPECT_EQ(v.as_int(), 42);
    EXPECT_NEAR(v.as_double(), 42.0, EPSILON);
}

TEST(JsonValue, DoubleValues) {
    json::Value v(3.14);
    EXPECT_TRUE(v.is_double());
    EXPECT_TRUE(v.is_number());
    EXPECT_NEAR(v.as_double(), 3.14, EPSILON);
    EXPECT_EQ(v.as_int(), 3);  // truncation
}

TEST(JsonValue, StringValues) {
    json::Value v("hello");
    EXPECT_TRUE(v.is_string());
    EXPECT_EQ(v.as_string(), "hello");
}

TEST(JsonValue, ArrayValues) {
    json::Array arr = {json::Value(1), json::Value(2), json::Value(3)};
    json::Value v(arr);
    EXPECT_TRUE(v.is_array());
    EXPECT_EQ(v.size(), 3u);
    EXPECT_EQ(v[static_cast<size_t>(0)].as_int(), 1);
    EXPECT_EQ(v[static_cast<size_t>(2)].as_int(), 3);
}

TEST(JsonValue, ObjectValues) {
    json::Object obj = {{"key", json::Value("val")}, {"num", json::Value(99)}};
    json::Value v(obj);
    EXPECT_TRUE(v.is_object());
    EXPECT_TRUE(v.contains("key"));
    EXPECT_FALSE(v.contains("missing"));
    EXPECT_EQ(v["key"].as_string(), "val");
    EXPECT_EQ(v["num"].as_int(), 99);
}

TEST(JsonValue, MissingKeyReturnsNull) {
    json::Object obj = {{"a", json::Value(1)}};
    json::Value v(obj);
    EXPECT_TRUE(v["nonexistent"].is_null());
    // Chained access on null returns null
    EXPECT_TRUE(v["nonexistent"]["deep"].is_null());
}

TEST(JsonValue, DefaultFallbacks) {
    json::Value null_v;
    EXPECT_EQ(null_v.as_bool(true), true);
    EXPECT_EQ(null_v.as_int(42), 42);
    EXPECT_NEAR(null_v.as_double(1.5), 1.5, EPSILON);
    EXPECT_EQ(null_v.as_string("fb"), "fb");
}

// ===========================================================================
// JSON parser
// ===========================================================================

TEST(JsonParser, EmptyObject) {
    auto v = json::parse("{}");
    EXPECT_TRUE(v.is_object());
    EXPECT_EQ(v.size(), 0u);
}

TEST(JsonParser, EmptyArray) {
    auto v = json::parse("[]");
    EXPECT_TRUE(v.is_array());
    EXPECT_EQ(v.size(), 0u);
}

TEST(JsonParser, SimpleObject) {
    auto v = json::parse(R"({"name": "FlexAID", "version": 2, "active": true})");
    EXPECT_EQ(v["name"].as_string(), "FlexAID");
    EXPECT_EQ(v["version"].as_int(), 2);
    EXPECT_TRUE(v["active"].as_bool());
}

TEST(JsonParser, NestedObject) {
    auto v = json::parse(R"({
        "scoring": {
            "function": "VCT",
            "weight": 1.5
        }
    })");
    EXPECT_EQ(v["scoring"]["function"].as_string(), "VCT");
    EXPECT_NEAR(v["scoring"]["weight"].as_double(), 1.5, EPSILON);
}

TEST(JsonParser, ArrayOfNumbers) {
    auto v = json::parse("[1.0, 0.5, 1.0, 0.5]");
    EXPECT_EQ(v.size(), 4u);
    EXPECT_NEAR(v[static_cast<size_t>(0)].as_double(), 1.0, EPSILON);
    EXPECT_NEAR(v[static_cast<size_t>(1)].as_double(), 0.5, EPSILON);
}

TEST(JsonParser, NullAndBooleans) {
    auto v = json::parse(R"({"a": null, "b": true, "c": false})");
    EXPECT_TRUE(v["a"].is_null());
    EXPECT_TRUE(v["b"].as_bool());
    EXPECT_FALSE(v["c"].as_bool());
}

TEST(JsonParser, NegativeAndScientific) {
    auto v = json::parse(R"({"neg": -42, "sci": 1.5e-3, "sciE": 2E10})");
    EXPECT_EQ(v["neg"].as_int(), -42);
    EXPECT_NEAR(v["sci"].as_double(), 0.0015, 1e-7);
    EXPECT_NEAR(v["sciE"].as_double(), 2e10, 1.0);
}

TEST(JsonParser, StringEscapes) {
    auto v = json::parse(R"({"msg": "hello\nworld\t\"quoted\""})");
    EXPECT_EQ(v["msg"].as_string(), "hello\nworld\t\"quoted\"");
}

TEST(JsonParser, LineComments) {
    // Extension: // comments are allowed in config files
    auto v = json::parse(R"({
        // This is a comment
        "key": 42
        // Another comment
    })");
    EXPECT_EQ(v["key"].as_int(), 42);
}

TEST(JsonParser, TrailingContentThrows) {
    EXPECT_THROW(json::parse("{} extra"), std::runtime_error);
}

TEST(JsonParser, InvalidJsonThrows) {
    EXPECT_THROW(json::parse("{invalid}"), std::runtime_error);
    EXPECT_THROW(json::parse(""), std::runtime_error);
}

// ===========================================================================
// File parsing
// ===========================================================================

TEST(JsonParser, ParseFile) {
    // Write a temp file
    const char* path = "/tmp/test_flexaid_config.json";
    {
        std::ofstream f(path);
        f << R"({
            "thermodynamics": {"temperature": 310},
            "ga": {"num_chromosomes": 2000}
        })";
    }
    auto v = json::parse_file(path);
    EXPECT_EQ(v["thermodynamics"]["temperature"].as_int(), 310);
    EXPECT_EQ(v["ga"]["num_chromosomes"].as_int(), 2000);
    std::remove(path);
}

TEST(JsonParser, ParseFileMissing) {
    EXPECT_THROW(json::parse_file("/nonexistent/path.json"), std::runtime_error);
}

// ===========================================================================
// Merge
// ===========================================================================

TEST(JsonMerge, ShallowOverride) {
    auto a = json::parse(R"({"x": 1, "y": 2})");
    auto b = json::parse(R"({"y": 99, "z": 3})");
    auto m = json::merge(a, b);
    EXPECT_EQ(m["x"].as_int(), 1);
    EXPECT_EQ(m["y"].as_int(), 99);  // overridden
    EXPECT_EQ(m["z"].as_int(), 3);   // added
}

TEST(JsonMerge, DeepMerge) {
    auto a = json::parse(R"({"scoring": {"function": "VCT", "weight": 1.0}})");
    auto b = json::parse(R"({"scoring": {"weight": 2.5}})");
    auto m = json::merge(a, b);
    EXPECT_EQ(m["scoring"]["function"].as_string(), "VCT");  // preserved
    EXPECT_NEAR(m["scoring"]["weight"].as_double(), 2.5, EPSILON);  // overridden
}

TEST(JsonMerge, NonObjectBOverridesA) {
    auto a = json::parse(R"({"x": {"nested": true}})");
    auto b = json::parse(R"({"x": 42})");
    auto m = json::merge(a, b);
    EXPECT_EQ(m["x"].as_int(), 42);  // object replaced by scalar
}

// ===========================================================================
// Mutable access (for building configs programmatically)
// ===========================================================================

TEST(JsonValue, MutableSet) {
    json::Value v;
    v.set("temp", json::Value(300));
    v.set("name", json::Value("test"));
    EXPECT_TRUE(v.is_object());
    EXPECT_EQ(v["temp"].as_int(), 300);
    EXPECT_EQ(v["name"].as_string(), "test");
}

// ===========================================================================
// Config defaults integration
// ===========================================================================

TEST(ConfigDefaults, DefaultsAreValid) {
    auto cfg = flexaid_default_config();
    EXPECT_TRUE(cfg.is_object());
    EXPECT_EQ(cfg["thermodynamics"]["temperature"].as_int(), 300);
    EXPECT_EQ(cfg["ga"]["num_chromosomes"].as_int(), 1000);
    EXPECT_EQ(cfg["scoring"]["function"].as_string(), "VCT");
    EXPECT_TRUE(cfg["flexibility"]["ligand_torsions"].as_bool());
    EXPECT_FALSE(cfg["flexibility"]["use_flexdee"].as_bool());
    EXPECT_NEAR(cfg["flexibility"]["soft_wall_cutoff"].as_double(), 0.40, EPSILON);

    // Rigid overrides
    auto rigid = flexaid_rigid_overrides();
    auto merged = json::merge(cfg, rigid);
    EXPECT_FALSE(merged["flexibility"]["ligand_torsions"].as_bool());
    EXPECT_EQ(merged["thermodynamics"]["temperature"].as_int(), 0);
    // Non-overridden values should persist
    EXPECT_EQ(merged["ga"]["num_chromosomes"].as_int(), 1000);
}
