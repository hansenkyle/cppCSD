// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef CLI_H
#define CLI_H

#include <filesystem>
#include <optional>

// Parses ldcsd's command-line arguments: a single required positional
// argument, the path to the input deck YAML file. Doesn't include CLI11's
// own header -- that's an implementation detail of cli.cpp -- so this stays
// a plain function callers (main, tests) can use without depending on CLI11
// themselves.
//
// Returns the parsed path on success. On --help or a parse error
// (missing/extra arguments, etc.), CLI11 has already printed its usage or
// error message to stdout/stderr; this returns std::nullopt and sets
// exit_code to what the caller should return from main() in that case,
// rather than calling std::exit() itself -- so the parsing logic stays
// testable without ending the test process.
std::optional<std::filesystem::path> parseArgs(int argc, char** argv, int& exit_code);

#endif
