// First line
// Second line
// WHAT IS ABOVE THIS LINE IS USED BY TESTS AND MUST NOT BE CHANGED !
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "external/doctest.h"

#include <cstdio>

#include "command_line.h"
#include "input.h"
#include "utils.h"

namespace doctest {
template <typename T> struct StringMaker<std::vector<T>> {
	static String convert (const std::vector<T> & v) {
		using doctest::toString;
		doctest::String s = "Vec{";
		if (v.size () > 0) {
			s += toString (v[0]);
		}
		for (std::size_t i = 1; i < v.size (); ++i) {
			s += toString (", ") + toString (v[i]);
		}
		return s + toString ("}");
	}
};
template <> struct StringMaker<string_view> {
	static String convert (string_view sv) {
		using doctest::toString;
		return toString ("\"") + toString (to_string (sv)) + toString ("\"");
	}
};
} // namespace doctest

/******************************************************************************
 * Command line tests.
 */
TEST_SUITE ("command_line") {
	TEST_CASE ("usage") {
		// No actual checks, but useful to manually look at the formatting for different cases

		auto parser = CommandLineParser ();

		SUBCASE ("no arg, no opt") {}
		SUBCASE ("args, no opt") {
			parser.positional ("arg1", "Argument 1", [](string_view) {});
			parser.positional ("arg_______2", "Argument 2", [](string_view) {});
		}
		SUBCASE ("no args, opts") {
			parser.flag ({"h", "help"}, "Shows help", []() {});
			parser.option ({"f", "file"}, "file", "A file argument", [](string_view) {});
			parser.option2 ({"p"}, "first", "second", "A pair of values", [](string_view, string_view) {});
		}
		SUBCASE ("args and opts") {
			parser.flag ({"h", "help"}, "Shows help", []() {});
			parser.option ({"f", "file"}, "file", "A file argument", [](string_view) {});
			parser.option2 ({"p"}, "first", "second", "A pair of values", [](string_view, string_view) {});
			parser.positional ("arg1", "Argument 1", [](string_view) {});
			parser.positional ("arg_______2", "Argument 2", [](string_view) {});
		}

		fmt::print (stdout, "#################################\n");
		parser.usage (stdout, "test");
	}

	TEST_CASE ("construction_errors") {
		CommandLineParser parser;
		parser.flag ({"h", "help"}, "Shows help", []() {});

		// No name
		CHECK_THROWS_AS (parser.flag ({}, "blah", []() {}), CommandLineParser::Exception);
		// Empty name
		CHECK_THROWS_AS (parser.flag ({""}, "blah", []() {}), CommandLineParser::Exception);
		// Name collision
		CHECK_THROWS_AS (parser.flag ({"h"}, "blah", []() {}), CommandLineParser::Exception);
	}

	TEST_CASE ("parsing") {
		{
			const char * argv_data[] = {"prog_name", "-f", "--f"};
			auto argv = CommandLineView (sizeof (argv_data) / sizeof (*argv_data), argv_data);

			// f as a flag should match both
			auto flag_parser = CommandLineParser ();
			int flag_parser_f_seen = 0;
			flag_parser.flag ({"f"}, "", [&]() { flag_parser_f_seen++; });
			flag_parser.parse (argv);
			CHECK (flag_parser_f_seen == 2);

			// f as value opt eats the second --f
			auto value_parser = CommandLineParser ();
			value_parser.option ({"f"}, "", "", [](string_view value) { CHECK (value == "--f"); });
			value_parser.parse (argv);

			// Fails because args look like opts that are not defined
			auto nothing_parser = CommandLineParser ();
			CHECK_THROWS_AS (nothing_parser.parse (argv), CommandLineParser::Exception);

			// Success, no option declared
			auto arg_arg_parser = CommandLineParser ();
			arg_arg_parser.positional ("1", "", [](string_view v) { CHECK (v == "-f"); });
			arg_arg_parser.positional ("2", "", [](string_view v) { CHECK (v == "--f"); });
			arg_arg_parser.parse (argv);

			// Fails, with options parsing enabled '-f' is unknown.
			auto arg_opt_parser = CommandLineParser ();
			arg_opt_parser.flag ({"z"}, "", []() {});
			arg_opt_parser.positional ("1", "", [](string_view) {});
			arg_opt_parser.positional ("2", "", [](string_view) {});
			CHECK_THROWS_AS (arg_opt_parser.parse (argv), CommandLineParser::Exception);
		}
		{
			const char * argv_data[] = {"prog_name", "1a", "a2"};
			auto argv = CommandLineView (sizeof (argv_data) / sizeof (*argv_data), argv_data);

			// One unexpected arg
			auto arg_parser = CommandLineParser ();
			arg_parser.positional ("1", "", [](string_view) {});
			CHECK_THROWS_AS (arg_parser.parse (argv), CommandLineParser::Exception);

			// Eats both args
			auto arg_arg_parser = CommandLineParser ();
			arg_arg_parser.positional ("1", "", [](string_view v) { CHECK (v == "1a"); });
			arg_arg_parser.positional ("2", "", [](string_view v) { CHECK (v == "a2"); });
			arg_arg_parser.parse (argv);

			// Missing one arg
			auto arg_arg_arg_parser = CommandLineParser ();
			arg_arg_arg_parser.positional ("1", "", [](string_view) {});
			arg_arg_arg_parser.positional ("2", "", [](string_view) {});
			arg_arg_arg_parser.positional ("3", "", [](string_view) {});
			CHECK_THROWS_AS (arg_arg_arg_parser.parse (argv), CommandLineParser::Exception);
		}
		{
			// Value arg parsing
			const char * argv_data[] = {"prog_name", "-opt=value", "--opt", "value", "-opt", "value"};
			auto argv = CommandLineView (sizeof (argv_data) / sizeof (*argv_data), argv_data);

			auto parser = CommandLineParser ();
			parser.option ({"opt"}, "value", "desc", [](string_view v) { CHECK (v == "value"); });
			parser.parse (argv);
		}
		{
			// Value2 arg parsing
			const char * argv_data[] = {"prog_name", "-opt=value1", "value2", "--opt", "value1", "value2"};
			auto argv = CommandLineView (sizeof (argv_data) / sizeof (*argv_data), argv_data);

			auto parser = CommandLineParser ();
			parser.option2 ({"opt"}, "value", "value2", "desc", [](string_view v1, string_view v2) {
				CHECK (v1 == "value1");
				CHECK (v2 == "value2");
			});
			parser.parse (argv);
		}
	}
}

/******************************************************************************
 * File reading.
 * For testing, we read this exact file.
 */
constexpr auto & this_file_name = __FILE__;
extern const int this_file_nb_lines; // Initialized at the bottom of this file.

TEST_SUITE ("file_reading") {
	TEST_CASE ("open_file") {
		auto f = open_file (this_file_name, "r");
		CHECK (f != nullptr);

		CHECK_THROWS (open_file ("", "r"));
	}

	TEST_CASE ("LineByLineReader") {
		auto f = open_file (this_file_name, "r");
		LineByLineReader reader (f.get ());

		CHECK (reader.read_next_line () == true);
		CHECK (reader.current_line () == "// First line\n");
		CHECK (reader.current_line_number () == 0);

		CHECK (reader.read_next_line () == true);
		CHECK (reader.current_line () == "// Second line\n");
		CHECK (reader.current_line_number () == 1);

		while (reader.read_next_line ()) {
			// Get all lines
		}
		CHECK (reader.current_line_number () == this_file_nb_lines - 1);
	}
}

/******************************************************************************
 * Test utilities.
 */
TEST_SUITE ("utils") {
	TEST_CASE ("split") {
		CHECK (split (',', "") == std::vector<string_view>{""});
		CHECK (split (',', ",") == std::vector<string_view>{"", ""});
		CHECK (split (',', ",,") == std::vector<string_view>{"", "", ""});
		CHECK (split (',', "a,b,c") == std::vector<string_view>{"a", "b", "c"});
		CHECK (split (',', "a,b") == std::vector<string_view>{"a", "b"});
		CHECK (split (',', "a,") == std::vector<string_view>{"a", ""});
		CHECK (split (',', ",b") == std::vector<string_view>{"", "b"});
		CHECK (split (',', " ,b ") == std::vector<string_view>{" ", "b "});
	}

	TEST_CASE ("trim_ws") {
		CHECK (trim_ws_left ("") == "");
		CHECK (trim_ws_right ("") == "");
		CHECK (trim_ws ("") == "");

		CHECK (trim_ws_left (" \t\n") == "");
		CHECK (trim_ws_right (" \t\n") == "");
		CHECK (trim_ws (" \t\n") == "");

		CHECK (trim_ws_left (" a ") == "a ");
		CHECK (trim_ws_right (" a ") == " a");
		CHECK (trim_ws (" a ") == "a");
	}
}

// __LINE__ counts from 1. Thus the line number of the last line is the number of lines of the file.
const int this_file_nb_lines = __LINE__; // MUST BE THE LAST LINE !