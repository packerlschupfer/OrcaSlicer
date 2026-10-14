#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include "libslic3r/PrintConfig.hpp"
#include "slic3r/Utils/CalibCLI.hpp"

using namespace Slic3r;

namespace {

// Each generator is reachable only if three lists agree: the enum, the string parser, and
// the --calibrate-type help that tells the user the name. They sit in two libraries, so
// nothing but a test spanning both can hold them together -- flow-yolo-coarse was accepted
// by the parser and named in neither user-facing list.
const std::vector<std::string> calibrate_type_names = {
    "flow-yolo-recommended",
    "flow-yolo-perfectionist",
    "flow-yolo-coarse",
    "flow-pass1",
    "flow-pass2",
    "z-offset-pattern",
    "temp-tower",
    "vol-speed-tower",
    "pa-tower",
    "retraction-tower",
    "vfa-tower",
    "z-ladder-banded",
    "z-ladder-ramp",
};

} // namespace

TEST_CASE("Every calibration generator name round-trips through the CLI parser", "[CalibCLI]")
{
    for (const std::string &name : calibrate_type_names) {
        CAPTURE(name);
        const CLICalibType type = cli_calib_type_from_string(name);
        REQUIRE(type != CLICalibType::NoCalib);
        CHECK(std::string(cli_calib_type_to_string(type)) == name);
    }
}

TEST_CASE("An unknown calibration type is rejected rather than defaulted", "[CalibCLI]")
{
    CHECK(cli_calib_type_from_string("") == CLICalibType::NoCalib);
    CHECK(cli_calib_type_from_string("flow-yolo") == CLICalibType::NoCalib);
    CHECK(cli_calib_type_from_string("Flow-Pass1") == CLICalibType::NoCalib);
}

TEST_CASE("No calibration type is left without a name", "[CalibCLI]")
{
    // Walks the enum by value, so a type dropped from the parser or the name table shows up
    // as "none" here even though its enumerator still exists.
    for (int i = 1; i <= int(calibrate_type_names.size()); ++i) {
        const CLICalibType type = CLICalibType(i);
        CAPTURE(i);
        const std::string name = cli_calib_type_to_string(type);
        CHECK(name != "none");
        CHECK(std::find(calibrate_type_names.begin(), calibrate_type_names.end(), name) !=
              calibrate_type_names.end());
    }
}

TEST_CASE("The calibrate-type help names every generator the parser accepts", "[CalibCLI]")
{
    const ConfigOptionDef *def = CLIActionsConfigDef().get("calibrate_type");
    REQUIRE(def != nullptr);

    for (const std::string &name : calibrate_type_names) {
        CAPTURE(name);
        CHECK(def->tooltip.find(name) != std::string::npos);
    }
}
