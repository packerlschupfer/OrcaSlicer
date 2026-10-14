#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "libslic3r/PrintConfig.hpp"

using namespace Slic3r;

namespace {

// An option's registration and the code that acts on it live in different files, so losing
// the registration while keeping the handler leaves a flag that parses as unknown and a
// feature with no way to reach it. Nothing else fails when that happens. Listing the fork's
// own CLI surface here is what turns it into a test failure.

const std::vector<std::string> fork_cli_actions = {
    "calibrate_type",
    "inspect_config",
    "list_calibrate_types",
    "paint",
    "render_paint",
};

const std::vector<std::string> fork_cli_params = {
    // Flow-rate calibration.
    "flow_blocks", "flow_brim", "flow_brim_gap", "flow_brim_width", "flow_height",
    "flow_pattern", "flow_range",
    // Tower generators.
    "pa_tower_end", "pa_tower_start", "pa_tower_step",
    "retraction_tower_end", "retraction_tower_start", "retraction_tower_step",
    "temp_tower_end", "temp_tower_start",
    "vfa_tower_end", "vfa_tower_start", "vfa_tower_step",
    "vs_tower_end", "vs_tower_start", "vs_tower_step",
    // Paint rendering.
    "render_clip_z_above", "render_clip_z_below", "render_height", "render_layer",
    "render_view", "render_width",
    // Z-offset calibration pattern.
    "zcal_fiducial_ids", "zcal_fiducials", "zcal_frame", "zcal_frame_layers",
    "zcal_frame_margin", "zcal_frame_width", "zcal_grid", "zcal_grid_over_zones",
    "zcal_grid_pitch", "zcal_scale_bar", "zcal_size", "zcal_struts",
    "zcal_zone_ids", "zcal_zone_labels", "zcal_zone_size",
    // Z-ladder.
    "zladder_band_height", "zladder_end", "zladder_fiducials", "zladder_height",
    "zladder_id_dots", "zladder_scale_bar", "zladder_start", "zladder_steps",
    "zladder_width",
};

// Upstream carries these now, having taken them from this fork. The fork's tooling still
// depends on them, and a rebase can drop them as readily as one of its own.
const std::vector<std::string> adopted_cli_actions = {
    "inspect_mesh", "inspect_paint", "strict",
};

const std::vector<std::string> adopted_cli_params = {
    "allow_rotations",
};

std::string dashed(const std::string &key)
{
    std::string out = key;
    boost::replace_all(out, "_", "-");
    return out;
}

// Every key the fork relies on, paired with the definition that has to carry it.
void require_registered(const ConfigDef &def, const std::vector<std::string> &keys)
{
    for (const std::string &key : keys) {
        CAPTURE(key);
        const ConfigOptionDef *opt = def.get(key);
        REQUIRE(opt != nullptr);

        // read_cli() matches on cli_args(), so this is the name the user actually types.
        const std::vector<std::string> args = opt->cli_args(key);
        const std::string registered = boost::join(args, ",");
        CAPTURE(registered);
        CHECK(std::find(args.begin(), args.end(), dashed(key)) != args.end());
    }
}

} // namespace

TEST_CASE("Every CLI action this fork adds is registered under its dashed name", "[CLIOptions]")
{
    require_registered(CLIActionsConfigDef(), fork_cli_actions);
}

TEST_CASE("Every CLI parameter this fork adds is registered under its dashed name", "[CLIOptions]")
{
    require_registered(CLIMiscConfigDef(), fork_cli_params);
}

TEST_CASE("CLI options upstream adopted from this fork are still registered", "[CLIOptions]")
{
    require_registered(CLIActionsConfigDef(), adopted_cli_actions);
    require_registered(CLIMiscConfigDef(), adopted_cli_params);
}

TEST_CASE("An action and a parameter definition do not claim the same CLI name", "[CLIOptions]")
{
    // --paint is an action and --paint-* were nearly parameters; a name served by two
    // definitions resolves to whichever read_cli() caches first, silently ignoring the other.
    const CLIActionsConfigDef   actions;
    const CLITransformConfigDef transform;
    const CLIMiscConfigDef      misc;

    std::vector<std::string> seen;
    for (const ConfigDef *def : {static_cast<const ConfigDef *>(&actions),
                                 static_cast<const ConfigDef *>(&transform),
                                 static_cast<const ConfigDef *>(&misc)})
        for (const auto &entry : def->options)
            for (const std::string &arg : entry.second.cli_args(entry.first))
                seen.push_back(arg);

    std::sort(seen.begin(), seen.end());
    const auto        it        = std::adjacent_find(seen.begin(), seen.end());
    const std::string duplicate = it == seen.end() ? "none" : *it;
    CAPTURE(duplicate);
    CHECK(duplicate == "none");
}
