# Command-Line Input Paths: High Level Design

## Why it exists

`GUI_App::init_app_config()` changes the process working directory to
`<data_dir>/log`, early in startup and on every platform (`_chdir` on Windows). The files
named on the command line are not opened there and then: the GUI opens them much later,
from `post_init()`, once a window and a canvas exist.

Anything relative therefore resolves against the log directory rather than the directory
the user typed the command in. The file is not found, and the failure surfaces far from
its cause.

`CLI::setup()` closes that gap by resolving the input files immediately after
`read_cli()`, while the process is still where the user started it. By the time the
working directory moves, every input file is already an absolute path.

## What resolution does

`resolve_cli_input_path()` returns its argument unchanged when it is empty, already
absolute, or a supported open protocol (the `orcaslicer://open`, `prusaslicer://open`,
`bambustudio://open`, `bambustudioopen://` and `cura://open` URLs the GUI can be handed). Otherwise it makes the path
absolute against the current directory and normalises it lexically. A path that cannot be
resolved is returned as it came in, so a malformed argument still reaches the code that
reports it properly.

## Why the failure is hard to recognise

Nothing reports a path problem. The archive simply does not open, so the importer returns
false and the user sees two dialogs in a row:

    Loading of a model file failed.
    The file does not contain any geometry data.

Both blame the model. The file is intact, and opening the same file by an absolute path —
or by double-clicking it, because the desktop entry passes `%F`, which is absolute — works.
The log line that identifies it is `Unable to open the file<path>`, with the raw relative
path still attached; that string is the reliable way to tell this apart from a genuinely
damaged 3MF.

The CLI itself never had the bug, because it does not change directory. A relative path
that works with `--inspect-mesh` and fails when the GUI opens the same file is this, and
that asymmetry is the quickest way to recognise it.

## Scope

Only `input_files` is resolved. `GUI_InitParams::load_configs`, which carries
`--load-settings`, is not — the GUI has no live reader for it at all, since the only call
that consumed it sits inside a commented-out block in `post_init()`. A preset named on
the command line is therefore ignored in GUI mode regardless of how its path is spelled.
