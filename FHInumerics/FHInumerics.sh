#!/usr/bin/env bash
# Start FHInumerics (Linux/Mac). On Windows, double-click FHInumerics.bat.
exec Rscript "$(dirname "${BASH_SOURCE[0]}")/launcher/launch.R"
