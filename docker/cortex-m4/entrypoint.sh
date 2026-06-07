#!/bin/bash
# cfx ARM Cortex-M4 container entrypoint.
#
# All build/test logic lives in scripts/build-cortex-m4.sh; this just forwards the command.
exec /cfx/scripts/build-cortex-m4.sh "$@"
