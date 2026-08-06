#!/usr/bin/env bash
set -euo pipefail

# Compile the code into a build folder
cmake -S . -B build \
  --log-level=NOTICE \
  -Wno-deprecated \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.12

cmake --build build 2>&1 | tee build.log

# Run our C executable

setsid ./build/simulation &
CPP_PID0=$! # This starts the executable with ID CPP_PID0
# If we have multiple executables this allows us to keep track of them all


cleanup() {
    kill "$CPP_PID0" 2>/dev/null || true
}

# Wait until all processes we care about exit
wait -n "$CPP_PID0"

# Then shut down the other ones
cleanup
wait || true