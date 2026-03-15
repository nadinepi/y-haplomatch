#!/usr/bin/env bash

# start only the backend server
set -e
cd "$(dirname "$0")"

# find a free port starting from 5174 up to 5199
find_free_port() {
    local port=$1
    local max_port=$2
    while [ $port -le $max_port ]; do
        if ! ss -ltn | awk '{print $4}' | grep -q ":$port$"; then
            echo $port
            return 0
        fi
        port=$((port+1))
    done
    return 1
}

BACKEND_PORT=$(find_free_port 5174 5199)
if [ -z "$BACKEND_PORT" ]; then
    echo "No free port found in range 5174-5199. Please free up a port!"
    exit 1
fi
echo "Starting backend on port $BACKEND_PORT..."

# Prefer conda env "popgen" when available.
if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    if conda env list | awk '{print $1}' | grep -qx "popgen"; then
        conda activate popgen
        echo "Using conda environment: popgen"
    else
        echo "Conda found, but env 'popgen' does not exist."
        if [ -n "${VIRTUAL_ENV:-}" ]; then
            echo "Using current virtual environment: $VIRTUAL_ENV"
        else
            echo "Using system Python. If imports fail, install backend dependencies manually first."
        fi
    fi
else
    if [ -n "${VIRTUAL_ENV:-}" ]; then
        echo "Conda not found. Using current virtual environment: $VIRTUAL_ENV"
    else
        echo "Conda not found. Using system Python. If imports fail, install backend dependencies manually first."
    fi
fi

export PORT="$BACKEND_PORT"
python -u backend/app.py
