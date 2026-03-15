#!/usr/bin/env bash
# Start only the frontend dev server
set -e
cd "$(dirname "$0")/frontend"

# find a free port
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

FRONTEND_PORT=$(find_free_port 5173 5190)
if [ -z "$FRONTEND_PORT" ]; then
    echo "No free port found in range 5173-5190. Exiting."
    exit 1
fi
echo "Starting frontend on port $FRONTEND_PORT..."
npm run dev -- --host 0.0.0.0 --port "$FRONTEND_PORT"
