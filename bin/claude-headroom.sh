#!/usr/bin/env bash
set -euo pipefail

PORT="${HEADROOM_PORT:-8787}"

# Headroom起動済みなら再利用
if ! curl -sf "http://127.0.0.1:${PORT}/health" >/dev/null 2>&1; then
    echo "Starting Headroom on port ${PORT}..."
    headroom proxy --port "${PORT}" >/tmp/headroom.log 2>&1 &
    HEADROOM_PID=$!

    # 起動待ち
    for i in {1..20}; do
        if curl -sf "http://127.0.0.1:${PORT}/health" >/dev/null 2>&1; then
            break
        fi
        sleep 1
    done
fi

export ANTHROPIC_BASE_URL="http://127.0.0.1:${PORT}"

# 会社のHTTPプロキシ設定がある場合の保険
export NO_PROXY="localhost,127.0.0.1"
export no_proxy="localhost,127.0.0.1"

echo "Using Headroom proxy: ${ANTHROPIC_BASE_URL}"

exec claude "$@"