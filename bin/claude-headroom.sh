#!/usr/bin/env bash
set -euo pipefail

PORT="${HEADROOM_PORT:-8787}"
AGENT="${HEADROOM_AGENT:-${HEADROOM_CLIENT:-claude}}"

case "${1:-}" in
    claude|codex)
        AGENT="$1"
        shift
        ;;
esac

case "${AGENT}" in
    claude|codex)
        ;;
    *)
        echo "Unsupported Headroom agent: ${AGENT}" >&2
        echo "Usage: $0 [claude|codex] [args...]" >&2
        exit 2
        ;;
esac

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

# 会社のHTTPプロキシ設定がある場合の保険
export NO_PROXY="localhost,127.0.0.1"
export no_proxy="localhost,127.0.0.1"

case "${AGENT}" in
    claude)
        export ANTHROPIC_BASE_URL="http://127.0.0.1:${PORT}"
        echo "Using Headroom proxy for Claude: ${ANTHROPIC_BASE_URL}"
        exec claude "$@"
        ;;
    codex)
        export OPENAI_BASE_URL="http://127.0.0.1:${PORT}/v1"
        echo "Using Headroom proxy for Codex: ${OPENAI_BASE_URL}"
        exec codex "$@"
        ;;
esac
