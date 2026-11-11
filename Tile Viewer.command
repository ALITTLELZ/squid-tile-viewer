#!/bin/bash
# Double-click to start Tile Viewer (macOS). 双击启动。
cd "$(dirname "$0")"

UV="$(command -v uv)"
for c in "$HOME/.local/bin/uv" "$HOME/.cargo/bin/uv" /opt/homebrew/bin/uv /usr/local/bin/uv; do
  [ -z "$UV" ] && [ -x "$c" ] && UV="$c"
done

if [ -z "$UV" ]; then
  echo "Tile Viewer uses 'uv' (a small Python installer from astral.sh) to set up its dependencies."
  echo "Tile Viewer 需要 uv（astral.sh 出品的 Python 环境管理工具）来安装依赖。"
  read -r -p "Install uv now? / 现在安装 uv？ [y/N] " ans
  if [[ "$ans" =~ ^[Yy] ]]; then
    curl -LsSf https://astral.sh/uv/install.sh | sh || { read -r -p "uv install failed. Press Enter to close."; exit 1; }
    UV="$HOME/.local/bin/uv"
  else
    exit 0
  fi
fi

echo "Starting Tile Viewer ... (first launch downloads ~300 MB of packages and may take a few minutes)"
echo "正在启动 ...（第一次启动需要下载约 300 MB 的依赖，可能要几分钟）"
"$UV" run --python 3.12 tile_viewer.py "$@"
status=$?
if [ $status -ne 0 ]; then
  read -r -p "Tile Viewer exited with an error (see above). Press Enter to close."
fi
