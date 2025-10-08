#!/usr/bin/env bash

set -euo pipefail

# Absolute project root
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

PYTHON_BIN="python3"
VENV_DIR="$PROJECT_ROOT/.venv"

echo "--- Setting up virtual environment ---"
if [ ! -d "$VENV_DIR" ]; then
  "$PYTHON_BIN" -m venv "$VENV_DIR"
fi
source "$VENV_DIR/bin/activate"

echo "--- Upgrading pip/setuptools/wheel ---"
pip install --upgrade pip setuptools wheel

echo "--- Installing build tools and PyInstaller ---"
pip install pyinstaller cmake

echo "--- Installing pybind11 (required for C++ extension build) ---"
pip install "pybind11>=2.6"
echo "--- Installing build helpers needed by libbdsg (no isolation) ---"
pip install setuptools_scm setuptools_git_ls_files

echo "--- Injecting temporary cmake shim to add policy flag on configure ---"
REAL_CMAKE_PATH="$(command -v cmake)"
SHIM_DIR="$VENV_DIR/cmake_shim"
mkdir -p "$SHIM_DIR"
cat > "$SHIM_DIR/cmake" << 'EOF'
#!/usr/bin/env bash
set -euo pipefail
# cmake shim: append policy minimum only on configure (no --build)
if echo " $* " | grep -q " --build "; then
  exec "$REAL_CMAKE_PATH" "$@"
fi
if echo " $* " | grep -q "CMAKE_POLICY_VERSION_MINIMUM"; then
  exec "$REAL_CMAKE_PATH" "$@"
fi
exec "$REAL_CMAKE_PATH" "$@" -DCMAKE_POLICY_VERSION_MINIMUM=3.5
EOF
chmod +x "$SHIM_DIR/cmake"
export REAL_CMAKE_PATH
export PATH="$SHIM_DIR:$PATH"

echo "--- Installing Python runtime dependencies ---"
if [ -f "$PROJECT_ROOT/requirements.txt" ]; then
  pip install -r "$PROJECT_ROOT/requirements.txt"
fi

echo "--- Building sdust (third_party) if needed ---"
make sdust

echo "--- Building and installing libbdsg ---"
(
  set -e
  cd "$PROJECT_ROOT/libbdsg"
  mkdir -p build
  cd build
  # If libbdsg's build does not pass the policy flag, inject it here without modifying libbdsg
  cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=3.5
  make -j"$(command -v nproc >/dev/null 2>&1 && nproc || echo 8)"
  cd ..
  PIP_NO_BUILD_ISOLATION=1 pip install --no-build-isolation .
)

echo "--- Building local C++ extension in-place (avoid PEP517 isolation) ---"
PIP_NO_BUILD_ISOLATION=1 "$PYTHON_BIN" setup.py build_ext --inplace

echo "--- Building one-file executable with PyInstaller using spec ---"
VERSION="0.1.0"
pyinstaller "$PROJECT_ROOT/vg-anchors-$VERSION.spec" --noconfirm

OUT_DIR="$PROJECT_ROOT/dist"
EXE_PATH="$OUT_DIR/vg-anchors-$VERSION"

if [ -f "$EXE_PATH" ]; then
  echo "--- Executable created at: $EXE_PATH ---"
else
  echo "ERROR: Executable not found. Check PyInstaller output in dist/." >&2
  exit 1
fi

echo "--- Done ---"


