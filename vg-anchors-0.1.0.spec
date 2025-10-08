# -*- mode: python ; coding: utf-8 -*-

from PyInstaller.utils.hooks import collect_all

# Collect bdsg package (pure files, binaries, hidden imports)
bdsg_datas, bdsg_binaries, bdsg_hiddenimports = collect_all('bdsg')

a = Analysis(
    ['assembler/cli.py'],
    pathex=[],
    binaries=[
        ('bin/sdust', 'bin'),
        ('shasta2/shasta2.so', '.')
    ] + bdsg_binaries,
    datas=[('libbdsg/lib', 'lib'), ('config.ini', '.')] + bdsg_datas,
    hiddenimports=['shasta2'] + bdsg_hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=['pyi_rth_vg_anchors_libpath.py'],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='vg-anchors-0.1.0',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
