#!/usr/bin/env python3
"""Smoke tests for the CLI adapter; Python is test-only, not a render dependency."""
import argparse
from pathlib import Path
import struct
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    args = parser.parse_args()
    binary = str(args.binary.resolve())
    checks = 0

    def run(*options, success=True):
        nonlocal checks
        result = subprocess.run([binary, *map(str, options)], capture_output=True, text=True, timeout=60)
        checks += 1
        if (result.returncode == 0) != success:
            raise AssertionError(f"Unexpected exit {result.returncode}: {options}\n{result.stdout}{result.stderr}")
        return result

    run("--help")
    listing = run("--list").stdout
    ids = [line.split("\t", 1)[0] for line in listing.splitlines()]
    assert len(ids) == len(set(ids)) == 32
    assert run("--list").stdout == listing
    with tempfile.TemporaryDirectory(prefix="fractal-render-test-") as temporary:
        directory = Path(temporary)
        for renderer in ids:
            image = directory / "image.bmp"
            serial = directory / "serial.bmp"
            restored = directory / "restored.bmp"
            view = directory / "state.view"
            common = ("--renderer", renderer, "--width", "24", "--height", "16", "--iterations", "32")
            run(*common, "--threads", "2", "--output", image, "--save-view", view)
            data = image.read_bytes()
            assert data[:2] == b"BM"
            assert struct.unpack_from("<I", data, 2)[0] == len(data) == 54 + 24 * 16 * 4
            assert struct.unpack_from("<ii", data, 18) == (24, -16)
            assert struct.unpack_from("<HH", data, 26) == (1, 32)
            run(*common, "--threads", "1", "--output", serial)
            assert serial.read_bytes() == data, renderer + ": worker counts differ"
            run("--view", view, "--width", "24", "--height", "16", "--output", restored)
            assert restored.read_bytes() == data, renderer + ": view round trip differs"

        for renderer in ("curvature", "extra_26", "extra_30", "extra_31"):
            for preset in range(1, 5):
                run("--renderer", renderer, "--preset", preset, "--width", 12, "--height", 8,
                    "--iterations", 32, "--samples", 2, "--output", directory / "preset.bmp")

        invalid = [
            ("--unknown",), ("--renderer",), ("--renderer", "invalid", "--output", directory / "bad.bmp"),
            ("--renderer", "extra_31"), ("--list", "--width", "20"),
            ("--renderer", "extra_31", "--view", directory / "state.view", "--output", directory / "bad.bmp"),
        ]
        for option, value in (("--width", "0"), ("--height", "8193"), ("--iterations", "31"),
                              ("--iterations", "1025"), ("--width", "2x"), ("--width", "-1"),
                              ("--threads", "0"), ("--samples", "5"), ("--preset", "0"),
                              ("--palette", "3"), ("--width", "9999999999999999999999999999")):
            invalid.append(("--renderer", "extra_31", "--output", directory / "bad.bmp", option, value))
        invalid += [
            ("--renderer", "extra_31", "--output", directory / "bad.png"),
            ("--renderer", "extra_31", "--width", "8192", "--height", "8192", "--output", directory / "bad.bmp"),
            ("--renderer", "extra_31", "--output", directory / "bad.bmp", "--width", "12", "--width", "12"),
            ("--renderer", "extra_31", "--output", directory / "bad.bmp", "--save-view", directory / "bad.bmp"),
            ("--view", directory / "missing.view", "--output", directory / "bad.bmp"),
        ]
        for options in invalid:
            run(*options, success=False)
            assert not (directory / "bad.bmp").exists()
    print(f"CLI_SMOKE=PASS RENDERERS=32 PRESET_RENDERS=16 INVALID_CASES={len(invalid)} INVOCATIONS={checks}")


if __name__ == "__main__":
    main()
