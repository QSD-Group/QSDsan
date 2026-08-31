from __future__ import annotations

import re
from pathlib import Path


DEPENDENCIES = ("thermosteam", "biosteam")
PIN_RE = re.compile(r"^(thermosteam|biosteam)==([0-9]+)\.([0-9]+)\.([0-9]+)\s*$")


def parse_release_pins(requirements_path: Path) -> dict[str, tuple[int, int, int]]:
    pins = {}
    for line in requirements_path.read_text().splitlines():
        match = PIN_RE.match(line.strip())
        if match:
            name, major, minor, patch = match.groups()
            pins[name] = (int(major), int(minor), int(patch))
    missing = sorted(set(DEPENDENCIES) - set(pins))
    if missing:
        raise ValueError(f"Missing release pins for: {', '.join(missing)}")
    return pins


def sync_lower_bounds(pyproject_path: Path, requirements_path: Path) -> bool:
    pins = parse_release_pins(requirements_path)
    text = pyproject_path.read_text()
    changed = False

    for name, release_version in pins.items():
        pattern = re.compile(
            rf'("{re.escape(name)}>=)([0-9]+)\.([0-9]+)\.([0-9]+)(")'
        )

        def replace(match: re.Match[str]) -> str:
            nonlocal changed
            current = tuple(int(match.group(i)) for i in range(2, 5))
            if current[:2] == release_version[:2]:
                return match.group(0)
            changed = True
            version = ".".join(str(part) for part in release_version)
            return f"{match.group(1)}{version}{match.group(5)}"

        text, count = pattern.subn(replace, text, count=1)
        if count == 0:
            raise ValueError(f"Could not find lower bound for {name!r}")

    if changed:
        pyproject_path.write_text(text)
    return changed


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[2]
    sync_lower_bounds(root / "pyproject.toml", root / ".github" / "requirements-release.txt")
