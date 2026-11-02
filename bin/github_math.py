#!/usr/bin/env python3
r"""
Convert LaTeX math delimiters in Markdown to GitHub-style dollar delimiters.

REQUIREMENTS
    Python 3.8 or newer. No packages need to be installed.

WHAT IT DOES
    Display math:  \[ ... \]  becomes  $$ on separate lines around the math.
    Inline math:   \( ... \)  becomes  $...$.
    Fenced code blocks and inline backtick code are left unchanged.
    Existing dollar-delimited math is not deliberately rewritten.
    Unmatched delimiters are left unchanged.

USAGE (run these commands in your terminal)
    Convert one file, keeping the original:
        python3 github_math.py study_guide.md
        # Writes study_guide.github.md

    Choose the output filename:
        python3 github_math.py study_guide.md -o converted.md

    Convert several files:
        python3 github_math.py chapter1.md chapter2.md

    Convert the .md files directly inside a folder:
        python3 github_math.py ./guides

    Include subfolders:
        python3 github_math.py ./guides --recursive


    Update originals inplace while saving a backup for each changed file:
        python3 github_math.py ./guides --recursive --in-place
        # Saves the original study_guide.md as study_guide.md.bak
        # Refuses to overwrite an existing backup.

    Show all options:
        python3 github_math.py --help

NOTES
    On Windows, use "py" instead of "python3" if needed.
    Quote paths that contain spaces: "My Study Guide.md".
    Files must be UTF-8. Existing line endings and a UTF-8 BOM are preserved.
    Existing output files are never overwritten.
    Directory scans skip generated *.github.md files.
    This is a lightweight Markdown converter, not a complete Markdown parser.
    It protects ordinary backtick/tilde fences (up to three leading spaces)
    and inline backtick code, but not indented or blockquote/list-nested code.
    It changes delimiters only; it does not repair invalid LaTeX expressions.
"""

import argparse
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple


FENCE_OPEN = re.compile(r"^ {0,3}(`{3,}|~{3,})([^\r\n]*)")
BACKTICKS = re.compile(r"`+")


def is_escaped(text: str, index: int) -> bool:
    """Whether an odd number of backslashes precedes this character."""
    count = 0
    index -= 1
    while index >= 0 and text[index] == "\\":
        count += 1
        index -= 1
    return count % 2 == 1


def find_closing(text: str, delimiter: str, start: int) -> int:
    position = text.find(delimiter, start)
    while position != -1 and is_escaped(text, position):
        position = text.find(delimiter, position + len(delimiter))
    return position


def convert_prose(text: str, newline: str) -> Tuple[str, int]:
    """Convert a region outside fenced code, honoring inline code spans."""
    result: List[str] = []
    changes = 0
    index = 0
    while index < len(text):
        if text[index] == "`" and not is_escaped(text, index):
            opening = BACKTICKS.match(text, index)
            assert opening is not None
            width = len(opening.group())
            closing = next(
                (run for run in BACKTICKS.finditer(text, opening.end())
                 if len(run.group()) == width),
                None,
            )
            if closing is not None:
                result.append(text[index:closing.end()])
                index = closing.end()
                continue
            # An unmatched backtick run is ordinary text, not a code span.
            result.append(opening.group())
            index = opening.end()
            continue

        if text[index:index + 2] in (r"\[", r"\(") and not is_escaped(text, index):
            display = text[index + 1] == "["
            closing = find_closing(text, r"\]" if display else r"\)", index + 2)
            if closing != -1:
                body = text[index + 2:closing]
                end = closing + 2
                if display:
                    # Put display delimiters on their own lines, even when the
                    # original expression was written as \[x\] within a line.
                    before = newline if index and text[index - 1] not in "\r\n" else ""
                    after = newline if end < len(text) and text[end] not in "\r\n" else ""
                    result.append(
                        before + "$$" + newline + body.strip("\r\n")
                        + newline + "$$" + after
                    )
                else:
                    result.append("$" + body + "$")
                changes += 1
                index = end
                continue

        result.append(text[index])
        index += 1
    return "".join(result), changes


def convert_markdown(text: str) -> Tuple[str, int]:
    """Return converted text and the number of converted math expressions."""
    newline = "\r\n" if "\r\n" in text else "\n"
    output: List[str] = []
    prose: List[str] = []
    fence_char: Optional[str] = None
    fence_width = 0
    changes = 0

    def flush_prose() -> None:
        nonlocal changes
        converted, count = convert_prose("".join(prose), newline)
        output.append(converted)
        changes += count
        prose.clear()

    for line in text.splitlines(keepends=True):
        if fence_char is not None:
            output.append(line)
            closing_pattern = (
                r" {0,3}" + re.escape(fence_char)
                + "{" + str(fence_width) + r",}[ \t]*"
            )
            if re.fullmatch(closing_pattern, line.rstrip("\r\n")):
                fence_char = None
            continue

        opening = FENCE_OPEN.match(line)
        # Backtick fence info strings cannot themselves contain backticks.
        if opening and not (
            opening.group(1)[0] == "`" and "`" in opening.group(2)
        ):
            flush_prose()
            fence_char = opening.group(1)[0]
            fence_width = len(opening.group(1))
            output.append(line)
        else:
            prose.append(line)

    flush_prose()
    return "".join(output), changes


def collect_files(paths: List[str], recursive: bool) -> List[Path]:
    files: List[Path] = []
    seen = set()
    for value in paths:
        path = Path(value).expanduser()
        if path.is_dir():
            candidates = sorted(path.rglob("*") if recursive else path.iterdir())
            candidates = [
                item for item in candidates
                if item.is_file() and item.suffix.lower() == ".md"
                and not item.name.lower().endswith(".github.md")
            ]
        elif path.is_file() and path.suffix.lower() == ".md":
            candidates = [path]
        else:
            raise ValueError("Not a Markdown file or directory: {}".format(path))
        for candidate in candidates:
            if candidate.is_symlink():
                raise ValueError("Refusing to modify or convert a symlink: {}".format(candidate))
            identity = candidate.resolve()
            if identity not in seen:
                files.append(candidate)
                seen.add(identity)
    return files


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert Markdown math from \\[...\\] / \\(...\\) to dollar delimiters.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="By default, guide.md becomes guide.github.md; the original is preserved.",
    )
    parser.add_argument("paths", nargs="+", help="Markdown files or directories")
    parser.add_argument("-r", "--recursive", action="store_true", help="Include subfolders")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("-o", "--output", type=Path, help="Output path; requires exactly one input file")
    mode.add_argument("-i", "--in-place", action="store_true", help="Edit originals; save .bak backups")
    args = parser.parse_args()

    try:
        files = collect_files(args.paths, args.recursive)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if not files:
        parser.error("No Markdown files found.")
    if args.output is not None and len(files) != 1:
        parser.error("--output requires exactly one input file.")

    failures = 0
    for path in files:
        try:
            original = path.read_bytes()
            bom = b"\xef\xbb\xbf" if original.startswith(b"\xef\xbb\xbf") else b""
            text = original[len(bom):].decode("utf-8")
            converted, count = convert_markdown(text)
            if count == 0:
                print("Unchanged: {} (no convertible math found)".format(path))
                continue
            encoded = bom + converted.encode("utf-8")
            if args.in_place:
                backup = path.with_name(path.name + ".bak")
                with backup.open("xb") as handle:
                    handle.write(original)
                path.write_bytes(encoded)
                print("Updated: {} ({} expressions; backup: {})".format(path, count, backup))
            else:
                destination = (
                    args.output.expanduser() if args.output is not None
                    else path.with_name(path.stem + ".github" + path.suffix)
                )
                with destination.open("xb") as handle:
                    handle.write(encoded)
                print("Created: {} ({} expressions)".format(destination, count))
        except (OSError, UnicodeError) as error:
            failures += 1
            print("ERROR: {}: {}".format(path, error), file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
