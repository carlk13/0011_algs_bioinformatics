import gzip
from collections.abc import Generator
from pathlib import Path

from tqdm import tqdm


def stream_lines(
    file_path: str | Path, limit: int | None = None, show_progress: bool = True
) -> Generator[str]:
    """
    Streams lines from a .txt or .gz file transparently with optional progress reporting.

    Args:
        file_path (str | Path): Path to the target text or gzip file.
        show_progress (bool): Whether to display a progress bar.

    Yields:
        str: Next line from the file (stripped of trailing newlines, or kept raw depending on preference).
    """
    path = Path(file_path)

    if path.suffix == ".gz":
        file_obj = gzip.open(path, "rt", encoding="utf-8")
    else:
        file_obj = open(path, "r", encoding="utf-8", buffering=1024 * 1024)

    line_count = 0

    with file_obj as f:
        stream = tqdm(f, desc="Streaming file", disable=not show_progress)

        for line in stream:
            if limit is not None and line_count >= limit:
                break

            # Strip trailing newlines (\r, \n)
            yield line.rstrip("\r\n")
            line_count += 1
