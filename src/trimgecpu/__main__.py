"""Allow ``python -m trimgecpu``."""

from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
