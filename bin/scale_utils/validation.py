"""Validation functions"""

import sys


def validateName(name: str, display_name: str = "Name", other_chars: str = "-."):
    """
    Check name for invalid characters
    Print error and exit for invalid names

    Args:
        name: str to check
        display_name: name to display in error message
        other_chars: what other characters are allowed beside alpha-numeric
    """
    for n in name:
        if not (n.isalnum() or n in other_chars):
            print(
                f"{display_name} should only contain [a-z],[A-Z],[0-9], or {', '.join([f'[{char}]' for char in other_chars])}: '{name}'",
                file=sys.stderr,
            )
            sys.exit(1)
    if not name[0].isalpha():
        print(f"Name should start with a letter: {name}", file=sys.stderr)
        sys.exit(1)
