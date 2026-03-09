#!/usr/bin/env python3

"""
A Python implementation of thisutils::log_message for use in the CLI and as a standalone function.
References:
 - https://github.com/mengxu98/thisutils/blob/main/R/log_message.R
 - https://github.com/mengxu98/thisutils/blob/main/scripts/log_message.py
"""

import argparse
import datetime as dt
import inspect
import os
import re
import sys
from typing import List, Optional, Sequence


class LogMessageError(RuntimeError):
    pass


def _is_hex_color(color: str) -> bool:
    return bool(re.fullmatch(r"#([0-9a-fA-F]{3}|[0-9a-fA-F]{6})", color))


def check_color(color: Optional[str]) -> bool:
    if color is None:
        return True
    if not isinstance(color, str):
        return False
    return color in {
        "black",
        "red",
        "green",
        "yellow",
        "blue",
        "magenta",
        "cyan",
        "white",
        "grey",
        "silver",
        "none",
        "br_black",
        "br_red",
        "br_green",
        "br_yellow",
        "br_blue",
        "br_magenta",
        "br_cyan",
        "br_white",
    } or _is_hex_color(color)


def hex_to_rgb(hex_value: str) -> tuple[int, int, int]:
    hex_value = hex_value.lstrip("#")
    if len(hex_value) == 3:
        hex_value = "".join(ch * 2 for ch in hex_value)
    return int(hex_value[0:2], 16), int(hex_value[2:4], 16), int(hex_value[4:6], 16)


def make_color_code(color: str, bg: bool = False) -> Optional[str]:
    if color is None:
        return None
    if _is_hex_color(color):
        r, g, b = hex_to_rgb(color)
        return f"{'48' if bg else '38'};2;{r};{g};{b}"

    if bg:
        color_codes = {
            "black": "40",
            "red": "41",
            "green": "42",
            "yellow": "43",
            "blue": "44",
            "magenta": "45",
            "cyan": "46",
            "white": "47",
            "none": "49",
            "br_black": "100",
            "br_red": "101",
            "br_green": "102",
            "br_yellow": "103",
            "br_blue": "104",
            "br_magenta": "105",
            "br_cyan": "106",
            "br_white": "107",
        }
    else:
        color_codes = {
            "black": "30",
            "red": "31",
            "green": "32",
            "yellow": "33",
            "blue": "34",
            "magenta": "35",
            "cyan": "36",
            "white": "37",
            "grey": "90",
            "silver": "37",
            "none": "39",
            "br_black": "90",
            "br_red": "91",
            "br_green": "92",
            "br_yellow": "93",
            "br_blue": "94",
            "br_magenta": "95",
            "br_cyan": "96",
            "br_white": "97",
        }

    return color_codes.get(color)


def style_formatting(
    msg: str,
    text_color: Optional[str],
    back_color: Optional[str],
    text_style: Optional[Sequence[str]],
) -> str:
    if text_color is None and back_color is None and not text_style:
        return msg

    codes: List[str] = []
    if text_color is not None:
        code = make_color_code(text_color, bg=False)
        if code is not None:
            codes.append(code)
    if back_color is not None:
        code = make_color_code(back_color, bg=True)
        if code is not None:
            codes.append(code)
    if text_style:
        style_codes = {
            "bold": "1",
            "italic": "3",
            "underline": "4",
            "strikethrough": "9",
            "dim": "2",
            "inverse": "7",
        }
        for style in text_style:
            style_code = style_codes.get(style)
            if style_code:
                codes.append(style_code)

    if not codes:
        return msg
    return f"\033[{';'.join(codes)}m{msg}\033[0m"


def capitalize(text: str, force_tolower: bool = False) -> str:
    if text is None:
        return text
    if len(text) == 0:
        return text

    if force_tolower:
        return text[:1].upper() + text[1:].lower()

    first_word = re.split(r"[\s-]+", text)[0] if text else ""
    if first_word == first_word.lower() and text:
        return text[:1].upper() + text[1:]
    return text


def build_message(*args: object) -> str:
    if len(args) == 0:
        return ""

    processed: List[str] = []
    for arg in args:
        if isinstance(arg, str):
            processed.append(arg)
        elif isinstance(arg, (list, tuple)):
            processed.append("".join(str(x) for x in arg))
        else:
            processed.append(str(arg))

    msg = "".join(processed)
    return capitalize(msg)


def apply_inline_format(text: str, format_type: str, evaluated: bool = False) -> str:
    inline_formats = {
        ".pkg": {"kind": "blue_plain"},
        ".code": {"kind": "backtick"},
        ".val": {"kind": "blue_quote_or_eval"},
        ".arg": {"kind": "backtick"},
        ".fun": {"kind": "function"},
        ".fn": {"kind": "function"},
        ".file": {"kind": "blue_single_quote"},
        ".path": {"kind": "blue_plain"},
        ".field": {"kind": "plain"},
        ".emph": {"kind": "italic"},
        ".strong": {"kind": "bold"},
        ".run": {"kind": "backtick"},
        ".str": {"kind": "double_quote"},
        ".help": {"kind": "help"},
        ".topic": {"kind": "topic"},
        ".vignette": {"kind": "vignette"},
        ".url": {"kind": "blue_italic_url"},
        ".email": {"kind": "blue_single_quote"},
        ".kbd": {"kind": "blue_kbd"},
        ".key": {"kind": "blue_kbd"},
        ".envvar": {"kind": "backtick"},
        ".dt": {"kind": "dt"},
        ".cls": {"kind": "class"},
        ".href": {"kind": "href_url_blue_italic"},
    }
    fmt = inline_formats.get(format_type)
    if fmt is None:
        return text

    kind = fmt.get("kind", "plain")
    if kind == "plain":
        return text
    if kind == "blue_plain":
        return style_formatting(text, "blue", None, None)
    if kind == "backtick":
        return f"`{text}`"
    if kind == "double_quote":
        return f'"{text}"'
    if kind == "single_quote":
        return f"'{text}'"
    if kind == "blue_single_quote":
        return style_formatting(f"'{text}'", "blue", None, None)
    if kind == "function":
        return f"`{text}()`"
    if kind == "help":
        return f"`help({text})`"
    if kind == "topic":
        return f"`help({text})`"
    if kind == "vignette":
        return f"`vignette({text})`"
    if kind == "url":
        return f"<{text}>"
    if kind == "blue_italic_url":
        return style_formatting(f"<{text}>", "blue", None, ["italic"])
    if kind == "kbd":
        return f"[{text}]"
    if kind == "blue_kbd":
        return style_formatting(f"[{text}]", "blue", None, None)
    if kind == "dt":
        return f"{text}: "
    if kind == "class":
        return f"<{text}>"
    if kind == "href_url_blue_italic":
        return text
    if kind == "italic":
        return style_formatting(text, None, None, ["italic"])
    if kind == "bold":
        return style_formatting(text, None, None, ["bold"])
    if kind == "quote_or_eval":
        return text if evaluated else f'"{text}"'
    if kind == "blue_quote_or_eval":
        content = text if evaluated else f'"{text}"'
        return style_formatting(content, "blue", None, None)
    return text


def parse_inline_expressions(message: str, caller_frame=None) -> str:
    max_iterations = 15
    iteration = 0
    allowed_tags = {
        "pkg",
        "code",
        "val",
        "arg",
        "fun",
        "fn",
        "file",
        "path",
        "field",
        "emph",
        "strong",
        "run",
        "str",
        "help",
        "topic",
        "vignette",
        "url",
        "email",
        "kbd",
        "key",
        "envvar",
        "dt",
        "cls",
        "href",
    }

    def eval_expression(expr: str) -> str:
        expr = expr.strip()
        if not expr:
            return ""
        if caller_frame is None:
            return expr
        try:
            value = eval(expr, caller_frame.f_globals, caller_frame.f_locals)
            return str(value)
        except Exception:
            return expr

    while iteration < max_iterations:
        iteration += 1

        href_match = re.search(r"\{\.href\s+\[([^\]]+)\]\(([^\)]+)\)\}", message)
        if href_match:
            title = href_match.group(1)
            url = href_match.group(2)
            styled_url = style_formatting(f"<{url}>", "blue", None, ["italic"])
            replacement = f"{title} ({styled_url})"
            message = (
                message[: href_match.start()]
                + replacement
                + message[href_match.end() :]
            )
            continue

        format_match = re.search(
            r"\{\.?([A-Za-z_][A-Za-z0-9_]*)(\s+)(\{[^{}]*\}|[^{}]+)\}", message
        )
        if format_match:
            tag = format_match.group(1)
            separator = format_match.group(2)
            content = format_match.group(3)

            if tag in allowed_tags:
                if content.startswith("{") and content.endswith("}"):
                    evaluated = eval_expression(content[1:-1])
                    was_evaluated = True
                else:
                    evaluated = content
                    was_evaluated = False

                if len(separator) > 1:
                    evaluated = separator[1:] + evaluated

                formatted = apply_inline_format(
                    evaluated, f".{tag}", evaluated=was_evaluated
                )
                message = (
                    message[: format_match.start()]
                    + formatted
                    + message[format_match.end() :]
                )
                continue

        bare_match = re.search(r"\{([^{}]+)\}", message)
        if bare_match:
            inner = bare_match.group(1).strip()
            if re.match(r"^\.?[A-Za-z_][A-Za-z0-9_]*\s+", inner):
                break
            evaluated = eval_expression(inner)
            message = (
                message[: bare_match.start()] + evaluated + message[bare_match.end() :]
            )
            continue

        break

    return message


def parse_bool(value: object) -> Optional[bool]:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "t", "1", "yes", "y"}:
            return True
        if lowered in {"false", "f", "0", "no", "n"}:
            return False
    return None


def get_verbose(verbose: Optional[bool] = None) -> bool:
    env_raw = os.getenv("LOG_MESSAGE_VERBOSE")

    if env_raw is None:
        if verbose is None:
            return True
        parsed_local = parse_bool(verbose)
        if parsed_local is None:
            print(
                "WARNING: verbose is not a logical value, set to TRUE", file=sys.stderr
            )
            return True
        return parsed_local

    parsed_env = parse_bool(env_raw)
    if parsed_env is None:
        print(
            "WARNING: LOG_MESSAGE_VERBOSE is not a logical value, treated as NULL",
            file=sys.stderr,
        )
        print("WARNING: or use LOG_MESSAGE_VERBOSE=true/false", file=sys.stderr)
        return False

    return parsed_env


def get_indent_part(symbol: str, level: int) -> str:
    if symbol != "  ":
        return f"{symbol * level} "
    if level > 1:
        return "  " * (level - 1)
    return ""


def validate_params(
    level: int,
    symbol: str,
    text_color: Optional[str],
    back_color: Optional[str],
    text_style: Optional[Sequence[str]],
    timestamp_style: bool,
) -> None:
    if not isinstance(level, int) or level < 1:
        raise LogMessageError("level must be a positive integer")

    if not isinstance(symbol, str):
        raise LogMessageError("symbol must be a single character string")

    if text_color is not None and not check_color(text_color):
        raise LogMessageError(
            "text_color must be a valid color name or hexadecimal color code"
        )

    if back_color is not None and not check_color(back_color):
        raise LogMessageError(
            "back_color must be a valid color name or hexadecimal color code"
        )

    if text_color is not None and back_color is not None and text_color == back_color:
        raise LogMessageError("text_color and back_color cannot be the same color")

    if text_style is not None and any(
        style not in {"bold", "italic", "underline", "strikethrough", "dim", "inverse"}
        for style in text_style
    ):
        raise LogMessageError(
            "text_style must be one or more of: bold, italic, underline, strikethrough, dim, inverse"
        )

    if not isinstance(timestamp_style, bool):
        raise LogMessageError("timestamp_style must be a single logical value")


def output_cli_message(message: str, message_type: str) -> None:
    symbol = {
        "info": "ℹ",
        "success": "✔",
        "warning": "!",
        "running": "◌",
        "ask": "?",
    }.get(message_type, "")
    symbol_color = {
        "info": "cyan",
        "success": "green",
        "warning": "yellow",
        "running": "#FFA500",
        "ask": "magenta",
    }.get(message_type)
    if symbol and symbol_color:
        symbol = style_formatting(symbol, symbol_color, None, None)

    if symbol:
        print(f"{symbol} {message}", file=sys.stderr)
    else:
        print(message, file=sys.stderr)


def plain_text_output(
    text: str,
    cli_model: bool,
    text_color: Optional[str],
    back_color: Optional[str],
    text_style: Optional[Sequence[str]],
) -> None:
    if text_color is not None or back_color is not None or text_style:
        text = style_formatting(text, text_color, back_color, text_style)
    print(text, file=sys.stderr)


def format_line_with_style(
    line: str,
    prefix: str,
    text_color: Optional[str],
    back_color: Optional[str],
    text_style: Optional[Sequence[str]],
    timestamp_style: bool,
) -> str:
    if text_color is None and back_color is None and not text_style:
        return f"{prefix}{line}"

    if timestamp_style:
        return style_formatting(f"{prefix}{line}", text_color, back_color, text_style)

    styled_line = style_formatting(line, text_color, back_color, text_style)
    return f"{prefix}{styled_line}"


def output_message(
    msg: str,
    message_type: str,
    cli_model: bool,
    text_color: Optional[str],
    back_color: Optional[str],
    text_style: Optional[Sequence[str]],
    timestamp: bool,
    timestamp_format: str,
    level: int,
    symbol: str,
    multiline_indent: bool,
    timestamp_style: bool,
    plain_text: bool,
) -> None:
    if plain_text:
        if "\n" in msg:
            for line in msg.split("\n"):
                plain_text_output(line, cli_model, text_color, back_color, text_style)
        else:
            plain_text_output(msg, cli_model, text_color, back_color, text_style)
        return

    if cli_model and "\n" in msg:
        lines = msg.split("\n")
        for i, line in enumerate(lines, start=1):
            if i == 1 or multiline_indent:
                timestamp_part = timestamp_format if timestamp else ""
                indent_part = get_indent_part(symbol, level)
                prefix = f"{timestamp_part}{indent_part}"
            else:
                indent_part = get_indent_part(symbol, level)
                alignment_spaces = " " * len(timestamp_format) if timestamp else ""
                prefix = f"{alignment_spaces}{indent_part}"

            formatted_line = format_line_with_style(
                line=line,
                prefix=prefix,
                text_color=text_color,
                back_color=back_color,
                text_style=text_style,
                timestamp_style=timestamp_style,
            )
            output_cli_message(formatted_line, message_type)
        return

    if cli_model:
        timestamp_part = timestamp_format if timestamp else ""
        indent_part = get_indent_part(symbol, level)

        if symbol != "  ":
            final_msg = f"{timestamp_part}{symbol * level} {msg}"
        else:
            final_msg = f"{timestamp_part}{indent_part}{msg}"

        if text_color is not None or back_color is not None or text_style:
            if timestamp_style:
                final_msg = style_formatting(
                    final_msg, text_color, back_color, text_style
                )
            else:
                styled_msg = style_formatting(msg, text_color, back_color, text_style)
                symbol_part = (symbol * level + " ") if symbol != "  " else indent_part
                final_msg = f"{timestamp_part}{symbol_part}{styled_msg}"

        output_cli_message(final_msg, message_type)
        return

    formatted_msg = msg
    if text_color is not None or back_color is not None or text_style:
        formatted_msg = style_formatting(
            formatted_msg, text_color, back_color, text_style
        )

    prefix = {
        "info": "",
        "success": "SUCCESS: ",
        "warning": "WARNING: ",
        "running": "RUNNING: ",
        "ask": "? ",
    }.get(message_type, "")

    print(f"{prefix}{formatted_msg}", file=sys.stderr)


def ask_yes_no_cancel() -> Optional[bool]:
    while True:
        try:
            answer = input("[y/n/c]: ").strip().lower()
        except EOFError:
            return None

        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        if answer in {"c", "cancel", ""}:
            return None


def log_message(
    *args: object,
    verbose: Optional[bool] = True,
    message_type: str = "info",
    cli_model: bool = True,
    level: int = 1,
    symbol: str = "  ",
    text_color: Optional[str] = None,
    back_color: Optional[str] = None,
    text_style: Optional[Sequence[str]] = None,
    multiline_indent: bool = False,
    timestamp: bool = True,
    timestamp_format: Optional[str] = None,
    timestamp_style: bool = False,
    plain_text: bool = False,
):
    if message_type not in {"info", "success", "warning", "error", "running", "ask"}:
        raise LogMessageError(
            "message_type must be one of: info, success, warning, error, running, ask"
        )

    msg = build_message(*args)

    current_frame = inspect.currentframe()
    caller_frame = current_frame.f_back if current_frame else None
    msg = parse_inline_expressions(msg, caller_frame)

    if message_type == "error":
        raise LogMessageError(msg)

    verbose_value = get_verbose(verbose)
    if not verbose_value:
        return None

    validate_params(level, symbol, text_color, back_color, text_style, timestamp_style)

    if timestamp_format is None:
        timestamp_format = dt.datetime.now().strftime("[%Y-%m-%d %H:%M:%S] ")

    if message_type == "ask":
        output_message(
            msg=msg,
            message_type=message_type,
            cli_model=cli_model,
            text_color=text_color,
            back_color=back_color,
            text_style=text_style,
            timestamp=timestamp,
            timestamp_format=timestamp_format,
            level=level,
            symbol=symbol,
            multiline_indent=multiline_indent,
            timestamp_style=timestamp_style,
            plain_text=plain_text,
        )
        return ask_yes_no_cancel()

    output_message(
        msg=msg,
        message_type=message_type,
        cli_model=cli_model,
        text_color=text_color,
        back_color=back_color,
        text_style=text_style,
        timestamp=timestamp,
        timestamp_format=timestamp_format,
        level=level,
        symbol=symbol,
        multiline_indent=multiline_indent,
        timestamp_style=timestamp_style,
        plain_text=plain_text,
    )
    return None


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Print formatted message similar to thisutils::log_message"
    )
    parser.add_argument("message", nargs="*", help="Message parts")
    parser.add_argument("--verbose", default="true", help="true/false")
    parser.add_argument(
        "--message-type",
        default="info",
        choices=["ask", "error", "info", "running", "success", "warning"],
    )
    parser.add_argument("--cli-model", default="true", help="true/false")
    parser.add_argument("--level", type=int, default=1)
    parser.add_argument("--symbol", default="  ")
    parser.add_argument("--text-color", default=None)
    parser.add_argument("--back-color", default=None)
    parser.add_argument("--text-style", nargs="*", default=None)
    parser.add_argument("--multiline-indent", action="store_true")
    parser.add_argument("--timestamp", default="true", help="true/false")
    parser.add_argument("--timestamp-format", default=None)
    parser.add_argument("--timestamp-style", action="store_true")
    parser.add_argument("--plain-text", action="store_true")
    return parser.parse_args(argv)


def _to_bool(value: object, default: bool) -> bool:
    parsed = parse_bool(value)
    return default if parsed is None else parsed


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    try:
        result = log_message(
            *args.message,
            verbose=_to_bool(args.verbose, True),
            message_type=args.message_type,
            cli_model=_to_bool(args.cli_model, True),
            level=args.level,
            symbol=args.symbol,
            text_color=args.text_color,
            back_color=args.back_color,
            text_style=args.text_style,
            multiline_indent=args.multiline_indent,
            timestamp=_to_bool(args.timestamp, True),
            timestamp_format=args.timestamp_format,
            timestamp_style=args.timestamp_style,
            plain_text=args.plain_text,
        )
    except LogMessageError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if args.message_type == "ask":
        if result is True:
            print("TRUE")
        elif result is False:
            print("FALSE")
        else:
            print("NA")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
