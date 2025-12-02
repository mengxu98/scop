"""
Python implementation of log_message function
Provides formatted logging with timestamps, colors, and message types
"""

import sys
import os
import re
import inspect
from datetime import datetime
from typing import Optional, List


class LogMessage:
    """Log message formatter with color and style support"""

    # ANSI color codes
    COLORS = {
        "black": "\033[30m",
        "red": "\033[31m",
        "green": "\033[32m",
        "yellow": "\033[33m",
        "blue": "\033[34m",
        "magenta": "\033[35m",
        "cyan": "\033[36m",
        "white": "\033[37m",
        "grey": "\033[90m",
        "orange": "\033[38;2;255;165;0m",
        "br_red": "\033[91m",
        "br_green": "\033[92m",
        "br_yellow": "\033[93m",
        "br_blue": "\033[94m",
        "br_magenta": "\033[95m",
        "br_cyan": "\033[96m",
        "br_white": "\033[97m",
        "none": "\033[0m",
    }

    # Background colors
    BG_COLORS = {
        "black": "\033[40m",
        "red": "\033[41m",
        "green": "\033[42m",
        "yellow": "\033[43m",
        "blue": "\033[44m",
        "magenta": "\033[45m",
        "cyan": "\033[46m",
        "white": "\033[47m",
        "none": "\033[0m",
    }

    # Text styles
    STYLES = {
        "bold": "\033[1m",
        "dim": "\033[2m",
        "italic": "\033[3m",
        "underline": "\033[4m",
        "strikethrough": "\033[9m",
        "inverse": "\033[7m",
    }

    # Message type symbols and colors
    MESSAGE_TYPES = {
        "info": {"symbol": "ℹ", "color": "blue"},
        "success": {"symbol": "✓", "color": "green"},
        "warning": {"symbol": "!", "color": "yellow"},
        "error": {"symbol": "✗", "color": "red"},
        "running": {"symbol": "◌", "color": "orange"},
    }

    RESET = "\033[0m"

    # Inline format styles (cli-style)
    INLINE_FORMATS = {
        ".pkg": {"color": "blue", "style": []},
        ".code": {"color": "grey", "style": []},
        ".val": {"color": "blue", "style": []},
        ".arg": {"color": "none", "style": []},
        ".fun": {"color": "none", "style": []},
        ".file": {"color": "blue", "style": []},
        ".path": {"color": "blue", "style": []},
        ".field": {"color": "blue", "style": []},
        ".emph": {"color": "none", "style": ["italic"]},
        ".strong": {"color": "none", "style": ["bold"]},
    }

    def __init__(self):
        self._check_color_support()

    def _check_color_support(self):
        """Check if colors should be enabled (prefer enabled unless NO_COLOR is set)."""
        # disable colors (e.g., R's reticulate output), unless explicitly set NO_COLOR
        self.color_support = os.environ.get("NO_COLOR") is None

    def _hex_to_rgb(self, hex_color: str) -> tuple:
        """Convert hex color to RGB tuple"""
        hex_color = hex_color.lstrip("#")
        if len(hex_color) != 6:
            raise ValueError(f"Invalid hex color: #{hex_color}")
        return tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))

    def _rgb_to_ansi(self, rgb: tuple, bg: bool = False) -> str:
        """Convert RGB to ANSI color code"""
        r, g, b = rgb
        if bg:
            return f"\033[48;2;{r};{g};{b}m"
        else:
            return f"\033[38;2;{r};{g};{b}m"

    def _apply_color(self, text: str, color: Optional[str], bg: bool = False) -> str:
        """Apply color to text"""
        if not self.color_support or not color:
            return text

        # Handle hex colors
        if color.startswith("#"):
            try:
                rgb = self._hex_to_rgb(color)
                ansi_color = self._rgb_to_ansi(rgb, bg)
                return f"{ansi_color}{text}{self.RESET}"
            except ValueError:
                return text

        # Handle named colors
        color_map = self.BG_COLORS if bg else self.COLORS
        if color in color_map:
            return f"{color_map[color]}{text}{self.RESET}"

        return text

    def _apply_style(self, text: str, styles: Optional[List[str]]) -> str:
        """Apply text styles"""
        if not self.color_support or not styles:
            return text

        style_codes = []
        for style in styles:
            if style in self.STYLES:
                style_codes.append(self.STYLES[style])

        if style_codes:
            return f"{''.join(style_codes)}{text}{self.RESET}"

        return text

    def _get_indent(self, level: int, symbol: str) -> str:
        """Generate indentation string"""
        if symbol != "  ":
            return symbol * level + " "
        elif level > 1:
            return "  " * (level - 1)
        else:
            return ""

    def _format_message(
        self,
        message: str,
        message_type: str = "info",
        timestamp: bool = True,
        timestamp_format: str = "%Y-%m-%d %H:%M:%S",
        level: int = 1,
        symbol: str = "  ",
        text_color: Optional[str] = None,
        back_color: Optional[str] = None,
        text_style: Optional[List[str]] = None,
        multiline_indent: bool = False,
        timestamp_style: bool = True,
    ) -> str:
        """Format the complete message"""

        # Get message type info
        msg_info = self.MESSAGE_TYPES.get(message_type, self.MESSAGE_TYPES["info"])
        msg_symbol = msg_info["symbol"]
        msg_color = msg_info["color"]

        # Build timestamp
        timestamp_str = ""
        if timestamp:
            timestamp_str = f"[{datetime.now().strftime(timestamp_format)}] "

        # Build indentation
        indent = self._get_indent(level, symbol)

        # Prepare message-type symbol (with color)
        symbol_colored = (
            self._apply_color(msg_symbol, msg_color)
            if self.color_support
            else msg_symbol
        )
        symbol_part = f"{symbol_colored} "

        # Handle multiline messages
        if "\n" in message:
            lines = message.split("\n")
            formatted_lines = []

            for i, line in enumerate(lines):
                if i == 0 or multiline_indent:
                    # First line or multiline_indent=True: full formatting (symbol before timestamp)
                    prefix = symbol_part + timestamp_str + indent
                else:
                    # Subsequent lines: alignment spaces + indent
                    alignment_spaces = (
                        (" " * (len(symbol_part) + len(timestamp_str)))
                        if timestamp
                        else (" " * len(symbol_part))
                    )
                    prefix = alignment_spaces + indent

                # Apply formatting to line
                formatted_line = self._apply_formatting(
                    line, text_color, back_color, text_style, timestamp_style, prefix
                )
                formatted_lines.append(formatted_line)

            return "\n".join(formatted_lines)

        # Single line message (symbol before timestamp)
        prefix = symbol_part + timestamp_str + indent

        # Apply formatting
        formatted_msg = self._apply_formatting(
            message, text_color, back_color, text_style, timestamp_style, prefix
        )

        # Already added colored symbol in prefix
        return formatted_msg

    def _apply_formatting(
        self,
        text: str,
        text_color: Optional[str],
        back_color: Optional[str],
        text_style: Optional[List[str]],
        timestamp_style: bool,
        prefix: str,
    ) -> str:
        """Apply all formatting to text"""

        # Apply styles first
        if text_style:
            text = self._apply_style(text, text_style)

        # Apply text color
        if text_color:
            text = self._apply_color(text, text_color)

        # Apply background color
        if back_color:
            text = self._apply_color(text, back_color, bg=True)

        return prefix + text

    def _parse_inline_expressions(self, message: str, caller_frame=None) -> str:
        """
        Parse inline expressions with cli-style formatting

        Supports:
        - {.pkg package_name} / {pkg package_name}
        - {.pkg {variable}} / {pkg {variable}}
        - {.code some_code} / {code some_code}
        - {.val {expression}} / {val {expression}}
        - {expression}  # bare expression evaluation, no formatting
        etc.
        """
        max_iterations = 15
        iteration = 0

        allowed_tags = set(t.lstrip(".") for t in self.INLINE_FORMATS.keys())

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

            # 1) first parse the inline format: {.tag content} or {tag content}
            # content can be non-curly brace text or single layer {expr}
            format_match = re.search(
                r"\{\.?([A-Za-z_][A-Za-z0-9_]*)\s+(\{[^{}]*\}|[^{}]+)\}", message
            )

            if format_match:
                tag = format_match.group(1)
                content = format_match.group(2)

                # only process allowed tags
                if tag in allowed_tags:
                    if content.startswith("{") and content.endswith("}"):
                        evaluated = eval_expression(content[1:-1])
                    else:
                        evaluated = content

                    formatted = self._apply_inline_format(evaluated, f".{tag}")
                    message = (
                        message[: format_match.start()]
                        + formatted
                        + message[format_match.end() :]
                    )
                    # continue to next iteration
                    continue
                else:
                    # non-supported tag, skip to bare expression stage
                    pass

            # 2) then parse the bare expression: {expr} (not starting with .tag or tag)
            # to avoid conflict with format syntax, here exclude {.xxx ...} and {xxx ...}
            bare_match = re.search(r"\{([^{}]+)\}", message)
            if bare_match:
                inner = bare_match.group(1).strip()
                # if it looks like a format prefix (.tag or tag followed by space), skip this iteration
                if re.match(r"^\.?[A-Za-z_][A-Za-z0-9_]*\s+", inner):
                    # no bare expression to process, end loop
                    break
                evaluated = eval_expression(inner)
                message = (
                    message[: bare_match.start()]
                    + evaluated
                    + message[bare_match.end() :]
                )
                continue

            # no content to parse, end loop
            break

        return message

    def _apply_inline_format(self, text: str, format_type: str) -> str:
        """Apply formatting to inline content"""
        if format_type not in self.INLINE_FORMATS:
            return text

        fmt = self.INLINE_FORMATS[format_type]

        # apply color
        if fmt["color"] and fmt["color"] != "none":
            text = self._apply_color(text, fmt["color"])

        # apply style
        if fmt["style"]:
            text = self._apply_style(text, fmt["style"])

        return text

    def _format_traceback(self, depth: int = 1, skip_frames: int = 3) -> str:
        """
        Format traceback information for error messages

        Parameters:
        - depth: number of stack frames to show
        - skip_frames: number of internal frames to skip

        Returns:
        - formatted stack information string
        """
        stack = inspect.stack()

        # skip log_message internal call frames
        # skip_frames: _format_traceback, log_message, _format_message, etc.
        start_idx = skip_frames
        end_idx = min(start_idx + depth, len(stack))

        lines = []
        for i in range(start_idx, end_idx):
            frame_info = stack[i]
            # format file location information
            location = f'  File "{os.path.basename(frame_info.filename)}", line {frame_info.lineno}, in {frame_info.function}'
            lines.append(location)

            # show code line if available
            if frame_info.code_context:
                code_line = frame_info.code_context[0].strip()
                lines.append(f"    {code_line}")

        return "\n".join(lines) if lines else ""


# Global instance
_logger = LogMessage()


def log_message(
    *args,
    verbose: bool = True,
    message_type: str = "info",
    timestamp: bool = True,
    timestamp_format: str = "%Y-%m-%d %H:%M:%S",
    level: int = 1,
    symbol: str = "  ",
    text_color: Optional[str] = None,
    back_color: Optional[str] = None,
    text_style: Optional[List[str]] = None,
    multiline_indent: bool = False,
    timestamp_style: bool = True,
    show_traceback: bool = True,
    traceback_depth: int = 1,
) -> None:
    """
    Print formatted message with timestamp, colors, styling, and inline expressions

    Parameters:
    -----------
    *args : str
        Message parts to concatenate. Supports cli-style inline expressions:
        - {.pkg package_name} - Package names (cyan + bold)
        - {.code code_snippet} - Code snippets (grey)
        - {.val variable_name} - Variable values (blue)
        - {.arg parameter_name} - Function parameters (green)
        - {.fun function_name} - Function names (magenta)
        - {.file file_path} - File paths (underline)
        - {.path directory_path} - Directory paths (underline)
        - {.field field_name} - Object fields (cyan)
        - {.emph text} - Emphasized text (italic)
        - {.strong text} - Strong text (bold)

        Expressions can be nested: {.pkg {package_name}}
    verbose : bool, default True
        Whether to print the message
    message_type : str, default "info"
        Type of message: "info", "success", "warning", "error", "running"
    timestamp : bool, default True
        Whether to show timestamp
    timestamp_format : str, default "%Y-%m-%d %H:%M:%S"
        Timestamp format string
    level : int, default 1
        Indentation level
    symbol : str, default "  "
        Symbol used for indentation
    text_color : str, optional
        Text color (named color or hex code)
    back_color : str, optional
        Background color (named color or hex code)
    text_style : list, optional
        Text styles: ["bold", "italic", "underline", "dim", "strikethrough", "inverse"]
    multiline_indent : bool, default False
        Whether to apply formatting to each line in multiline messages
    timestamp_style : bool, default True
        Whether to apply styling to timestamp
    show_traceback : bool, default True
        Whether to show traceback for error messages
    traceback_depth : int, default 1
        Number of stack frames to show in traceback

    Returns:
    --------
    None
    """

    if not verbose:
        return

    # Validate message_type
    valid_types = ["info", "success", "warning", "error", "running"]
    if message_type not in valid_types:
        message_type = "info"

    # Get caller frame for inline expression evaluation and traceback
    # Need to handle both direct calls and calls through convenience functions
    current_frame = inspect.currentframe()
    caller_frame = current_frame.f_back

    # If called through convenience function (log_info, log_error, etc.), skip one more frame
    if caller_frame and caller_frame.f_code.co_name in [
        "log_info",
        "log_success",
        "log_warning",
        "log_error",
        "log_running",
    ]:
        caller_frame = caller_frame.f_back

    # Build message from args
    if not args:
        message = ""
    else:
        message = "".join(str(arg) for arg in args)

    # Parse inline expressions
    message = _logger._parse_inline_expressions(message, caller_frame)

    # Format and print message
    formatted_message = _logger._format_message(
        message=message,
        message_type=message_type,
        timestamp=timestamp,
        timestamp_format=timestamp_format,
        level=level,
        symbol=symbol,
        text_color=text_color,
        back_color=back_color,
        text_style=text_style,
        multiline_indent=multiline_indent,
        timestamp_style=timestamp_style,
    )

    # Add traceback for error messages
    if message_type == "error" and show_traceback:
        traceback_info = _logger._format_traceback(depth=traceback_depth)
        if traceback_info:
            formatted_message += "\n" + traceback_info

    print(formatted_message)


# Convenience functions for common message types
def log_info(*args, **kwargs):
    """Log info message"""
    log_message(*args, message_type="info", **kwargs)


def log_success(*args, **kwargs):
    """Log success message"""
    log_message(*args, message_type="success", **kwargs)


def log_warning(*args, **kwargs):
    """Log warning message"""
    log_message(*args, message_type="warning", **kwargs)


def log_error(*args, **kwargs):
    """Log error message"""
    log_message(*args, message_type="error", **kwargs)


def log_running(*args, **kwargs):
    """Log running message"""
    log_message(*args, message_type="running", **kwargs)


# SCVELO analysis function
def SCVELO(
    adata=None,
    h5ad=None,
    group_by=None,
    palette=None,
    linear_reduction=None,
    nonlinear_reduction=None,
    basis=None,
    mode=["deterministic", "stochastic", "dynamical"],
    fitting_by="stochastic",
    magic_impute=False,
    knn=5,
    t=2,
    min_shared_counts=30,
    n_pcs=30,
    n_neighbors=30,
    filter_genes=True,
    min_counts=3,
    min_counts_u=3,
    normalize_per_cell=True,
    log_transform=True,
    use_raw=False,
    diff_kinetics=False,
    stream_smooth=None,
    stream_density=2,
    arrow_length=5,
    arrow_size=5,
    arrow_density=0.5,
    denoise=False,
    denoise_topn=3,
    kinetics=False,
    kinetics_topn=100,
    calculate_velocity_genes=False,
    compute_velocity_confidence=True,
    compute_terminal_states=True,
    compute_pseudotime=True,
    compute_paga=True,
    top_n=6,
    n_jobs=1,
    show_plot=True,
    save=False,
    dpi=300,
    dirpath="./",
    fileprefix="",
    verbose=True,
):
    # Configure OpenMP settings to prevent conflicts
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    # Check if running on Apple silicon
    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"

    if is_apple_silicon:
        log_message(
            "Apple silicon detected: Applying specific configurations",
            message_type="info",
            verbose=verbose,
        )
        os.environ["PYTHONHASHSEED"] = "0"
        os.environ["PYTHONUNBUFFERED"] = "1"
        os.environ["SCANPY_SETTINGS"] = "scanpy_settings"
        os.environ["MPLBACKEND"] = "Agg"
        os.environ["DISPLAY"] = ""

        import numba

        numba.config.DISABLE_JIT = True  # Disable JIT compilation
        numba.set_num_threads(1)  # Force single thread
        log_message(
            "NUMBA configured for Apple silicon",
            message_type="success",
            verbose=verbose,
        )

    import matplotlib

    matplotlib.use("Agg")  # Use non-interactive backend
    import matplotlib.pyplot as plt
    import scvelo as scv
    import scanpy as sc
    import pandas as pd
    import numpy as np
    from scipy import sparse
    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(dirpath))

    log_message(
        "Starting {.pkg scVelo} analysis...", message_type="running", verbose=verbose
    )
    try:
        # Input validation
        if adata is None and h5ad is None:
            raise ValueError("Either {.arg adata} or {.arg h5ad} must be provided")

        if adata is None:
            adata = scv.read(h5ad)

        if group_by is None:
            raise ValueError("{.arg group_by} must be provided")

        if linear_reduction is None and nonlinear_reduction is None:
            raise ValueError(
                "At least one of {.arg linear_reduction} or {.arg nonlinear_reduction} must be provided"
            )

        # Setup basis
        if basis is None:
            if nonlinear_reduction is not None:
                # Check if the nonlinear reduction exists in obsm
                if nonlinear_reduction in adata.obsm:
                    basis = nonlinear_reduction
                elif f"X_{nonlinear_reduction}" in adata.obsm:
                    basis = f"X_{nonlinear_reduction}"
                else:
                    log_message(
                        "{.val {nonlinear_reduction}} not found in adata.obsm. Available keys: {.val {list(adata.obsm.keys())}}",
                        message_type="warning",
                        verbose=verbose,
                    )
                    basis = (
                        linear_reduction
                        if linear_reduction in adata.obsm
                        else f"X_{linear_reduction}"
                    )
            else:
                basis = (
                    linear_reduction
                    if linear_reduction in adata.obsm
                    else f"X_{linear_reduction}"
                )

        # Ensure the basis exists in obsm
        if basis not in adata.obsm:
            log_message(
                "basis '{.val {basis}}' not found in adata.obsm. Available keys: {.val {list(adata.obsm.keys())}}",
                message_type="warning",
                verbose=verbose,
            )
            # Try to find alternative basis
            if linear_reduction in adata.obsm:
                basis = linear_reduction
                log_message(
                    "Using {.val {linear_reduction}} as basis instead",
                    message_type="info",
                    verbose=verbose,
                )
            elif f"X_{linear_reduction}" in adata.obsm:
                basis = f"X_{linear_reduction}"
                log_message(
                    "Using {.val {linear_reduction}} as basis instead",
                    message_type="info",
                    verbose=verbose,
                )
            else:
                # Create a 2D basis from linear reduction
                if linear_reduction in adata.obsm:
                    adata.obsm["basis"] = adata.obsm[linear_reduction][:, 0:2]
                    basis = "basis"
                elif f"X_{linear_reduction}" in adata.obsm:
                    adata.obsm["basis"] = adata.obsm[f"X_{linear_reduction}"][:, 0:2]
                    basis = "basis"
                else:
                    raise ValueError(
                        "Cannot find suitable basis. Available obsm keys: {.val {list(adata.obsm.keys())}}"
                    )

        log_message("Using basis: {.val {basis}}", message_type="info", verbose=verbose)
        log_message(
            "Available embeddings in adata.obsm: {.val {list(adata.obsm.keys())}}",
            message_type="info",
            verbose=verbose,
        )

        # Ensure group_by is categorical
        adata.obs[group_by] = adata.obs[group_by].astype("category")

        # PREPROCESSING PHASE
        log_message("Starting preprocessing", message_type="running", verbose=verbose)

        # 1. Gene filtering (optional)
        if filter_genes:
            log_message("Filtering genes...", message_type="info", verbose=verbose)
            scv.pp.filter_genes(adata, min_counts=min_counts)
            scv.pp.filter_genes(adata, min_counts_u=min_counts_u)

        # 2. Normalization and transformation
        if normalize_per_cell:
            log_message("Normalizing per cell...", message_type="info", verbose=verbose)
            scv.pp.normalize_per_cell(adata)

        if log_transform:
            log_message("Log transforming...", message_type="info", verbose=verbose)
            sc.pp.log1p(adata)

        # 3. Magic imputation (if requested)
        if magic_impute:
            log_message(
                "Performing {.pkg magic-impute} imputation...",
                message_type="info",
                verbose=verbose,
            )
            try:
                import magic

                magic_operator = magic.MAGIC(knn=knn, t=t, random_state=42)
                adata.layers["spliced_raw"] = adata.layers["spliced"].copy()
                adata.layers["unspliced_raw"] = adata.layers["unspliced"].copy()
                adata.layers["spliced"] = magic_operator.fit_transform(
                    adata.layers["spliced"]
                )
                adata.layers["unspliced"] = magic_operator.transform(
                    adata.layers["unspliced"]
                )
            except ImportError:
                log_message(
                    "{.pkg magic-impute} not installed. Skipping imputation",
                    message_type="warning",
                    verbose=verbose,
                )

        # 4. Compute neighbors and moments with version compatibility
        log_message(
            "Computing {.pkg scvelo} neighbors and moments...",
            message_type="info",
            verbose=verbose,
        )

        # Find the best representation for neighbors computation
        # Priority: PCA > other linear reductions > fallback to raw data
        use_rep = None
        if linear_reduction in adata.obsm:
            rep_dims = adata.obsm[linear_reduction].shape[1]
            if rep_dims >= n_pcs:
                use_rep = linear_reduction
                log_message(
                    "Using {.val {linear_reduction}} with {.val {rep_dims}} dimensions",
                    message_type="info",
                    verbose=verbose,
                )
            else:
                log_message(
                    "{.val {linear_reduction}} has only {.val {rep_dims}} dimensions, need {.val {n_pcs}}",
                    message_type="warning",
                    verbose=verbose,
                )
        elif f"X_{linear_reduction}" in adata.obsm:
            rep_dims = adata.obsm[f"X_{linear_reduction}"].shape[1]
            if rep_dims >= n_pcs:
                use_rep = f"X_{linear_reduction}"
                log_message(
                    "Using {.val {linear_reduction}} with {.val {rep_dims}} dimensions",
                    message_type="info",
                    verbose=verbose,
                )
            else:
                log_message(
                    "{.val {linear_reduction}} has only {.val {rep_dims}} dimensions, need {.val {n_pcs}}",
                    message_type="warning",
                    verbose=verbose,
                )

        # Try to find a suitable PCA representation
        if use_rep is None:
            pca_candidates = [
                key
                for key in adata.obsm.keys()
                if "pca" in key.lower() and "umap" not in key.lower()
            ]
            for candidate in pca_candidates:
                rep_dims = adata.obsm[candidate].shape[1]
                if rep_dims >= n_pcs:
                    use_rep = candidate
                    log_message(
                        "Using PCA representation {.val {candidate}} with {.val {rep_dims}} dimensions",
                        message_type="info",
                        verbose=verbose,
                    )
                    break

        # If still no suitable representation, reduce n_pcs to available dimensions
        if use_rep is None:
            # Find the representation with the most dimensions
            max_dims = 0
            for key in adata.obsm.keys():
                if "pca" in key.lower() and "umap" not in key.lower():
                    dims = adata.obsm[key].shape[1]
                    if dims > max_dims:
                        max_dims = dims
                        use_rep = key

            if use_rep and max_dims > 0:
                n_pcs = min(n_pcs, max_dims)
                log_message(
                    "Reducing n_pcs to {.val {n_pcs}} to match available dimensions in {.val {use_rep}}",
                    message_type="info",
                    verbose=verbose,
                )
            else:
                # Fallback to raw data
                use_rep = None
                n_pcs = min(n_pcs, adata.X.shape[1])
                log_message(
                    "Using raw data with n_pcs={.val {n_pcs}}",
                    message_type="info",
                    verbose=verbose,
                )

        try:
            # Method 1: Try using scVelo's workflow
            scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors, use_rep=use_rep)
        except Exception as e:
            log_message(
                "{.pkg scvelo} moments failed ({.val {e}}), using manual computation...",
                message_type="warning",
                verbose=verbose,
            )
            # Method 2: Manual computation for compatibility
            sc.pp.neighbors(
                adata, n_pcs=n_pcs, n_neighbors=n_neighbors, use_rep=use_rep
            )

            # Manual moments calculation
            connectivities = adata.obsp["connectivities"]
            if sparse.issparse(adata.layers["spliced"]):
                Ms = connectivities @ adata.layers["spliced"]
                Mu = connectivities @ adata.layers["unspliced"]
            else:
                Ms = connectivities @ sparse.csr_matrix(adata.layers["spliced"])
                Mu = connectivities @ sparse.csr_matrix(adata.layers["unspliced"])

            adata.layers["Ms"] = Ms
            adata.layers["Mu"] = Mu

        log_message(
            "Starting {.pkg velocity} estimation",
            message_type="running",
            verbose=verbose,
        )

        for m in mode:
            log_message(f"Processing mode: {m}", message_type="info", verbose=verbose)

            if m == "dynamical":
                log_message(
                    "Performing {.pkg velocity} dynamical modeling...",
                    message_type="info",
                    verbose=verbose,
                )
                gene_subset = (
                    adata.var[fitting_by + "_genes"]
                    if (fitting_by + "_genes") in adata.var.columns
                    else None
                )

                if gene_subset is not None and gene_subset.sum() > 0:
                    scv.tl.recover_dynamics(
                        adata,
                        var_names=gene_subset,
                        use_raw=use_raw,
                        n_jobs=n_jobs,
                    )
                else:
                    log_message(
                        "No genes found for dynamical modeling. Using all genes",
                        message_type="warning",
                        verbose=verbose,
                    )
                    scv.tl.recover_dynamics(adata, use_raw=use_raw, n_jobs=n_jobs)

            log_message(
                "Computing {.pkg velocity} {.val {m}} velocity...",
                message_type="info",
                verbose=verbose,
            )
            scv.tl.velocity(adata, mode=m, diff_kinetics=diff_kinetics)

            log_message(
                "Computing {.pkg velocity} {.val {m}} graph...",
                message_type="info",
                verbose=verbose,
            )
            scv.tl.velocity_graph(adata, vkey=m, n_neighbors=n_neighbors, n_jobs=n_jobs)

            log_message(
                "Downstream analysis for {.pkg velocity} {.val {m}}",
                message_type="running",
                verbose=verbose,
            )

            log_message(
                "Computing {.pkg velocity} embedding...",
                message_type="info",
                verbose=verbose,
            )
            scv.tl.velocity_embedding(adata, basis=basis, vkey=m)

            # Velocity confidence (with error handling)
            if compute_velocity_confidence:
                log_message(
                    "Computing {.pkg velocity} confidence...",
                    message_type="info",
                    verbose=verbose,
                )
                try:
                    scv.tl.velocity_confidence(adata, vkey=m)
                except Exception as e:
                    log_message(
                        "velocity confidence failed ({.val {e}}), using default values",
                        message_type="warning",
                        verbose=verbose,
                    )
                    n_obs = adata.n_obs
                    adata.obs[m + "_length"] = np.ones(n_obs) * 0.5
                    adata.obs[m + "_confidence"] = np.ones(n_obs) * 0.5

            # Terminal states
            if compute_terminal_states:
                log_message(
                    "Computing terminal states...", message_type="info", verbose=verbose
                )
                try:
                    scv.tl.terminal_states(adata, vkey=m)
                    # Rename for consistency
                    for term in ["root_cells", "end_points"]:
                        if term in adata.obs.columns:
                            adata.obs[m + "_" + term] = adata.obs[term]
                            adata.obs.drop(term, axis=1, inplace=True)
                except Exception as e:
                    log_message(
                        "Terminal states computation failed: {.val {e}}",
                        message_type="warning",
                        verbose=verbose,
                    )

            # Pseudotime
            if compute_pseudotime:
                log_message(
                    "Computing {.pkg velocity} pseudotime...",
                    message_type="info",
                    verbose=verbose,
                )
                try:
                    root_key = (
                        m + "_root_cells"
                        if (m + "_root_cells") in adata.obs.columns
                        else None
                    )
                    end_key = (
                        m + "_end_points"
                        if (m + "_end_points") in adata.obs.columns
                        else None
                    )
                    scv.tl.velocity_pseudotime(
                        adata, vkey=m, root_key=root_key, end_key=end_key
                    )
                except Exception as e:
                    log_message(
                        "Pseudotime computation failed: {.val {e}}",
                        message_type="warning",
                        verbose=verbose,
                    )

            # PAGA
            if compute_paga:
                log_message("Computing PAGA...", message_type="info", verbose=verbose)
                try:
                    if "neighbors" not in adata.uns:
                        adata.uns["neighbors"] = {}
                    adata.uns["neighbors"]["distances"] = adata.obsp["distances"]
                    adata.uns["neighbors"]["connectivities"] = adata.obsp[
                        "connectivities"
                    ]

                    root_key = (
                        m + "_root_cells"
                        if (m + "_root_cells") in adata.obs.columns
                        else None
                    )
                    end_key = (
                        m + "_end_points"
                        if (m + "_end_points") in adata.obs.columns
                        else None
                    )
                    scv.tl.paga(
                        adata,
                        groups=group_by,
                        vkey=m,
                        root_key=root_key,
                        end_key=end_key,
                    )
                except Exception as e:
                    log_message(
                        "PAGA computation failed ({.val {e}})",
                        message_type="warning",
                        verbose=verbose,
                    )

            if calculate_velocity_genes:
                log_message(
                    "Ranking {.pkg velocity} genes...",
                    message_type="info",
                    verbose=verbose,
                )
                try:
                    if m != "dynamical":
                        scv.tl.rank_velocity_genes(adata, vkey=m, groupby=group_by)
                        if "spearmans_score" in adata.var.columns:
                            adata.var[m + "_score"] = adata.var["spearmans_score"]
                    else:
                        scv.tl.rank_dynamical_genes(adata, groupby=group_by)
                except Exception as e:
                    log_message(
                        "velocity genes ranking failed ({.val {e}})",
                        message_type="warning",
                        verbose=verbose,
                    )

            if show_plot:
                log_message(
                    "Generating plots for {.pkg {m}}...",
                    message_type="info",
                    verbose=verbose,
                )

                # Setup palette
                groups = (
                    adata.obs[group_by].cat.categories
                    if hasattr(adata.obs[group_by], "cat")
                    else adata.obs[group_by].unique()
                )
                if palette is None:
                    palette = dict(
                        zip(groups, plt.cm.tab10(np.linspace(0, 1, len(groups))))
                    )

                try:
                    scv.pl.velocity_embedding_stream(
                        adata,
                        vkey=m,
                        basis=basis,
                        title=f"{m} velocity",
                        color=group_by,
                        palette=palette,
                        smooth=stream_smooth,
                        density=stream_density,
                        legend_loc="right margin",
                        save=f"{fileprefix}_{m}_stream.pdf" if save else False,
                        show=show_plot,
                    )
                except Exception as e:
                    log_message(
                        "stream plot failed ({.val {e}})",
                        message_type="warning",
                        verbose=verbose,
                    )

                try:
                    scv.pl.velocity_embedding(
                        adata,
                        vkey=m,
                        basis=basis,
                        title=f"{m} velocity",
                        color=group_by,
                        palette=palette,
                        arrow_length=arrow_length,
                        arrow_size=arrow_size,
                        density=arrow_density,
                        linewidth=0.3,
                        save=f"{fileprefix}_{m}_arrow.pdf" if save else False,
                        show=show_plot,
                    )
                except Exception as e:
                    log_message(
                        "arrow plot failed ({.val {e}})",
                        message_type="warning",
                        verbose=verbose,
                    )

                if compute_velocity_confidence:
                    for metric in ["length", "confidence"]:
                        color_key = f"{m}_{metric}"
                        if color_key in adata.obs.columns:
                            try:
                                scv.pl.scatter(
                                    adata,
                                    basis=basis,
                                    color=color_key,
                                    title=f"{m} {metric}",
                                    cmap="viridis",
                                    save=f"{fileprefix}_{m}_{metric}.pdf"
                                    if save
                                    else False,
                                    show=show_plot,
                                )
                            except Exception as e:
                                log_message(
                                    "{.pkg {metric}} plot failed ({.val {e}})",
                                    message_type="warning",
                                    verbose=verbose,
                                )

        log_message(
            "{.pkg scVelo} analysis completed", message_type="success", verbose=verbose
        )

    except Exception as e:
        log_message(
            "Error in {.pkg scVelo} analysis: {.val {e}}",
            message_type="error",
            verbose=verbose,
        )
        raise
    finally:
        os.chdir(prevdir)

    try:
        if hasattr(adata, "_raw") and adata._raw is not None:
            if hasattr(adata._raw, "_var"):
                adata._raw._var = adata._raw._var.rename(columns={"_index": "features"})
    except Exception:
        pass

    return adata


def CellRank(
    adata=None,
    h5ad=None,
    group_by=None,
    palette=None,
    linear_reduction=None,
    nonlinear_reduction=None,
    basis=None,
    mode=["deterministic", "stochastic", "dynamical"],
    fitting_by="stochastic",
    magic_impute=False,
    knn=5,
    t=2,
    min_shared_counts=30,
    n_pcs=30,
    n_neighbors=30,
    stream_smooth=None,
    stream_density=2,
    arrow_size=5,
    arrow_length=5,
    arrow_density=0.5,
    calculate_velocity_genes=False,
    denoise=False,
    kinetics=False,
    n_jobs=1,
    show_plot=True,
    dpi=300,
    save=False,
    dirpath="./",
    fileprefix="",
    verbose=True,
):
    # Configure OpenMP settings to prevent conflicts
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    # Check if running on Apple silicon
    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"

    if is_apple_silicon:
        log_message(
            "Apple silicon detected: Applying specific configurations",
            message_type="info",
            verbose=verbose,
        )
        os.environ["PYTHONHASHSEED"] = "0"
        os.environ["PYTHONUNBUFFERED"] = "1"
        os.environ["SCANPY_SETTINGS"] = "scanpy_settings"
        os.environ["MPLBACKEND"] = "Agg"
        os.environ["DISPLAY"] = ""

        import numba

        numba.config.DISABLE_JIT = True # Disable JIT compilation
        numba.set_num_threads(1) # Force single thread
        log_message(
            "NUMBA configured for Apple silicon",
            message_type="success",
            verbose=verbose,
        )

    import matplotlib.pyplot as plt
    import scvelo as scv
    import cellrank as cr
    import pandas as pd
    import numpy as np
    import scanpy as sc

    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(dirpath))

    if platform.system() == "Windows":
        import sys, multiprocessing, re

        if re.match(pattern=".*pythonw.exe$", string=sys.executable):
            pythonw = sys.executable
        else:
            pythonw = sys.executable.replace("python.exe", "pythonw.exe")
        sys.executable = pythonw
        sys._base_executable = pythonw
        multiprocessing.set_executable(pythonw)

    try:
        if adata is None and h5ad is None:
            log_message("adata or h5ad must be provided.", message_type="error")
            exit()

        if adata is None:
            adata = scv.read(h5ad)
        # del adata.uns

        if group_by is None:
            log_message("group_by must be provided.", message_type="error")
            exit()

        if linear_reduction is None and nonlinear_reduction is None:
            log_message(
                "linear_reduction or nonlinear_reduction must be provided at least one.",
                message_type="error",
            )
            exit()

        if linear_reduction is None:
            sc.pp.pca(adata, n_comps=n_pcs)
            linear_reduction = "X_pca"

        if basis is None:
            if nonlinear_reduction is not None:
                basis = nonlinear_reduction
            else:
                basis = "basis"
                adata.obsm["basis"] = adata.obsm[linear_reduction][:, 0:2]

        mode.append(fitting_by)
        if kinetics is True or denoise is True:
            mode.append("dynamical")

        mode = list(set(mode))
        if "dynamical" in mode:
            mode.sort(key="dynamical".__eq__)

        if not fitting_by in ["deterministic", "stochastic"]:
            log_message(
                "'fitting_by' must be one of 'deterministic' and 'stochastic'.",
                message_type="error",
            )
            exit()

        if not all([m in ["deterministic", "stochastic", "dynamical"] for m in mode]):
            log_message(
                "Invalid mode name! Must be the 'deterministic', 'stochastic' or 'dynamical'.",
                message_type="error",
            )
            exit()

        adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")

        if mode[-1] + "_graph" not in adata.obs.keys():
            adata = SCVELO(
                adata=adata,
                group_by=group_by,
                n_jobs=n_jobs,
                linear_reduction=linear_reduction,
                nonlinear_reduction=nonlinear_reduction,
                basis=basis,
                mode=mode,
                fitting_by=fitting_by,
                magic_impute=magic_impute,
                knn=knn,
                t=t,
                min_shared_counts=min_shared_counts,
                n_pcs=n_pcs,
                n_neighbors=n_neighbors,
                stream_smooth=stream_smooth,
                stream_density=stream_density,
                arrow_size=arrow_size,
                arrow_length=arrow_length,
                arrow_density=arrow_density,
                denoise=denoise,
                kinetics=kinetics,
                calculate_velocity_genes=calculate_velocity_genes,
                show_plot=show_plot,
                dpi=dpi,
                save=save,
                dirpath=dirpath,
                fileprefix=fileprefix,
            )
        adata.layers["velocity"] = adata.layers[mode[-1]]
        cr.tl.terminal_states(adata, cluster_key=group_by)
        cr.pl.terminal_states(adata)
        cr.tl.initial_states(adata, cluster_key=group_by)
        cr.pl.initial_states(adata)
        cr.tl.lineages(adata)
        cr.pl.lineages(adata, same_plot=False)
        cr.pl.lineages(adata, same_plot=True)

        # Check if dynamical mode was used, if not, run recover_dynamics for latent time
        if "dynamical" not in mode:
            log_message(
                "Running recover_dynamics for latent time computation...",
                message_type="info",
                verbose=verbose,
            )
            try:
                scv.tl.recover_dynamics(adata, use_raw=False, n_jobs=n_jobs)
            except Exception as e:
                log_message(
                    f"recover_dynamics failed ({e}), skipping latent time computation...",
                    message_type="warning",
                    verbose=verbose,
                )
                pass
            else:
                scv.tl.recover_latent_time(
                    adata,
                    root_key="initial_states_probs",
                    end_key="terminal_states_probs",
                )
        else:
            scv.tl.recover_latent_time(
                adata, root_key="initial_states_probs", end_key="terminal_states_probs"
            )
        scv.tl.paga(
            adata,
            groups=group_by,
            root_key="initial_states_probs",
            end_key="terminal_states_probs",
            use_time_prior="velocity_pseudotime",
        )
        cr.pl.cluster_fates(
            adata,
            mode="paga_pie",
            cluster_key=group_by,
            basis=basis,
            legend_kwargs={"loc": "top right out"},
            legend_loc="top left out",
            node_size_scale=5,
            edge_width_scale=1,
            max_edge_width=4,
            title="directed PAGA",
        )
        if show_plot is True:
            plt.show()

        cr.tl.lineage_drivers(adata, cluster_key=group_by)
        cr.pl.lineage_drivers(adata, lineage=adata.obs[group_by].unique()[1], n_genes=4)
        if show_plot is True:
            plt.show()

    finally:
        os.chdir(prevdir)

    try:
        adata.__dict__["_raw"].__dict__["_var"] = (
            adata.__dict__["_raw"]
            .__dict__["_var"]
            .rename(columns={"_index": "features"})
        )
    except:
        pass

    return adata


def PAGA(
    adata=None,
    h5ad=None,
    group_by=None,
    palette=None,
    linear_reduction=None,
    nonlinear_reduction=None,
    basis=None,
    n_pcs=30,
    n_neighbors=30,
    use_rna_velocity=False,
    vkey="stochastic",
    embedded_with_PAGA=False,
    paga_layout="fr",
    threshold=0.1,
    point_size=20,
    infer_pseudotime=False,
    root_cell=None,
    root_group=None,
    n_dcs=10,
    n_branchings=0,
    min_group_size=0.01,
    show_plot=True,
    dpi=300,
    save=False,
    dirpath="./",
    fileprefix="",
    verbose=True,
):
    # Configure OpenMP settings to prevent conflicts
    import os
    import platform
    import sys

    # Enhanced thread configuration for Apple silicon MacBooks
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    # Additional Apple silicon specific configurations
    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"
    if is_apple_silicon:
        os.environ["PYTHONHASHSEED"] = "0"
        os.environ["PYTHONUNBUFFERED"] = "1"
        os.environ["SCANPY_SETTINGS"] = "scanpy_settings"
        os.environ["SCANPY_SETTINGS_VERBOSITY"] = "1"
        # Disable problematic features for Apple silicon
        os.environ["MPLBACKEND"] = "Agg"
        os.environ["DISPLAY"] = ""
        # Force single-threaded execution - set early
        os.environ["NUMBA_NUM_THREADS"] = "1"
        os.environ["NUMBA_DISABLE_JIT"] = "1"
        # Additional thread control
        os.environ["NUMBA_THREADING_LAYER"] = "tbb"
        os.environ["NUMBA_DEFAULT_NUM_THREADS"] = "1"

    # Import with error handling
    try:
        import matplotlib

        matplotlib.use("Agg")  # Use non-interactive backend
        import matplotlib.pyplot as plt

        # Configure matplotlib for Apple silicon
        if is_apple_silicon:
            plt.rcParams["figure.max_open_warning"] = 0
            plt.rcParams["figure.figsize"] = (4, 4)  # Small figures
            plt.rcParams["axes.linewidth"] = 0.5
            plt.rcParams["lines.linewidth"] = 0.5
            plt.rcParams["font.size"] = 8
    except Exception as e:
        log_message(
            f"matplotlib setup failed: {e}", message_type="warning", verbose=verbose
        )
        plt = None

    # Set NUMBA threads before importing scanpy
    if is_apple_silicon:
        try:
            import numba

            # Completely disable JIT compilation for Apple silicon
            numba.config.DISABLE_JIT = True
            numba.set_num_threads(1)
        except:
            pass

    try:
        import scanpy as sc
    except Exception as e:
        log_message(f"scanpy import failed: {e}", message_type="error")
        raise

    try:
        import numpy as np
    except Exception as e:
        log_message(f"numpy import failed: {e}", message_type="error")
        raise

    import statistics
    from math import hypot
    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    # Additional Apple silicon configurations are already set above

    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(dirpath))

    import platform

    if platform.system() == "Windows":
        import sys, multiprocessing, re

        if re.match(pattern=".*pythonw.exe$", string=sys.executable):
            pythonw = sys.executable
        else:
            pythonw = sys.executable.replace("python.exe", "pythonw.exe")
        sys.executable = pythonw
        sys._base_executable = pythonw
        multiprocessing.set_executable(pythonw)

    try:
        if adata is None and h5ad is None:
            log_message("adata or h5ad must be provided.", message_type="error")
            return None

        if adata is None:
            adata = sc.read(h5ad)

        if group_by is None:
            log_message("group_by must be provided.", message_type="error")
            return None

        if linear_reduction is None and nonlinear_reduction is None:
            log_message(
                "linear_reduction or nonlinear_reduction must be provided at least one.",
                message_type="error",
            )
            return None

        if linear_reduction is None:
            sc.pp.pca(adata, n_comps=n_pcs)
            linear_reduction = "X_pca"

        if basis is None:
            if nonlinear_reduction is not None:
                basis = nonlinear_reduction
            else:
                basis = "basis"
                adata.obsm["basis"] = adata.obsm[linear_reduction][:, 0:2]

        if point_size is None:
            point_size = min(100000 / adata.shape[0], 20)

        if infer_pseudotime is True and root_cell is None and root_group is None:
            log_message(
                "root_cell or root_group should be provided.", message_type="error"
            )
            return None

        if use_rna_velocity is True:
            adata.uns["velocity_graph"] = adata.uns[vkey + "_graph"]

        adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")

        # Choose representation key and make n_pcs safe for its dimensionality
        rep_key = linear_reduction
        try:
            obsm_keys = set(adata.obsm_keys())
            cand_keys = [
                linear_reduction,
                (
                    f"X_{linear_reduction}"
                    if linear_reduction is not None
                    and not str(linear_reduction).startswith("X_")
                    else None
                ),
            ]
            rep_key = next(
                (k for k in cand_keys if k is not None and k in obsm_keys),
                linear_reduction,
            )

            neighbors_n_pcs = n_pcs
            if rep_key in obsm_keys:
                rep_dims = adata.obsm[rep_key].shape[1]
                if neighbors_n_pcs is not None and isinstance(
                    neighbors_n_pcs, (int, float)
                ):
                    # If the representation has fewer dims than requested PCs, let Scanpy use full dims
                    if rep_dims < int(neighbors_n_pcs):
                        neighbors_n_pcs = None
            else:
                neighbors_n_pcs = n_pcs
        except Exception:
            neighbors_n_pcs = n_pcs

        if "X_diffmap" in adata.obsm_keys():
            X_diffmap = adata.obsm["X_diffmap"]
            del adata.obsm["X_diffmap"]
            sc.pp.neighbors(
                adata, n_pcs=neighbors_n_pcs, use_rep=rep_key, n_neighbors=n_neighbors
            )
            adata.obsm["X_diffmap"] = X_diffmap
        else:
            sc.pp.neighbors(
                adata, n_pcs=neighbors_n_pcs, use_rep=rep_key, n_neighbors=n_neighbors
            )

        # Core PAGA computation - this is the essential part
        sc.tl.paga(adata, groups=group_by, use_rna_velocity=use_rna_velocity)

        # PAGA compare plot with error handling
        try:
            if use_rna_velocity is True:
                sc.pl.paga_compare(
                    adata,
                    basis=basis,
                    palette=palette,
                    threshold=threshold,
                    size=point_size,
                    min_edge_width=1,
                    node_size_scale=1,
                    dashed_edges="connectivities",
                    solid_edges="transitions_confidence",
                    transitions="transitions_confidence",
                    title=basis,
                    frameon=False,
                    edges=True,
                    save=False,
                    show=show_plot,
                )
            else:
                sc.pl.paga_compare(
                    adata,
                    basis=basis,
                    palette=palette,
                    threshold=threshold,
                    size=point_size,
                    title=basis,
                    frameon=False,
                    edges=True,
                    save=False,
                    show=show_plot,
                )

            if save:
                plt.savefig(
                    ".".join(filter(None, [fileprefix, "paga_compare.pdf"])),
                    dpi=dpi,
                    bbox_inches="tight",
                    facecolor="white",
                )
        except Exception as e:
            log_message(
                f"PAGA compare plot failed ({e}), continuing...",
                message_type="warning",
            )
            if plt is not None:
                plt.clf()
                plt.close("all")

        # PAGA layout plot
        try:
            sc.pl.paga(
                adata,
                threshold=threshold,
                layout=paga_layout,
                title="PAGA layout: " + paga_layout,
                frameon=False,
                save=False,
                show=show_plot,
            )

            if save:
                plt.savefig(
                    ".".join(filter(None, [fileprefix, "paga_layout.pdf"])),
                    dpi=dpi,
                    bbox_inches="tight",
                    facecolor="white",
                )
        except Exception as e:
            log_message(
                f"PAGA layout plot failed ({e}), continuing...",
                message_type="warning",
            )
            if plt is not None:
                plt.clf()
                plt.close("all")

        # draw_graph computation and plotting
        try:
            sc.tl.draw_graph(adata, init_pos="paga", layout=paga_layout)
            log_message(
                "draw_graph computation completed successfully", message_type="success"
            )

            sc.pl.draw_graph(
                adata,
                color=group_by,
                palette=palette,
                title="PAGA layout: " + paga_layout,
                layout=paga_layout,
                frameon=False,
                legend_loc="on data",
                show=show_plot,
            )

            if save:
                plt.savefig(
                    ".".join(filter(None, [fileprefix, "paga_graph.pdf"])),
                    dpi=dpi,
                    bbox_inches="tight",
                    facecolor="white",
                )
        except Exception as e:
            log_message(
                f"PAGA draw_graph failed ({e}), continuing...", message_type="warning"
            )
            if plt is not None:
                plt.clf()
                plt.close("all")

        if embedded_with_PAGA is True:
            # UMAP computation and plotting
            try:
                umap2d = sc.tl.umap(adata, init_pos="paga", n_components=2, copy=True)
                adata.obsm["PAGAUMAP2D"] = umap2d.obsm["X_umap"]

                sc.pl.paga_compare(
                    adata,
                    basis="PAGAUMAP2D",
                    palette=palette,
                    threshold=threshold,
                    size=point_size,
                    title="PAGA-initialized UMAP",
                    edges=True,
                    save=False,
                    show=show_plot,
                )

                if save:
                    plt.savefig(
                        ".".join(filter(None, [fileprefix, "paga_umap.pdf"])),
                        dpi=dpi,
                        bbox_inches="tight",
                        facecolor="white",
                    )
            except Exception as e:
                log_message(
                    f"PAGA UMAP failed ({e}), continuing...", message_type="warning"
                )
                if plt is not None:
                    plt.clf()
                    plt.close("all")

        if infer_pseudotime is True:
            if root_group is not None and root_cell is None:
                cell = adata.obs[group_by].index.values[
                    adata.obs[group_by] == root_group
                ]
                root_group_cell = adata.obsm[basis][adata.obs[group_by] == root_group,][
                    :, [0, 1]
                ]
                x = statistics.median(root_group_cell[:, 0])
                y = statistics.median(root_group_cell[:, 1])
                diff = np.array((x - root_group_cell[:, 0], y - root_group_cell[:, 1]))
                dist = []
                for i in range(diff.shape[1]):
                    dist.append(hypot(diff[0, i], diff[1, i]))

                root_cell = cell[dist.index(min(dist))]

            sc.tl.diffmap(adata, n_comps=n_dcs)
            adata.uns["iroot"] = np.flatnonzero(adata.obs_names == root_cell)[0]
            sc.tl.dpt(
                adata,
                n_dcs=n_dcs,
                n_branchings=n_branchings,
                min_group_size=min_group_size,
            )

            # Pseudotime plotting
            try:
                sc.pl.embedding(
                    adata,
                    basis=basis,
                    color="dpt_pseudotime",
                    save=False,
                    show=show_plot,
                )

                if save:
                    plt.savefig(
                        ".".join(filter(None, [fileprefix, "dpt_pseudotime.pdf"])),
                        dpi=dpi,
                        bbox_inches="tight",
                        facecolor="white",
                    )
            except Exception as e:
                log_message(
                    f"DPT pseudotime plot failed ({e}), continuing...",
                    message_type="warning",
                )
                if plt is not None:
                    plt.clf()
                    plt.close("all")

    finally:
        os.chdir(prevdir)

        # Clean up matplotlib for Apple silicon MacBooks
        if is_apple_silicon and plt is not None:
            try:
                plt.clf()
                plt.close("all")
                # Force garbage collection
                import gc

                gc.collect()
            except:
                pass

    try:
        adata.__dict__["_raw"].__dict__["_var"] = (
            adata.__dict__["_raw"]
            .__dict__["_var"]
            .rename(columns={"_index": "features"})
        )
    except:
        pass

    return adata


def Palantir(
    adata=None,
    h5ad=None,
    group_by=None,
    palette=None,
    linear_reduction=None,
    nonlinear_reduction=None,
    basis=None,
    n_pcs=30,
    n_neighbors=30,
    dm_n_components=10,
    dm_alpha=0,
    dm_n_eigs=None,
    early_group=None,
    terminal_groups=None,
    early_cell=None,
    terminal_cells=None,
    num_waypoints=1200,
    scale_components=True,
    use_early_cell_as_start=False,
    adjust_early_cell=False,
    adjust_terminal_cells=False,
    max_iterations=25,
    n_jobs=1,
    point_size=20,
    show_plot=True,
    dpi=300,
    save=False,
    dirpath="./",
    fileprefix="",
    verbose=True,
):
    # Configure OpenMP settings to prevent conflicts
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    # Check if running on M-series machine
    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"

    if is_apple_silicon:
        log_message(
            "Apple silicon detected: Applying specific configurations",
            message_type="info",
            verbose=verbose,
        )

        os.environ["PYTHONHASHSEED"] = "0"
        os.environ["PYTHONUNBUFFERED"] = "1"
        os.environ["SCANPY_SETTINGS"] = "scanpy_settings"
        os.environ["MPLBACKEND"] = "Agg"
        os.environ["DISPLAY"] = ""

        import numba

        numba.config.DISABLE_JIT = True  # Disable JIT compilation
        numba.set_num_threads(1)  # Force single thread
        log_message(
            "NUMBA configured for Apple silicon",
            message_type="success",
            verbose=verbose,
        )

    import matplotlib.pyplot as plt
    import matplotlib
    import statistics
    from math import hypot
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import palantir

    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(dirpath))

    if platform.system() == "Windows":
        import sys, multiprocessing, re

        if re.match(pattern=".*pythonw.exe$", string=sys.executable):
            pythonw = sys.executable
        else:
            pythonw = sys.executable.replace("python.exe", "pythonw.exe")
        sys.executable = pythonw
        sys._base_executable = pythonw
        multiprocessing.set_executable(pythonw)

    try:
        if adata is None and h5ad is None:
            log_message("adata or h5ad must be provided.", message_type="error")
            exit()

        if adata is None:
            adata = sc.read(h5ad)
        # del adata.uns

        if group_by is None and (
            early_group is not None or terminal_groups is not None
        ):
            log_message("`group_by` must be provided", message_type="error")
            exit()

        if linear_reduction is None and nonlinear_reduction is None:
            log_message(
                "`linear_reduction` or `nonlinear_reduction` must be provided at least one",
                message_type="error",
            )
            exit()

        if linear_reduction is None:
            sc.pp.pca(adata, n_comps=n_pcs)
            linear_reduction = "X_pca"

        if basis is None:
            if nonlinear_reduction is not None:
                basis = nonlinear_reduction
            else:
                basis = linear_reduction

        if point_size is None:
            point_size = min(100000 / adata.shape[0], 20)

        if early_group is not None and early_cell is None:
            cell = adata.obs[group_by].index.values[adata.obs[group_by] == early_group]
            early_group_cell = adata.obsm[basis][adata.obs[group_by] == early_group,][
                :, [0, 1]
            ]
            x = statistics.median(early_group_cell[:, 0])
            y = statistics.median(early_group_cell[:, 1])
            diff = np.array((x - early_group_cell[:, 0], y - early_group_cell[:, 1]))
            dist = []
            for i in range(diff.shape[1]):
                dist.append(hypot(diff[0, i], diff[1, i]))

            early_cell = cell[dist.index(min(dist))]

        if early_cell is None:
            log_message("`early_cell` must be provided", message_type="error")
            exit()
        else:
            log_message(
                f"`early_cell`: {early_cell}", message_type="info", verbose=verbose
            )

        terminal_cells_dict = dict()
        if terminal_groups is not None and terminal_cells is None:
            for n in range(len(terminal_groups)):
                terminal_group = terminal_groups[n]
                cell = adata.obs[group_by].index.values[
                    adata.obs[group_by] == terminal_group
                ]
                terminal_group_cell = adata.obsm[basis][
                    adata.obs[group_by] == terminal_group,
                ][:, [0, 1]]
                x = statistics.median(terminal_group_cell[:, 0])
                y = statistics.median(terminal_group_cell[:, 1])
                diff = np.array(
                    (x - terminal_group_cell[:, 0], y - terminal_group_cell[:, 1])
                )
                dist = []
                for i in range(diff.shape[1]):
                    dist.append(hypot(diff[0, i], diff[1, i]))
                terminal_cells_dict[cell[dist.index(min(dist))]] = (
                    terminal_group.replace(" ", ".") + "_diff_potential"
                )

            terminal_cells = list(terminal_cells_dict.keys())

        if terminal_cells is None:
            log_message("`terminal_cells`: None", message_type="info", verbose=verbose)
        else:
            log_message(
                f"`terminal_cells`: {terminal_cells}",
                message_type="info",
                verbose=verbose,
            )

        # if HVF is None:
        #   sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        # else:
        #   df = pd.DataFrame([False] * adata.X.shape[1],columns=["highly_variable"])
        #   df = df.set_index(adata.var_names)
        #   df.highly_variable.iloc[:n_top_genes] = True
        #   df.loc[df['channel'].isin(['sale','fullprice'])]
        #   df.highly_variable.iloc[df.index.isin(HVF)] = True
        #   if "highly_variable" in adata.var.columns:
        #     adata.var.drop('highly_variable', axis=1, inplace=True)
        #   adata.var=adata.var.join(df)

        # adata.uns['pca']['variance_ratio']
        # pca_projections=n_comps = np.where(np.cumsum(ad.uns['pca']['variance_ratio']) > 0.85)[0][0]
        # pca_projections, _ = palantir.utils.run_pca(adata, use_hvg=True)

        pca_projections = pd.DataFrame(
            adata.obsm[linear_reduction][:, :n_pcs], index=adata.obs_names
        )
        log_message("Running diffusion maps", message_type="info", verbose=verbose)
        dm_res = palantir.utils.run_diffusion_maps(
            pca_projections,
            n_components=dm_n_components,
            knn=n_neighbors,
            alpha=dm_alpha,
        )
        ms_data = palantir.utils.determine_multiscale_space(dm_res, n_eigs=dm_n_eigs)
        log_message("Running palantir", message_type="info", verbose=verbose)
        pr_res = palantir.core.run_palantir(
            data=ms_data,
            early_cell=early_cell,
            terminal_states=terminal_cells,
            knn=n_neighbors,
            num_waypoints=num_waypoints,
            scale_components=scale_components,
            use_early_cell_as_start=use_early_cell_as_start,
            max_iterations=max_iterations,
            n_jobs=n_jobs,
        )

        if adjust_early_cell is True or adjust_terminal_cells is True:
            if adjust_early_cell is True:
                early_cell_group = adata.obs[group_by][early_cell]
                cells = adata.obs[group_by].index.values[
                    adata.obs[group_by] == early_cell_group
                ]
                early_cell = pr_res.pseudotime[cells].index.values[
                    pr_res.pseudotime[cells] == min(pr_res.pseudotime[cells])
                ][0]
            if adjust_terminal_cells is True:
                terminal_cells_dict = dict()
                for n in range(len(terminal_cells)):
                    terminal_cell = terminal_cells[n]
                    terminal_cell_group = adata.obs[group_by][terminal_cell]
                    cells = adata.obs[group_by].index.values[
                        adata.obs[group_by] == terminal_cell_group
                    ]
                    terminal_cells_dict[
                        pr_res.pseudotime[cells].index.values[
                            pr_res.pseudotime[cells] == max(pr_res.pseudotime[cells])
                        ][0]
                    ] = terminal_cell_group.replace(" ", ".") + "_diff_potential"
                terminal_cells = list(terminal_cells_dict.keys())

            pr_res = palantir.core.run_palantir(
                data=ms_data,
                early_cell=early_cell,
                terminal_states=terminal_cells,
                knn=n_neighbors,
                num_waypoints=num_waypoints,
                scale_components=scale_components,
                use_early_cell_as_start=use_early_cell_as_start,
                max_iterations=max_iterations,
                n_jobs=n_jobs,
            )

        adata.obsm["palantir_dm"] = dm_res["T"].toarray()
        adata.uns["dm_kernel"] = dm_res["kernel"]
        if len(terminal_cells_dict) > 0:
            pr_res.branch_probs = pr_res.branch_probs.rename(
                columns=terminal_cells_dict
            )
        for term in np.append(
            pr_res.branch_probs.columns.values,
            np.array(["palantir_pseudotime", "palantir_diff_potential"]),
        ):
            if term in adata.obs.columns:
                adata.obs.drop(term, axis=1, inplace=True)
        adata.obs = adata.obs.join(pr_res.pseudotime.to_frame("palantir_pseudotime"))
        adata.obs = adata.obs.join(pr_res.entropy.to_frame("palantir_diff_potential"))
        adata.obs = adata.obs.join(pr_res.branch_probs)

        sc.pl.embedding(
            adata,
            basis=basis,
            color="palantir_pseudotime",
            size=point_size,
            show=show_plot,
        )
        if save:
            plt.savefig(
                ".".join(filter(None, [fileprefix, "palantir_pseudotime.pdf"])), dpi=dpi
            )

        sc.pl.embedding(
            adata,
            basis=basis,
            color="palantir_diff_potential",
            size=point_size,
            show=show_plot,
        )
        if save:
            plt.savefig(
                ".".join(filter(None, [fileprefix, "palantir_diff_potential.pdf"])),
                dpi=dpi,
            )

        sc.pl.embedding(
            adata,
            basis=basis,
            color=pr_res.branch_probs.columns.values,
            size=point_size,
            show=show_plot,
        )
        if save:
            plt.savefig(
                ".".join(filter(None, [fileprefix, "palantir_probs.pdf"])), dpi=dpi
            )

        if group_by is not None:
            sc.pl.embedding(
                adata,
                basis=basis,
                color=group_by,
                size=point_size,
                palette=palette,
                show=show_plot,
            )
            if save:
                plt.savefig(
                    ".".join(filter(None, [fileprefix, "palantir_group_by.pdf"])),
                    dpi=dpi,
                )

    finally:
        os.chdir(prevdir)

    try:
        adata.__dict__["_raw"].__dict__["_var"] = (
            adata.__dict__["_raw"]
            .__dict__["_var"]
            .rename(columns={"_index": "features"})
        )
    except:
        pass

    return adata


def WOT(
    adata=None,
    h5ad=None,
    group_by=None,
    palette=None,
    time_field="Time",
    growth_iters=3,
    tmap_out="tmaps/tmap_out",
    time_from=None,
    time_to=None,
    get_coupling=False,
    recalculate=False,
    show_plot=True,
    dpi=300,
    save=False,
    dirpath="./",
    fileprefix="",
    verbose=True,
):
    # Configure OpenMP settings to prevent conflicts
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    # Check if running on Apple silicon
    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"

    if is_apple_silicon:
        log_message(
            "Apple silicon detected: Applying specific configurations",
            message_type="info",
            verbose=verbose,
        )

        os.environ["PYTHONHASHSEED"] = "0"
        os.environ["PYTHONUNBUFFERED"] = "1"
        os.environ["SCANPY_SETTINGS"] = "scanpy_settings"
        os.environ["MPLBACKEND"] = "Agg"
        os.environ["DISPLAY"] = ""

        import numba

        numba.config.DISABLE_JIT = True  # Disable JIT compilation
        numba.set_num_threads(1)  # Force single thread
        log_message(
            "NUMBA configured for Apple silicon",
            message_type="success",
            verbose=verbose,
        )

    import matplotlib.pyplot as plt
    import scanpy as sc
    import numpy as np
    import statistics
    import pandas as pd
    from math import hypot
    import wot

    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    import os

    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(dirpath))

    import platform

    if platform.system() == "Windows":
        import sys
        import re
        import multiprocessing

        if re.match(pattern=".*pythonw.exe$", string=sys.executable):
            pythonw = sys.executable
        else:
            pythonw = sys.executable.replace("python.exe", "pythonw.exe")
        sys.executable = pythonw
        sys._base_executable = pythonw
        multiprocessing.set_executable(pythonw)

    try:
        if adata is None and h5ad is None:
            log_message("adata or h5ad must be provided.", message_type="error")
            exit()

        if adata is None:
            adata = sc.read(h5ad)

        if group_by is None:
            log_message("group_by must be provided.", message_type="error")
            exit()

        if time_field is None:
            log_message("time_field must be provided.", message_type="error")
            exit()

        adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")
        if pd.api.types.is_categorical_dtype(adata.obs[time_field]):
            adata.obs["time_field"] = adata.obs[time_field].cat.codes
        elif not pd.api.types.is_numeric_dtype(adata.obs[time_field]):
            try:
                adata.obs["time_field"] = adata.obs[time_field].astype("float")
            except ValueError:
                log_message(
                    f"Unable to convert column '{time_field}' to float type.",
                    message_type="warning",
                )
        else:
            adata.obs["time_field"] = adata.obs[time_field]

        time_dict = dict(zip(adata.obs[time_field], adata.obs["time_field"]))
        if time_from not in time_dict.keys():
            log_message("'time_from' is incorrect", message_type="error")
            exit()

        ot_model = wot.ot.OTModel(
            adata, growth_iters=growth_iters, day_field="time_field"
        )

        if recalculate is True:
            ot_model.compute_all_transport_maps(tmap_out=tmap_out)
            tmap_model = wot.tmap.TransportMapModel.from_directory(tmap_out)
        else:
            try:
                tmap_model = wot.tmap.TransportMapModel.from_directory(tmap_out)
            except (FileNotFoundError, ValueError):
                ot_model.compute_all_transport_maps(tmap_out=tmap_out)
                tmap_model = wot.tmap.TransportMapModel.from_directory(tmap_out)

        cell_sets = {}
        for k, v in zip(adata.obs[group_by], adata.obs_names):
            if k not in cell_sets:
                cell_sets[k] = []
            cell_sets[k].append(v)

        from_populations = tmap_model.population_from_cell_sets(
            cell_sets, at_time=time_dict[time_from]
        )

        trajectory_ds = tmap_model.trajectories(from_populations)
        trajectory_df = pd.DataFrame(
            trajectory_ds.X,
            index=trajectory_ds.obs_names,
            columns=trajectory_ds.var_names,
        )
        adata.uns["trajectory_" + str(time_from)] = trajectory_df.reindex(
            adata.obs_names
        )

        fates_ds = tmap_model.fates(from_populations)
        fates_df = pd.DataFrame(
            fates_ds.X, index=fates_ds.obs_names, columns=fates_ds.var_names
        )
        existing_rows = fates_df.index.tolist()
        new_rows = list(set(adata.obs_names) - set(existing_rows))
        new_df = pd.DataFrame(0, index=new_rows, columns=fates_df.columns)
        fates_df = pd.concat([fates_df, new_df])
        adata.uns["fates_" + str(time_from)] = fates_df.reindex(adata.obs_names)

        # obs_list = wot.tmap.trajectory_trends_from_trajectory(
        #     trajectory_ds = trajectory_ds,
        #     expression_ds = adata
        # )

        if time_to is not None:
            if time_to not in time_dict.keys():
                log_message("`time_to` is incorrect", message_type="error")
                exit()

            to_populations = tmap_model.population_from_cell_sets(
                cell_sets, at_time=time_dict[time_to]
            )
            transition_table = tmap_model.transition_table(
                from_populations, to_populations
            )
            transition_df = pd.DataFrame(
                transition_table.X,
                index=transition_table.obs_names,
                columns=transition_table.var_names,
            )
            adata.uns["transition_" + str(time_from) + "_to_" + str(time_to)] = fates_df
            if get_coupling is True:
                coupling = tmap_model.get_coupling(
                    time_dict[time_from], time_dict[time_to]
                )
                coupling_df = pd.DataFrame(
                    coupling.X, index=coupling.obs_names, columns=coupling.var_names
                )
                adata.uns["coupling_" + str(time_from) + "_to_" + str(time_to)] = (
                    coupling_df
                )

    finally:
        os.chdir(prevdir)

    try:
        adata.__dict__["_raw"].__dict__["_var"] = (
            adata.__dict__["_raw"]
            .__dict__["_var"]
            .rename(columns={"_index": "features"})
        )
    except:
        pass

    return adata
