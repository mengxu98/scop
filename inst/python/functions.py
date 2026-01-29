# LogMessage class and log_message function

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
    save_plot=False,
    plot_format="png",
    plot_dpi=600,
    plot_prefix="scvelo",
    dirpath="./scvelo",
    save=False,
    dpi=300,
    fileprefix="",
    verbose=True,
    legend_loc="on data",
):
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

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

        numba.config.DISABLE_JIT = True
        numba.set_num_threads(1)
        log_message(
            "NUMBA configured for Apple silicon",
            message_type="success",
            verbose=verbose,
        )

    import matplotlib

    matplotlib.use("Agg")
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

    if save_plot:
        expanded_path = os.path.expanduser(dirpath)
        if not os.path.exists(expanded_path):
            os.makedirs(expanded_path, exist_ok=True)
            log_message(
                f"Created directory: {expanded_path}",
                message_type="info",
                verbose=verbose,
            )
        os.chdir(expanded_path)
        sc.settings.figdir = "."
    else:
        sc.settings.figdir = "."
    sc.settings.file_format_figs = plot_format
    sc.settings.autosave = save_plot

    log_message(
        "Starting {.pkg scVelo} analysis...", message_type="running", verbose=verbose
    )
    try:
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

        if basis is None:
            if nonlinear_reduction is not None:
                if nonlinear_reduction in adata.obsm:
                    basis = nonlinear_reduction
                else:
                    log_message(
                        "{.val {nonlinear_reduction}} not found in adata.obsm. Available keys: {.val {list(adata.obsm.keys())}}",
                        message_type="error",
                        verbose=verbose,
                    )
                    raise ValueError(
                        f"nonlinear_reduction '{nonlinear_reduction}' not found in adata.obsm"
                    )
            else:
                if linear_reduction in adata.obsm:
                    basis = linear_reduction
                else:
                    log_message(
                        "{.val {linear_reduction}} not found in adata.obsm. Available keys: {.val {list(adata.obsm.keys())}}",
                        message_type="error",
                        verbose=verbose,
                    )
                    raise ValueError(
                        f"linear_reduction '{linear_reduction}' not found in adata.obsm"
                    )

        if basis not in adata.obsm:
            log_message(
                "basis '{.val {basis}}' not found in adata.obsm. Available keys: {.val {list(adata.obsm.keys())}}",
                message_type="error",
                verbose=verbose,
            )
            if linear_reduction in adata.obsm:
                adata.obsm["basis"] = adata.obsm[linear_reduction][:, 0:2]
                basis = "basis"
            else:
                raise ValueError(
                    f"Cannot find suitable basis. Available obsm keys: {list(adata.obsm.keys())}"
                )

        log_message("Using basis: {.val {basis}}", message_type="info", verbose=verbose)
        log_message(
            "Available embeddings in adata.obsm: {.val {list(adata.obsm.keys())}}",
            message_type="info",
            verbose=verbose,
        )

        adata.obs[group_by] = adata.obs[group_by].astype("category")

        log_message("Starting preprocessing", message_type="running", verbose=verbose)

        if filter_genes:
            log_message("Filtering genes...", message_type="info", verbose=verbose)
            scv.pp.filter_genes(adata, min_counts=min_counts)
            scv.pp.filter_genes(adata, min_counts_u=min_counts_u)

        if normalize_per_cell:
            log_message("Normalizing per cell...", message_type="info", verbose=verbose)
            scv.pp.normalize_per_cell(adata)

        if log_transform:
            log_message("Log transforming...", message_type="info", verbose=verbose)
            sc.pp.log1p(adata)

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

        log_message(
            "Computing {.pkg scvelo} neighbors and moments...",
            message_type="info",
            verbose=verbose,
        )

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

        if use_rep is None:
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
                use_rep = None
                n_pcs = min(n_pcs, adata.X.shape[1])
                log_message(
                    "Using raw data with n_pcs={.val {n_pcs}}",
                    message_type="info",
                    verbose=verbose,
                )

        try:
            scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors, use_rep=use_rep)
        except Exception as e:
            log_message(
                "{.pkg scvelo} moments failed ({.val {e}}), using manual computation...",
                message_type="warning",
                verbose=verbose,
            )
            sc.pp.neighbors(
                adata, n_pcs=n_pcs, n_neighbors=n_neighbors, use_rep=use_rep
            )

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
            if basis not in adata.obsm:
                log_message(
                    f"Basis '{basis}' not found in adata.obsm. Available keys: {list(adata.obsm.keys())}",
                    message_type="error",
                    verbose=verbose,
                )
                raise ValueError(
                    f"Basis '{basis}' not found in adata.obsm. Available keys: {list(adata.obsm.keys())}"
                )

            basis_embedding = adata.obsm[basis]
            if basis_embedding.shape[1] < 2:
                log_message(
                    f"Basis '{basis}' has only {basis_embedding.shape[1]} dimensions, need at least 2 for velocity embedding",
                    message_type="error",
                    verbose=verbose,
                )
                raise ValueError(
                    f"Basis '{basis}' must have at least 2 dimensions for velocity embedding"
                )

            x_basis_key = f"X_{basis}"
            if x_basis_key not in adata.obsm:
                adata.obsm[x_basis_key] = basis_embedding
                log_message(
                    f"Created '{x_basis_key}' in adata.obsm for scvelo compatibility",
                    message_type="info",
                    verbose=verbose,
                )

            scv.tl.velocity_embedding(adata, basis=basis, vkey=m)

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

            if compute_terminal_states:
                log_message(
                    "Computing terminal states...", message_type="info", verbose=verbose
                )
                try:
                    scv.tl.terminal_states(adata, vkey=m)
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

            if compute_paga:
                log_message("Computing PAGA...", message_type="info", verbose=verbose)
                try:
                    if "neighbors" not in adata.uns:
                        adata.uns["neighbors"] = {}
                    adata.uns["neighbors"]["distances"] = adata.obsp["distances"]
                    adata.uns["neighbors"]["connectivities"] = adata.obsp[
                        "connectivities"
                    ]

                    if m + "_graph" in adata.uns:
                        adata.uns["velocity_graph"] = adata.uns[m + "_graph"]

                    adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")

                    sc.tl.paga(adata, groups=group_by, use_rna_velocity=True)
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

                groups = (
                    adata.obs[group_by].cat.categories
                    if hasattr(adata.obs[group_by], "cat")
                    else adata.obs[group_by].unique()
                )
                if palette is None:
                    palette = dict(
                        zip(groups, plt.cm.tab10(np.linspace(0, 1, len(groups))))
                    )

                def get_scvelo_save_path(name):
                    if save_plot:
                        ext = (
                            plot_format
                            if plot_format in ["png", "pdf", "svg"]
                            else "png"
                        )
                        save_path = os.path.abspath(f"{plot_prefix}_{m}_{name}.{ext}")
                        return save_path
                    return False

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
                        legend_loc=legend_loc,
                        save=False,
                        show=show_plot,
                        dpi=plot_dpi,
                    )
                    if save_plot:
                        save_path = get_scvelo_save_path("stream")
                        try:
                            plt.savefig(save_path, dpi=plot_dpi, bbox_inches="tight")
                        except Exception as pdf_error:
                            if os.path.exists(save_path):
                                try:
                                    os.remove(save_path)
                                except Exception:
                                    pass
                            save_path_png = save_path.replace(".pdf", ".png")
                            plt.savefig(
                                save_path_png, dpi=plot_dpi, bbox_inches="tight"
                            )
                            log_message(
                                f"PDF save failed, saved as PNG instead: {save_path_png}",
                                message_type="warning",
                                verbose=verbose,
                            )
                    plt.close()
                except Exception as e:
                    log_message(
                        "stream plot failed ({.val {e}})",
                        message_type="warning",
                        verbose=verbose,
                    )
                    if plt is not None:
                        plt.close()

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
                        legend_loc=legend_loc,
                        save=False,
                        show=show_plot,
                        dpi=plot_dpi,
                    )
                    if save_plot:
                        save_path = get_scvelo_save_path("arrow")
                        try:
                            plt.savefig(save_path, dpi=plot_dpi, bbox_inches="tight")
                        except Exception as pdf_error:
                            if os.path.exists(save_path):
                                try:
                                    os.remove(save_path)
                                except Exception:
                                    pass
                            save_path_png = save_path.replace(".pdf", ".png")
                            plt.savefig(
                                save_path_png, dpi=plot_dpi, bbox_inches="tight"
                            )
                            log_message(
                                f"PDF save failed, saved as PNG instead: {save_path_png}",
                                message_type="warning",
                                verbose=verbose,
                            )
                    plt.close()
                except Exception as e:
                    log_message(
                        "arrow plot failed ({.val {e}})",
                        message_type="warning",
                        verbose=verbose,
                    )
                    if plt is not None:
                        plt.close()

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
                                    legend_loc="right margin",
                                    save=get_scvelo_save_path(metric),
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
        try:
            figures_dir = os.path.join(os.getcwd(), "figures")
            if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
                if not os.listdir(figures_dir):
                    os.rmdir(figures_dir)
                    log_message(
                        f"Removed empty figures directory: {figures_dir}",
                        message_type="info",
                        verbose=verbose,
                    )
        except Exception:
            pass

        if save_plot:
            os.chdir(prevdir)

    try:
        if hasattr(adata, "_raw") and adata._raw is not None:
            if hasattr(adata._raw, "_var"):
                adata._raw._var = adata._raw._var.rename(columns={"_index": "features"})
    except Exception:
        pass

    return adata


def compute_transition_matrix(kernel, verbose=True):
    """
    Compute transition matrix

    Parameters
    ----------
    kernel : cellrank.kernels.Kernel
        The kernel object
    verbose : bool
        Whether to log messages

    Returns
    -------
    bool
        True if matrix was computed, False otherwise
    """
    try:
        kernel.compute_transition_matrix()
        return fix_transition_matrix(kernel, verbose=verbose)
    except ValueError as e:
        if "not row stochastic" in str(e):
            if verbose:
                log_message(
                    f"Transition matrix validation failed: {e}. Attempting to fix...",
                    message_type="warning",
                    verbose=verbose,
                )
            try:
                import numpy as np
                from scipy.sparse import issparse, csr_matrix

                if hasattr(kernel, "connectivity"):
                    conn = kernel.connectivity
                elif hasattr(kernel, "_conn"):
                    conn = kernel._conn
                else:
                    if (
                        hasattr(kernel, "adata")
                        and "connectivities" in kernel.adata.obsp
                    ):
                        conn = kernel.adata.obsp["connectivities"]
                    else:
                        raise ValueError("Cannot access connectivity matrix")

                if issparse(conn):
                    row_sums = np.array(conn.sum(axis=1)).flatten()
                    row_sums_safe = np.maximum(row_sums, 1e-10)
                    from scipy.sparse import diags

                    tmat = diags(1.0 / row_sums_safe) @ conn
                else:
                    row_sums = conn.sum(axis=1, keepdims=True)
                    row_sums_safe = np.maximum(row_sums, 1e-10)
                    tmat = conn / row_sums_safe

                kernel._transition_matrix = tmat

                fix_transition_matrix(kernel, verbose=verbose)

                if verbose:
                    log_message(
                        "Transition matrix fixed successfully",
                        message_type="success",
                        verbose=verbose,
                    )
                return True
            except Exception as fix_error:
                if verbose:
                    log_message(
                        f"Failed to fix transition matrix: {fix_error}",
                        message_type="error",
                        verbose=verbose,
                    )
                raise
        else:
            raise


def fix_transition_matrix(kernel, verbose=True):
    """
    Fix transition matrix to be row stochastic (rows sum to 1).

    This function ensures the transition matrix is valid for CellRank estimators
    by cleaning NaN/Inf values, clipping negatives, adding self-loops where needed,
    and normalizing rows to sum to 1.

    Parameters
    ----------
    kernel : cellrank.kernels.Kernel
        The kernel object with a transition_matrix attribute
    verbose : bool
        Whether to log messages

    Returns
    -------
    bool
        True if matrix was modified, False otherwise
    """
    try:
        import numpy as np
        from scipy.sparse import diags, issparse, csr_matrix

        tmat = kernel.transition_matrix
        matrix_modified = False

        if verbose:
            log_message(
                "Validating and fixing transition matrix...",
                message_type="info",
                verbose=verbose,
            )

        if issparse(tmat):
            if np.any(np.isnan(tmat.data)) or np.any(np.isinf(tmat.data)):
                if verbose:
                    log_message(
                        "Cleaning NaN/Inf values...",
                        message_type="warning",
                        verbose=verbose,
                    )
                tmat.data[np.isnan(tmat.data)] = 0
                tmat.data[np.isinf(tmat.data)] = 0
                tmat.eliminate_zeros()
                matrix_modified = True
        else:
            if np.any(np.isnan(tmat)) or np.any(np.isinf(tmat)):
                if verbose:
                    log_message(
                        "Cleaning NaN/Inf values...",
                        message_type="warning",
                        verbose=verbose,
                    )
                tmat[np.isnan(tmat)] = 0
                tmat[np.isinf(tmat)] = 0
                matrix_modified = True

        if issparse(tmat):
            if np.any(tmat.data < 0):
                if verbose:
                    log_message(
                        f"Clipping {(tmat.data < 0).sum()} negative values to zero...",
                        message_type="warning",
                        verbose=verbose,
                    )
                tmat.data = np.maximum(tmat.data, 0)
                matrix_modified = True
        else:
            if np.any(tmat < 0):
                if verbose:
                    log_message(
                        "Clipping negative values to zero...",
                        message_type="warning",
                        verbose=verbose,
                    )
                tmat = np.maximum(tmat, 0)
                matrix_modified = True

        if issparse(tmat):
            diag_values = tmat.diagonal()
            min_self_loop = 0.01
            needs_self_loop = diag_values < min_self_loop

            if np.any(needs_self_loop):
                if verbose:
                    log_message(
                        f"Adding self-loops to {needs_self_loop.sum()} cells for matrix primitivity...",
                        message_type="info",
                        verbose=verbose,
                    )

                tmat = tmat.tolil()
                for i in np.where(needs_self_loop)[0]:
                    tmat[i, i] = min_self_loop
                tmat = tmat.tocsr()
                matrix_modified = True

        row_sums = np.array(tmat.sum(axis=1)).flatten()
        zero_rows = row_sums < 1e-10

        if np.any(zero_rows):
            if verbose:
                log_message(
                    f"Found {zero_rows.sum()} zero rows. Setting to self-loops...",
                    message_type="warning",
                    verbose=verbose,
                )

            if issparse(tmat):
                tmat = tmat.tolil()
                for i in np.where(zero_rows)[0]:
                    tmat[i, i] = 1.0
                tmat = tmat.tocsr()
            else:
                for i in np.where(zero_rows)[0]:
                    tmat[i, :] = 0
                    tmat[i, i] = 1.0

            row_sums = np.array(tmat.sum(axis=1)).flatten()
            matrix_modified = True

        if not np.allclose(row_sums, 1.0, rtol=1e-6, atol=1e-8):
            log_message(
                f"Normalizing rows (current range: {row_sums.min():.8f} - {row_sums.max():.8f})...",
                message_type="info",
                verbose=verbose,
            )

            row_sums_safe = np.maximum(row_sums, 1e-10)

            if issparse(tmat):
                tmat = diags(1.0 / row_sums_safe) @ tmat
            else:
                tmat = tmat / row_sums_safe[:, np.newaxis]

            matrix_modified = True

        final_row_sums = np.array(tmat.sum(axis=1)).flatten()

        if issparse(tmat):
            has_negative = np.any(tmat.data < -1e-10)
        else:
            has_negative = np.any(tmat < -1e-10)

        rows_sum_to_one = np.allclose(final_row_sums, 1.0, rtol=1e-6, atol=1e-8)

        if has_negative and verbose:
            log_message(
                "Warning: Matrix still has negative values after fixing",
                message_type="warning",
                verbose=verbose,
            )

        if not rows_sum_to_one and verbose:
            log_message(
                f"Warning: Matrix rows still don't sum to 1 (range: {final_row_sums.min():.8f} - {final_row_sums.max():.8f})",
                message_type="warning",
                verbose=verbose,
            )

        if matrix_modified:
            kernel._transition_matrix = tmat
            if verbose:
                log_message(
                    f"Matrix fixed and validated (row sums: {final_row_sums.min():.8f} - {final_row_sums.max():.8f})",
                    message_type="success",
                    verbose=verbose,
                )
        elif verbose:
            log_message(
                f"Matrix validation passed (row sums: {final_row_sums.min():.8f} - {final_row_sums.max():.8f})",
                message_type="success",
                verbose=verbose,
            )

        return matrix_modified

    except Exception as e:
        if verbose:
            log_message(
                f"Matrix validation encountered error: {e}. Proceeding with original matrix...",
                message_type="warning",
                verbose=verbose,
            )
        return False


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
    dirpath="./cellrank",
    fileprefix="",
    verbose=True,
    kernel_type="velocity",
    time_key="dpt_pseudotime",
    estimator_type="GPCCA",
    use_connectivity_kernel=True,
    velocity_weight=0.8,
    connectivity_weight=0.2,
    softmax_scale=4,
    n_macrostates=None,
    schur_method="krylov",
    n_cells_terminal=10,
    save_plot=False,
    plot_format="png",
    plot_dpi=600,
    plot_prefix="cellrank",
    legend_loc="on data",
):
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

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

        numba.config.DISABLE_JIT = True
        numba.set_num_threads(1)
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

    if save_plot:
        expanded_path = os.path.expanduser(dirpath)
        if not os.path.exists(expanded_path):
            os.makedirs(expanded_path, exist_ok=True)
            log_message(
                f"Created directory: {expanded_path}",
                message_type="info",
                verbose=verbose,
            )
        os.chdir(expanded_path)
        sc.settings.figdir = "."
        cr.settings.figdir = "."
    else:
        sc.settings.figdir = "."
        cr.settings.figdir = "."
    sc.settings.file_format_figs = plot_format
    sc.settings.autosave = save_plot
    log_message(
        f"CellRank figdir set to: {cr.settings.figdir}",
        message_type="info",
        verbose=verbose,
    )

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
            log_message("adata or h5ad must be provided", message_type="error")
            exit()

        if adata is None:
            adata = scv.read(h5ad)

        if group_by is None:
            log_message("group_by must be provided", message_type="error")
            exit()

        if linear_reduction is None and nonlinear_reduction is None:
            log_message(
                "linear_reduction or nonlinear_reduction must be provided at least one",
                message_type="error",
            )
            exit()

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

        log_message(
            f"Using kernel_type: '{kernel_type}'", message_type="info", verbose=verbose
        )

        use_velocity = False
        use_pseudotime = False
        use_cytotrace = False

        if kernel_type == "velocity":
            has_velocity_data = (
                "spliced" in adata.layers and "unspliced" in adata.layers
            )

            if not has_velocity_data:
                log_message(
                    "No spliced/unspliced data found. Consider using kernel_type='cytotrace' or 'pseudotime' for RNA-only data",
                    message_type="warning",
                    verbose=verbose,
                )
                log_message(
                    "Falling back to ConnectivityKernel only",
                    message_type="info",
                    verbose=verbose,
                )
                use_velocity = False
            else:
                use_velocity = True

                if mode[-1] + "_graph" not in adata.obs.keys():
                    log_message(
                        "Running scVelo to compute RNA velocity...",
                        message_type="info",
                        verbose=verbose,
                    )
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
                        save_plot=save_plot,
                        plot_format=plot_format,
                        plot_dpi=plot_dpi,
                        plot_prefix="scvelo",
                        dirpath=".",
                        legend_loc=legend_loc,
                        dpi=dpi,
                        save=save,
                        fileprefix=fileprefix,
                    )
                adata.layers["velocity"] = adata.layers[mode[-1]]

        elif kernel_type == "pseudotime":
            use_pseudotime = True
            if time_key not in adata.obs:
                log_message(
                    f"Pseudotime '{time_key}' not found. Computing DPT pseudotime...",
                    message_type="info",
                    verbose=verbose,
                )
                if linear_reduction and linear_reduction in adata.obsm:
                    rep_key = linear_reduction
                else:
                    log_message(
                        f"linear_reduction '{linear_reduction}' not found in adata.obsm",
                        message_type="error",
                        verbose=verbose,
                    )
                    exit()

                if "connectivities" not in adata.obsp:
                    sc.pp.neighbors(
                        adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep_key
                    )

                sc.tl.diffmap(adata)
                adata.uns["iroot"] = 0
                sc.tl.dpt(adata)
                time_key = "dpt_pseudotime"
                adata.obs["cellrank_pseudotime"] = adata.obs["dpt_pseudotime"]
                log_message(
                    "DPT pseudotime computed and stored in adata.obs['cellrank_pseudotime']",
                    message_type="success",
                    verbose=verbose,
                )
            else:
                log_message(
                    f"Using existing pseudotime from adata.obs['{time_key}']",
                    message_type="info",
                    verbose=verbose,
                )

        elif kernel_type == "cytotrace":
            use_cytotrace = True
            log_message(
                "Using CytoTRACEKernel for RNA-only data...",
                message_type="info",
                verbose=verbose,
            )

        log_message(
            "Using CellRank kernel-estimator architecture...",
            message_type="info",
            verbose=verbose,
        )

        main_kernel = None

        def get_rep_key():
            if linear_reduction and linear_reduction in adata.obsm:
                return linear_reduction
            else:
                log_message(
                    f"linear_reduction '{linear_reduction}' not found in adata.obsm",
                    message_type="error",
                    verbose=verbose,
                )
                exit()

        if use_pseudotime:
            log_message(
                f"Creating PseudotimeKernel with time_key='{time_key}'...",
                message_type="info",
                verbose=verbose,
            )
            try:
                pk = cr.kernels.PseudotimeKernel(adata, time_key=time_key)
                compute_transition_matrix(pk, verbose=verbose)
                main_kernel = pk
                log_message(
                    "PseudotimeKernel created successfully",
                    message_type="success",
                    verbose=verbose,
                )
            except Exception as e:
                log_message(
                    f"PseudotimeKernel failed: {e}. Falling back to ConnectivityKernel.",
                    message_type="warning",
                    verbose=verbose,
                )
                use_pseudotime = False

        elif use_cytotrace:
            log_message(
                "Creating CytoTRACEKernel and computing CytoTRACE score...",
                message_type="info",
                verbose=verbose,
            )
            try:
                if "Ms" not in adata.layers:
                    log_message(
                        "Computing moments (required for CytoTRACEKernel)...",
                        message_type="info",
                        verbose=verbose,
                    )
                    rep_key = get_rep_key()
                    if "connectivities" not in adata.obsp:
                        sc.pp.neighbors(
                            adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep_key
                        )

                    import numpy as np
                    from scipy.sparse import issparse, csr_matrix

                    log_message(
                        "Creating smoothed expression layer (Ms) for RNA-only data...",
                        message_type="info",
                        verbose=verbose,
                    )

                    connectivities = adata.obsp["connectivities"].copy()
                    if issparse(connectivities):
                        connectivities = connectivities.toarray()
                    connectivities += np.eye(connectivities.shape[0])
                    row_sums = connectivities.sum(axis=1, keepdims=True)
                    connectivities = connectivities / row_sums

                    X = adata.X.toarray() if issparse(adata.X) else adata.X
                    Ms = connectivities @ X
                    adata.layers["Ms"] = csr_matrix(Ms)

                    log_message(
                        "Smoothed expression layer (Ms) created successfully",
                        message_type="success",
                        verbose=verbose,
                    )

                ctk = cr.kernels.CytoTRACEKernel(adata)
                ctk.compute_cytotrace()
                compute_transition_matrix(ctk, verbose=verbose)
                main_kernel = ctk

                if "ct_score" in adata.obs:
                    adata.obs["cytotrace_score"] = adata.obs.pop("ct_score")
                    adata.obs["cytotrace_pseudotime"] = 1 - adata.obs["cytotrace_score"]
                    adata.obs["cellrank_pseudotime"] = adata.obs["cytotrace_pseudotime"]
                if "ct_pseudotime" in adata.obs:
                    del adata.obs["ct_pseudotime"]
                if "ct_num_exp_genes" in adata.obs:
                    adata.obs["cytotrace_num_exp_genes"] = adata.obs.pop(
                        "ct_num_exp_genes"
                    )

                log_message(
                    "CytoTRACEKernel created successfully",
                    message_type="success",
                    verbose=verbose,
                )
            except Exception as e:
                log_message(
                    f"CytoTRACEKernel failed: {e}. Falling back to ConnectivityKernel.",
                    message_type="warning",
                    verbose=verbose,
                )
                use_cytotrace = False

        elif use_velocity:
            if velocity_weight <= 0 and connectivity_weight <= 0:
                log_message(
                    "Both kernel weights are <= 0. Using equal weights (0.5, 0.5).",
                    message_type="warning",
                    verbose=verbose,
                )
                velocity_weight = 0.5
                connectivity_weight = 0.5
            elif velocity_weight <= 0:
                log_message(
                    "velocity_weight <= 0. Using ConnectivityKernel only.",
                    message_type="info",
                    verbose=verbose,
                )
                use_velocity = False
                connectivity_weight = 1.0
            elif connectivity_weight <= 0:
                log_message(
                    "connectivity_weight <= 0. Using VelocityKernel only.",
                    message_type="info",
                    verbose=verbose,
                )
                use_connectivity_kernel = False
                velocity_weight = 1.0
            else:
                total_weight = velocity_weight + connectivity_weight
                if abs(total_weight - 1.0) > 0.01:
                    log_message(
                        f"Normalizing kernel weights ({velocity_weight}, {connectivity_weight}) to sum to 1.0",
                        message_type="info",
                        verbose=verbose,
                    )
                    velocity_weight = velocity_weight / total_weight
                    connectivity_weight = connectivity_weight / total_weight

        if use_velocity and velocity_weight > 0:
            log_message(
                f"Creating VelocityKernel with model='{mode[-1]}', softmax_scale={softmax_scale}...",
                message_type="info",
                verbose=verbose,
            )
            try:
                vk = cr.kernels.VelocityKernel(adata, backward=False)
                vk.compute_transition_matrix(
                    model=mode[-1], softmax_scale=softmax_scale
                )
                compute_transition_matrix(vk, verbose=verbose)
                main_kernel = vk
            except Exception as e:
                log_message(
                    f"VelocityKernel computation failed: {e}",
                    message_type="error",
                    verbose=verbose,
                )
                log_message(
                    "Trying with default parameters...",
                    message_type="info",
                    verbose=verbose,
                )
                try:
                    vk = cr.kernels.VelocityKernel(adata, backward=False)
                    vk.compute_transition_matrix(softmax_scale=softmax_scale)
                    compute_transition_matrix(vk, verbose=verbose)
                    main_kernel = vk
                except Exception as e2:
                    log_message(
                        f"VelocityKernel still failed: {e2}. Using ConnectivityKernel only.",
                        message_type="warning",
                        verbose=verbose,
                    )
                    use_velocity = False

        def ensure_neighbors():
            if "connectivities" not in adata.obsp:
                log_message(
                    "Neighbors not found. Computing neighbors...",
                    message_type="info",
                    verbose=verbose,
                )
                rep_key = get_rep_key()
                sc.pp.neighbors(
                    adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep_key
                )
                log_message(
                    f"Neighbors computed successfully (n_neighbors={n_neighbors})",
                    message_type="success",
                    verbose=verbose,
                )

        if (
            main_kernel is not None
            and use_connectivity_kernel
            and connectivity_weight > 0
        ):
            try:
                log_message(
                    "Creating ConnectivityKernel...",
                    message_type="info",
                    verbose=verbose,
                )
                ensure_neighbors()
                ck = cr.kernels.ConnectivityKernel(adata)
                compute_transition_matrix(ck, verbose=verbose)
                final_kernel = velocity_weight * main_kernel + connectivity_weight * ck
                log_message(
                    f"Combined kernels with weights: main={velocity_weight:.2f}, connectivity={connectivity_weight:.2f}",
                    message_type="success",
                    verbose=verbose,
                )
            except Exception as e:
                log_message(
                    f"Failed to combine kernels: {e}. Using main kernel only.",
                    message_type="warning",
                    verbose=verbose,
                )
                final_kernel = main_kernel
        elif main_kernel is not None:
            final_kernel = main_kernel
            log_message(
                f"Using {type(main_kernel).__name__} only",
                message_type="info",
                verbose=verbose,
            )
        else:
            log_message(
                "Creating ConnectivityKernel (fallback)...",
                message_type="info",
                verbose=verbose,
            )
            ensure_neighbors()
            ck = cr.kernels.ConnectivityKernel(adata)
            compute_transition_matrix(ck, verbose=verbose)
            final_kernel = ck
            log_message(
                "Using ConnectivityKernel only (connectivity-based transitions)",
                message_type="warning",
                verbose=verbose,
            )

        try:
            import numpy as np
            from scipy.sparse import diags, issparse, csr_matrix

            tmat = final_kernel.transition_matrix
            matrix_modified = False

            log_message(
                "Validating and fixing transition matrix for GPCCA compatibility...",
                message_type="info",
                verbose=verbose,
            )

            if issparse(tmat):
                if np.any(np.isnan(tmat.data)) or np.any(np.isinf(tmat.data)):
                    log_message(
                        "Cleaning NaN/Inf values...",
                        message_type="warning",
                        verbose=verbose,
                    )
                    tmat.data[np.isnan(tmat.data)] = 0
                    tmat.data[np.isinf(tmat.data)] = 0
                    tmat.eliminate_zeros()
                    matrix_modified = True
            else:
                if np.any(np.isnan(tmat)) or np.any(np.isinf(tmat)):
                    log_message(
                        "Cleaning NaN/Inf values...",
                        message_type="warning",
                        verbose=verbose,
                    )
                    tmat[np.isnan(tmat)] = 0
                    tmat[np.isinf(tmat)] = 0
                    matrix_modified = True

            if issparse(tmat):
                if np.any(tmat.data < 0):
                    log_message(
                        f"Clipping {(tmat.data < 0).sum()} negative values to zero...",
                        message_type="warning",
                        verbose=verbose,
                    )
                    tmat.data = np.maximum(tmat.data, 0)
                    matrix_modified = True
            else:
                if np.any(tmat < 0):
                    log_message(
                        "Clipping negative values to zero...",
                        message_type="warning",
                        verbose=verbose,
                    )
                    tmat = np.maximum(tmat, 0)
                    matrix_modified = True

            if issparse(tmat):
                diag_values = tmat.diagonal()
                min_self_loop = 0.01
                needs_self_loop = diag_values < min_self_loop

                if np.any(needs_self_loop):
                    log_message(
                        f"Adding self-loops to {needs_self_loop.sum()} cells for matrix primitivity...",
                        message_type="info",
                        verbose=verbose,
                    )

                    tmat = tmat.tolil()
                    for i in np.where(needs_self_loop)[0]:
                        tmat[i, i] = min_self_loop
                    tmat = tmat.tocsr()
                    matrix_modified = True

            row_sums = np.array(tmat.sum(axis=1)).flatten()
            zero_rows = row_sums < 1e-10

            if np.any(zero_rows):
                log_message(
                    f"Found {zero_rows.sum()} zero rows. Setting to self-loops...",
                    message_type="warning",
                    verbose=verbose,
                )

                if issparse(tmat):
                    tmat = tmat.tolil()
                    for i in np.where(zero_rows)[0]:
                        tmat[i, i] = 1.0
                    tmat = tmat.tocsr()
                else:
                    for i in np.where(zero_rows)[0]:
                        tmat[i, :] = 0
                        tmat[i, i] = 1.0

                row_sums = np.array(tmat.sum(axis=1)).flatten()
                matrix_modified = True

            if not np.allclose(row_sums, 1.0, rtol=1e-6, atol=1e-8):
                log_message(
                    f"Normalizing rows (current range: {row_sums.min():.8f} - {row_sums.max():.8f})...",
                    message_type="info",
                    verbose=verbose,
                )

                row_sums_safe = np.maximum(row_sums, 1e-10)

                if issparse(tmat):
                    tmat = diags(1.0 / row_sums_safe) @ tmat
                else:
                    tmat = tmat / row_sums_safe[:, np.newaxis]

                matrix_modified = True

            final_row_sums = np.array(tmat.sum(axis=1)).flatten()

            if issparse(tmat):
                has_negative = np.any(tmat.data < -1e-10)
            else:
                has_negative = np.any(tmat < -1e-10)

            rows_sum_to_one = np.allclose(final_row_sums, 1.0, rtol=1e-6, atol=1e-8)

            if has_negative:
                log_message(
                    "Warning: Matrix still has negative values after fixing",
                    message_type="warning",
                    verbose=verbose,
                )

            if not rows_sum_to_one:
                log_message(
                    f"Warning: Matrix rows still don't sum to 1 (range: {final_row_sums.min():.8f} - {final_row_sums.max():.8f})",
                    message_type="warning",
                    verbose=verbose,
                )

            if matrix_modified:
                final_kernel._transition_matrix = tmat
                log_message(
                    f"Matrix fixed and validated (row sums: {final_row_sums.min():.8f} - {final_row_sums.max():.8f})",
                    message_type="success",
                    verbose=verbose,
                )
            else:
                log_message(
                    f"Matrix validation passed (row sums: {final_row_sums.min():.8f} - {final_row_sums.max():.8f})",
                    message_type="success",
                    verbose=verbose,
                )

        except Exception as e:
            log_message(
                f"Matrix validation encountered error: {e}. Proceeding with original matrix...",
                message_type="warning",
                verbose=verbose,
            )

        log_message(
            f"Creating {estimator_type} estimator...",
            message_type="info",
            verbose=verbose,
        )

        gpcca_failed = False

        if estimator_type == "GPCCA":
            try:
                estimator = cr.estimators.GPCCA(final_kernel)

                log_message(
                    "Computing eigendecomposition...",
                    message_type="info",
                    verbose=verbose,
                )
                estimator.compute_eigendecomposition()

                if n_macrostates is None:
                    n_cells = adata.n_obs
                    if n_cells < 100:
                        n_states_schur = 5
                    elif n_cells < 500:
                        n_states_schur = 8
                    elif n_cells < 2000:
                        n_states_schur = 10
                    else:
                        n_states_schur = 15

                    log_message(
                        f"Auto-determined n_states={n_states_schur} for Schur (based on {n_cells} cells)",
                        message_type="info",
                        verbose=verbose,
                    )
                else:
                    n_states_schur = n_macrostates + 2
                    log_message(
                        f"Using n_states={n_states_schur} for Schur (n_macrostates={n_macrostates} + 2)",
                        message_type="info",
                        verbose=verbose,
                    )

                log_message(
                    f"Computing Schur decomposition (n_states={n_states_schur}, method='{schur_method}')...",
                    message_type="info",
                    verbose=verbose,
                )
                try:
                    estimator.compute_schur(n_states_schur, method=schur_method)
                except (ValueError, RuntimeError) as schur_error:
                    if "subspace_angles" in str(
                        schur_error
                    ) or "invariant subspace" in str(schur_error):
                        if schur_method == "brandts":
                            log_message(
                                f"Schur decomposition with '{schur_method}' method failed: {schur_error}",
                                message_type="warning",
                                verbose=verbose,
                            )
                            log_message(
                                "Trying with 'krylov' method instead...",
                                message_type="info",
                                verbose=verbose,
                            )
                            try:
                                estimator.compute_schur(n_states_schur, method="krylov")
                            except Exception as krylov_error:
                                log_message(
                                    f"Schur decomposition with 'krylov' method also failed: {krylov_error}",
                                    message_type="warning",
                                    verbose=verbose,
                                )
                                raise schur_error
                        else:
                            raise
                    else:
                        raise

                if n_macrostates is None:
                    n_macro = max(2, n_states_schur - 2)
                else:
                    n_macro = n_macrostates

                log_message(
                    f"Computing macrostates (n_states={n_macro}, n_cells={n_cells_terminal})...",
                    message_type="info",
                    verbose=verbose,
                )
                try:
                    estimator.compute_macrostates(
                        n_states=n_macro, n_cells=n_cells_terminal
                    )
                except ValueError as macro_error:
                    if "Discretizing leads to a cluster" in str(
                        macro_error
                    ) or "0 samples" in str(macro_error):
                        log_message(
                            f"Macrostates computation failed: {macro_error}",
                            message_type="warning",
                            verbose=verbose,
                        )
                        log_message(
                            "Trying with fewer macrostates...",
                            message_type="info",
                            verbose=verbose,
                        )
                        # Try with fewer macrostates
                        n_macro_reduced = max(2, n_macro - 2)
                        try:
                            estimator.compute_macrostates(
                                n_states=n_macro_reduced, n_cells=n_cells_terminal
                            )
                            log_message(
                                f"Macrostates computed successfully with n_states={n_macro_reduced}",
                                message_type="success",
                                verbose=verbose,
                            )
                        except Exception as e2:
                            log_message(
                                f"Reduced macrostates also failed: {e2}",
                                message_type="warning",
                                verbose=verbose,
                            )
                            log_message(
                                "Automatically switching to CFLARE estimator (more robust for problematic matrices)...",
                                message_type="warning",
                                verbose=verbose,
                            )
                            gpcca_failed = True
                            estimator_type = "CFLARE"
                            raise
                    else:
                        raise

                if not gpcca_failed:
                    log_message(
                        f"Setting terminal states (n_cells={n_cells_terminal})...",
                        message_type="info",
                        verbose=verbose,
                    )
                    estimator.set_terminal_states(n_cells=n_cells_terminal)

                log_message(
                    "GPCCA estimator completed successfully",
                    message_type="success",
                    verbose=verbose,
                )

            except (ValueError, RuntimeError) as e:
                if (
                    "not a transition matrix" in str(e).lower()
                    or "schur" in str(e).lower()
                    or "gpcca" in str(e).lower()
                    or "Discretizing leads to a cluster" in str(e)
                    or "0 samples" in str(e)
                    or "subspace_angles" in str(e).lower()
                    or "invariant subspace" in str(e).lower()
                ):
                    if not gpcca_failed:
                        log_message(
                            f"GPCCA estimator failed: {e}",
                            message_type="warning",
                            verbose=verbose,
                        )
                        log_message(
                            "Automatically switching to CFLARE estimator (more robust for problematic matrices)...",
                            message_type="warning",
                            verbose=verbose,
                        )
                        gpcca_failed = True
                        estimator_type = "CFLARE"
                    else:
                        # Already tried to switch, just raise
                        raise
                else:
                    log_message(
                        f"GPCCA failed with unexpected error: {e}",
                        message_type="error",
                        verbose=verbose,
                    )
                    raise

        if estimator_type == "CFLARE":
            estimator = cr.estimators.CFLARE(final_kernel)

            log_message(
                "Computing eigendecomposition...", message_type="info", verbose=verbose
            )
            estimator.compute_eigendecomposition()

            predict_method = "leiden"
            log_message(
                f"Predicting terminal states (use={n_macrostates}, method='{predict_method}')...",
                message_type="info",
                verbose=verbose,
            )

            try:
                estimator.predict(use=n_macrostates, method=predict_method)
            except Exception as e:
                log_message(
                    f"Prediction with '{predict_method}' failed: {e}. Trying 'kmeans'...",
                    message_type="warning",
                    verbose=verbose,
                )
                estimator.predict(use=n_macrostates, method="kmeans")

            if gpcca_failed:
                log_message(
                    "Successfully switched to CFLARE estimator",
                    message_type="success",
                    verbose=verbose,
                )

        log_message(
            "Computing fate probabilities...", message_type="info", verbose=verbose
        )
        estimator.compute_fate_probabilities()

        log_message(
            f"Computing lineage drivers for cluster_key='{group_by}'...",
            message_type="info",
            verbose=verbose,
        )
        try:
            estimator.compute_lineage_drivers(cluster_key=group_by, use_raw=False)
        except RuntimeError as e:
            if "Compute `.fate_probabilities`" in str(e):
                log_message(
                    "Skipping lineage drivers computation (fate probabilities not available)",
                    message_type="warning",
                    verbose=verbose,
                )
            else:
                raise

        log_message(
            "CellRan analysis completed!", message_type="success", verbose=verbose
        )

        if "cellrank_pseudotime" not in adata.obs:
            if "to_terminal_states" in adata.obsm or "lineages_fwd" in adata.obsm:
                try:
                    fate_probs_key = (
                        "to_terminal_states"
                        if "to_terminal_states" in adata.obsm
                        else "lineages_fwd"
                    )
                    fate_probs = adata.obsm[fate_probs_key]
                    if hasattr(fate_probs, "X"):
                        import numpy as np

                        max_probs = np.array(fate_probs.X.max(axis=1)).flatten()
                    else:
                        max_probs = np.array(fate_probs.max(axis=1)).flatten()
                    adata.obs["cellrank_pseudotime"] = max_probs
                    log_message(
                        "Created cellrank_pseudotime from fate probabilities",
                        message_type="info",
                        verbose=verbose,
                    )
                except Exception as e:
                    log_message(
                        f"Failed to create cellrank_pseudotime from fate probabilities: {e}",
                        message_type="warning",
                        verbose=verbose,
                    )

        if save_plot:
            log_message(
                "Generating visualizations...", message_type="info", verbose=verbose
            )

            def get_save_path(name):
                ext = plot_format if plot_format in ["png", "pdf", "svg"] else "png"
                save_path = os.path.abspath(f"{plot_prefix}_{name}.{ext}")
                log_message(
                    f"Save path for {name}: {save_path}",
                    message_type="info",
                    verbose=verbose,
                )
                return save_path

            try:
                log_message(
                    "Plotting eigenvalue spectrum...",
                    message_type="info",
                    verbose=verbose,
                )
                estimator.plot_spectrum(
                    real_only=True, dpi=plot_dpi, save=get_save_path("spectrum")
                )
                if show_plot:
                    plt.show()
                plt.close()
            except Exception as e:
                log_message(
                    f"Failed to plot spectrum: {e}",
                    message_type="warning",
                    verbose=verbose,
                )

            if estimator_type == "GPCCA":
                try:
                    log_message(
                        "Plotting Schur matrix...",
                        message_type="info",
                        verbose=verbose,
                    )
                    estimator.plot_schur_matrix(
                        dpi=plot_dpi, save=get_save_path("schur_matrix")
                    )
                    if show_plot:
                        plt.show()
                    plt.close()
                except Exception as e:
                    log_message(
                        f"Failed to plot Schur matrix: {e}",
                        message_type="warning",
                        verbose=verbose,
                    )

                try:
                    log_message(
                        "Plotting coarse-grained transition matrix...",
                        message_type="info",
                        verbose=verbose,
                    )
                    estimator.plot_coarse_T(
                        show_initial_dist=True,
                        show_stationary_dist=True,
                        dpi=plot_dpi,
                        save=get_save_path("coarse_T"),
                    )
                    if show_plot:
                        plt.show()
                    plt.close()
                except Exception as e:
                    log_message(
                        f"Failed to plot coarse T: {e}",
                        message_type="warning",
                        verbose=verbose,
                    )

            try:
                log_message(
                    "Plotting all macrostates...",
                    message_type="info",
                    verbose=verbose,
                )
                estimator.plot_macrostates(
                    which="all",
                    basis=basis,
                    dpi=plot_dpi,
                    save=get_save_path("macrostates_all"),
                )
                if show_plot:
                    plt.show()
                plt.close()

                log_message(
                    "Plotting terminal states...",
                    message_type="info",
                    verbose=verbose,
                )
                estimator.plot_macrostates(
                    which="terminal",
                    basis=basis,
                    dpi=plot_dpi,
                    save=get_save_path("macrostates_terminal"),
                )
                if show_plot:
                    plt.show()
                plt.close()
            except Exception as e:
                log_message(
                    f"Failed to plot macrostates: {e}",
                    message_type="warning",
                    verbose=verbose,
                )

            try:
                log_message(
                    "Plotting fate probabilities...",
                    message_type="info",
                    verbose=verbose,
                )
                estimator.plot_fate_probabilities(
                    basis=basis,
                    ncols=3,
                    dpi=plot_dpi,
                    save=get_save_path("fate_probabilities"),
                )
                if show_plot:
                    plt.show()
                plt.close()
            except Exception as e:
                log_message(
                    f"Failed to plot fate probabilities: {e}",
                    message_type="warning",
                    verbose=verbose,
                )

            try:
                fate_probs_key = (
                    "to_terminal_states"
                    if not hasattr(estimator, "backward") or not estimator.backward
                    else "from_terminal_states"
                )
                if fate_probs_key in adata.obsm:
                    lineages = adata.obsm[fate_probs_key].names
                    log_message(
                        f"Plotting lineage drivers for {len(lineages)} lineages...",
                        message_type="info",
                        verbose=verbose,
                    )

                    for lineage in lineages:
                        try:
                            estimator.plot_lineage_drivers(
                                lineage=lineage,
                                n_genes=8,
                                dpi=plot_dpi,
                                save=get_save_path(f"lineage_drivers_{lineage}"),
                            )
                            if show_plot:
                                plt.show()
                            plt.close()
                        except Exception as e:
                            log_message(
                                f"Failed to plot lineage drivers for {lineage}: {e}",
                                message_type="warning",
                                verbose=verbose,
                            )
            except Exception as e:
                log_message(
                    f"Failed to plot lineage drivers: {e}",
                    message_type="warning",
                    verbose=verbose,
                )

            try:
                log_message(
                    "Plotting aggregate fate probabilities (mode=paga_pie)...",
                    message_type="info",
                    verbose=verbose,
                )
                cr.pl.aggregate_fate_probabilities(
                    adata,
                    mode="paga_pie",
                    cluster_key=group_by,
                    basis=basis,
                    dpi=plot_dpi,
                    save=get_save_path("aggregate_fates_paga_pie"),
                )
                if show_plot:
                    plt.show()
                plt.close()
            except Exception as e:
                log_message(
                    f"Failed to plot aggregate fates: {e}",
                    message_type="warning",
                    verbose=verbose,
                )

            if "cellrank_pseudotime" in adata.obs:
                try:
                    log_message(
                        "Plotting gene expression trends...",
                        message_type="info",
                        verbose=verbose,
                    )

                    driver_key = "to_terminal_states_lineage_drivers"
                    if driver_key in adata.varm:
                        drivers_df = adata.varm[driver_key]
                        top_genes = drivers_df.nlargest(
                            20, columns=drivers_df.columns[0]
                        ).index.tolist()
                        log_message(
                            f"Using top {len(top_genes)} driver genes",
                            message_type="info",
                            verbose=verbose,
                        )

                        model = cr.models.GAM()
                        cr.pl.gene_trends(
                            adata,
                            model=model,
                            genes=top_genes[:20],
                            time_key="cellrank_pseudotime",
                            n_jobs=n_jobs,
                            dpi=plot_dpi,
                            save=get_save_path("gene_trends"),
                        )
                        if show_plot:
                            plt.show()
                        plt.close()
                    else:
                        log_message(
                            "No lineage drivers found, skipping gene trends",
                            message_type="warning",
                            verbose=verbose,
                        )
                except Exception as e:
                    log_message(
                        f"Failed to plot gene trends: {e}",
                        message_type="warning",
                        verbose=verbose,
                    )

            try:
                log_message(
                    "Plotting kernel projection...",
                    message_type="info",
                    verbose=verbose,
                )
                proj_basis = basis
                log_message(
                    f"Using basis: {proj_basis} for projection",
                    message_type="info",
                    verbose=verbose,
                )
                if final_kernel is not None:
                    log_message(
                        f"Final kernel type: {type(final_kernel).__name__}",
                        message_type="info",
                        verbose=verbose,
                    )

                    plot_palette = palette
                    if plot_palette is None and group_by is not None:
                        groups = (
                            adata.obs[group_by].cat.categories
                            if hasattr(adata.obs[group_by], "cat")
                            else adata.obs[group_by].unique()
                        )
                        plot_palette = dict(
                            zip(groups, plt.cm.tab10(np.linspace(0, 1, len(groups))))
                        )
                        log_message(
                            f"Created default palette for {len(groups)} groups",
                            message_type="info",
                            verbose=verbose,
                        )

                    plot_kwargs = {}
                    if save_plot:
                        ext = (
                            plot_format
                            if plot_format in ["png", "pdf", "svg"]
                            else "png"
                        )
                        save_path = os.path.abspath(
                            f"{plot_prefix}_projection_{proj_basis}.{ext}"
                        )
                        plot_kwargs["save"] = save_path
                        log_message(
                            f"Will save projection to: {save_path}",
                            message_type="info",
                            verbose=verbose,
                        )
                    if plot_dpi:
                        plot_kwargs["dpi"] = plot_dpi
                        log_message(
                            f"Using DPI: {plot_dpi}",
                            message_type="info",
                            verbose=verbose,
                        )

                    if group_by is not None:
                        plot_kwargs["color"] = group_by
                        if plot_palette is not None:
                            plot_kwargs["palette"] = plot_palette
                        plot_kwargs["legend_loc"] = legend_loc
                        log_message(
                            f"Adding color information: color={group_by}, palette={'provided' if palette is not None else 'default'}",
                            message_type="info",
                            verbose=verbose,
                        )

                    log_message(
                        "Calling final_kernel.plot_projection()...",
                        message_type="info",
                        verbose=verbose,
                    )
                    try:
                        final_kernel.plot_projection(
                            basis=proj_basis,
                            recompute=True,
                            **plot_kwargs,
                        )
                        if show_plot:
                            import matplotlib.pyplot as plt

                            plt.show()
                        if save_plot and plot_format == "pdf":
                            save_path = plot_kwargs.get("save")
                            if save_path and save_path.endswith(".pdf"):
                                png_path = save_path.replace(".pdf", ".png")
                                if os.path.exists(png_path) and not os.path.exists(
                                    save_path
                                ):
                                    if os.path.exists(save_path):
                                        try:
                                            os.remove(save_path)
                                        except Exception:
                                            pass
                                    log_message(
                                        f"PDF save failed, CellRank saved as PNG: {png_path}",
                                        message_type="warning",
                                        verbose=verbose,
                                    )
                        log_message(
                            "Projection plot completed successfully",
                            message_type="success",
                            verbose=verbose,
                        )
                    except Exception as e:
                        log_message(
                            f"Failed to plot projection: {e}",
                            message_type="warning",
                            verbose=verbose,
                        )
                else:
                    log_message(
                        "Cannot plot projection: kernel object is None",
                        message_type="warning",
                        verbose=verbose,
                    )
            except Exception as e:
                log_message(
                    f"Failed to plot projection: {e}",
                    message_type="warning",
                    verbose=verbose,
                )
                import traceback

                log_message(
                    f"Projection error traceback: {traceback.format_exc()}",
                    message_type="warning",
                    verbose=verbose,
                )

        log_message("Visualization completed!", message_type="success", verbose=verbose)

        if use_velocity:
            try:
                if "dynamical" not in mode:
                    log_message(
                        "Running recover_dynamics for latent time computation...",
                        message_type="info",
                        verbose=verbose,
                    )
                    try:
                        vkey = mode[-1] if isinstance(mode, list) else mode
                        velocity_graph_key = vkey + "_graph"

                        if velocity_graph_key in adata.uns:
                            adata.uns["velocity_graph"] = adata.uns[velocity_graph_key]
                            if velocity_graph_key + "_neg" in adata.uns:
                                adata.uns["velocity_graph_neg"] = adata.uns[
                                    velocity_graph_key + "_neg"
                                ]
                            else:
                                log_message(
                                    "velocity_graph_neg not found, computing velocity graph...",
                                    message_type="info",
                                    verbose=verbose,
                                )
                                scv.tl.velocity_graph(
                                    adata,
                                    vkey=vkey,
                                    n_neighbors=n_neighbors,
                                    n_jobs=n_jobs,
                                )
                        elif "velocity_graph" not in adata.uns:
                            log_message(
                                "velocity_graph not found, computing velocity graph...",
                                message_type="info",
                                verbose=verbose,
                            )
                            scv.tl.velocity_graph(
                                adata, vkey=vkey, n_neighbors=n_neighbors, n_jobs=n_jobs
                            )

                        scv.tl.recover_dynamics(adata, use_raw=False, n_jobs=n_jobs)
                        terminal_key = (
                            "to_terminal_states_probs"
                            if hasattr(estimator, "terminal_states_probabilities")
                            else "terminal_states_probs"
                        )
                        initial_key = (
                            "to_initial_states_probs"
                            if hasattr(estimator, "initial_states_probabilities")
                            else "initial_states_probs"
                        )

                        scv.tl.recover_latent_time(
                            adata,
                            root_key=initial_key,
                            end_key=terminal_key,
                        )
                        if "latent_time" in adata.obs:
                            adata.obs["cellrank_pseudotime"] = adata.obs["latent_time"]
                        log_message(
                            "Latent time computed successfully",
                            message_type="success",
                            verbose=verbose,
                        )
                    except Exception as e:
                        log_message(
                            f"recover_dynamics failed ({e}), skipping latent time computation...",
                            message_type="warning",
                            verbose=verbose,
                        )
                else:
                    terminal_key = (
                        "to_terminal_states_probs"
                        if hasattr(estimator, "terminal_states_probabilities")
                        else "terminal_states_probs"
                    )
                    initial_key = (
                        "to_initial_states_probs"
                        if hasattr(estimator, "initial_states_probabilities")
                        else "initial_states_probs"
                    )

                    scv.tl.recover_latent_time(
                        adata, root_key=initial_key, end_key=terminal_key
                    )
                    if "latent_time" in adata.obs:
                        adata.obs["cellrank_pseudotime"] = adata.obs["latent_time"]
                    log_message(
                        "Latent time computed successfully",
                        message_type="success",
                        verbose=verbose,
                    )
            except Exception as e:
                log_message(
                    f"Latent time computation failed: {e}",
                    message_type="warning",
                    verbose=verbose,
                )
        else:
            log_message(
                "Skipping latent time computation (no velocity data available)",
                message_type="info",
                verbose=verbose,
            )

    finally:
        try:
            figures_dir = os.path.join(os.getcwd(), "figures")
            if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
                if not os.listdir(figures_dir):
                    os.rmdir(figures_dir)
                    log_message(
                        f"Removed empty figures directory: {figures_dir}",
                        message_type="info",
                        verbose=verbose,
                    )
        except Exception:
            pass

        if save_plot:
            os.chdir(prevdir)

    try:
        adata.__dict__["_raw"].__dict__["_var"] = (
            adata.__dict__["_raw"]
            .__dict__["_var"]
            .rename(columns={"_index": "features"})
        )
    except:
        pass

    log_message(
        "Returning adata, estimator, and kernel objects",
        message_type="info",
        verbose=verbose,
    )
    return adata, estimator, final_kernel


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
    n_jobs=1,
    show_plot=True,
    save_plot=False,
    plot_format="png",
    plot_dpi=600,
    plot_prefix="paga",
    dirpath="./paga",
    dpi=300,
    save=False,
    fileprefix="",
    verbose=True,
    legend_loc="on data",
):
    import os
    import platform
    import sys

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"
    if is_apple_silicon:
        os.environ["PYTHONHASHSEED"] = "0"
        os.environ["PYTHONUNBUFFERED"] = "1"
        os.environ["SCANPY_SETTINGS"] = "scanpy_settings"
        os.environ["SCANPY_SETTINGS_VERBOSITY"] = "1"
        os.environ["MPLBACKEND"] = "Agg"
        os.environ["DISPLAY"] = ""
        os.environ["NUMBA_NUM_THREADS"] = "1"
        os.environ["NUMBA_DISABLE_JIT"] = "1"
        os.environ["NUMBA_THREADING_LAYER"] = "tbb"
        os.environ["NUMBA_DEFAULT_NUM_THREADS"] = "1"

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        log_message(
            f"matplotlib setup failed: {e}", message_type="warning", verbose=verbose
        )
        plt = None

    if is_apple_silicon:
        try:
            import numba

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

    prevdir = os.getcwd()

    expanded_path = os.path.expanduser(dirpath)
    os.makedirs(expanded_path, exist_ok=True)
    os.chdir(expanded_path)

    import scanpy as sc

    sc.settings.figdir = "."

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
            log_message("adata or h5ad must be provided", message_type="error")
            return None

        if adata is None:
            adata = sc.read(h5ad)

        if group_by is None:
            log_message("group_by must be provided", message_type="error")
            return None

        if linear_reduction is None and nonlinear_reduction is None:
            log_message(
                "linear_reduction or nonlinear_reduction must be provided at least one",
                message_type="error",
            )
            return None

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
                "root_cell or root_group should be provided", message_type="error"
            )
            return None

        if use_rna_velocity is True:
            adata.uns["velocity_graph"] = adata.uns[vkey + "_graph"]

        adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")

        rep_key = linear_reduction
        try:
            obsm_keys = set(adata.obsm_keys())
            if linear_reduction not in obsm_keys:
                log_message(
                    f"linear_reduction '{linear_reduction}' not found in adata.obsm",
                    message_type="error",
                    verbose=verbose,
                )
                return None
            rep_key = linear_reduction

            neighbors_n_pcs = n_pcs
            if rep_key in obsm_keys:
                rep_dims = adata.obsm[rep_key].shape[1]
                if neighbors_n_pcs is not None and isinstance(
                    neighbors_n_pcs, (int, float)
                ):
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

        sc.tl.paga(adata, groups=group_by, use_rna_velocity=use_rna_velocity)

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
                    "./" + ".".join(filter(None, [fileprefix, "paga_compare.pdf"])),
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
                    "./" + ".".join(filter(None, [fileprefix, "paga_layout.pdf"])),
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
                legend_loc=legend_loc,
                show=show_plot,
            )

            if save:
                plt.savefig(
                    "./" + ".".join(filter(None, [fileprefix, "paga_graph.pdf"])),
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
                        "./" + ".".join(filter(None, [fileprefix, "paga_umap.pdf"])),
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
                        "./"
                        + ".".join(filter(None, [fileprefix, "dpt_pseudotime.pdf"])),
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
        try:
            figures_dir = os.path.join(os.getcwd(), "figures")
            if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
                if not os.listdir(figures_dir):
                    os.rmdir(figures_dir)
                    log_message(
                        f"Removed empty figures directory: {figures_dir}",
                        message_type="info",
                        verbose=verbose,
                    )
        except Exception:
            pass

        os.chdir(prevdir)

        if is_apple_silicon and plt is not None:
            try:
                plt.clf()
                plt.close("all")
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
    legend_loc="on data",
):
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

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

        numba.config.DISABLE_JIT = True
        numba.set_num_threads(1)
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
    expanded_path = os.path.expanduser(dirpath)
    os.makedirs(expanded_path, exist_ok=True)
    os.chdir(expanded_path)

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
            log_message("adata or h5ad must be provided", message_type="error")
            exit()

        if adata is None:
            adata = sc.read(h5ad)

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

        # if linear_reduction is None:
        #     sc.pp.pca(adata, n_comps=n_pcs)
        #     linear_reduction = "X_pca"

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
                "./" + ".".join(filter(None, [fileprefix, "palantir_pseudotime.pdf"])),
                dpi=dpi,
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
                "./"
                + ".".join(filter(None, [fileprefix, "palantir_diff_potential.pdf"])),
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
                "./" + ".".join(filter(None, [fileprefix, "palantir_probs.pdf"])),
                dpi=dpi,
            )

        if group_by is not None:
            sc.pl.embedding(
                adata,
                basis=basis,
                color=group_by,
                size=point_size,
                palette=palette,
                legend_loc=legend_loc,
                show=show_plot,
            )
            if save:
                plt.savefig(
                    "./"
                    + ".".join(filter(None, [fileprefix, "palantir_group_by.pdf"])),
                    dpi=dpi,
                )

    finally:
        try:
            figures_dir = os.path.join(os.getcwd(), "figures")
            if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
                if not os.listdir(figures_dir):
                    os.rmdir(figures_dir)
                    log_message(
                        f"Removed empty figures directory: {figures_dir}",
                        message_type="info",
                        verbose=verbose,
                    )
        except Exception:
            pass

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
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

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

        numba.config.DISABLE_JIT = True
        numba.set_num_threads(1)
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
    expanded_path = os.path.expanduser(dirpath)
    os.makedirs(expanded_path, exist_ok=True)
    os.chdir(expanded_path)

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
            log_message("`adata` or `h5ad` must be provided", message_type="error")
            exit()

        if adata is None:
            adata = sc.read(h5ad)

        if group_by is None:
            log_message("`group_by` must be provided", message_type="error")
            exit()

        if time_field is None:
            log_message("`time_field` must be provided", message_type="error")
            exit()

        adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")
        if pd.api.types.is_categorical_dtype(adata.obs[time_field]):
            adata.obs["time_field"] = adata.obs[time_field].cat.codes
        elif not pd.api.types.is_numeric_dtype(adata.obs[time_field]):
            try:
                adata.obs["time_field"] = adata.obs[time_field].astype("float")
            except ValueError:
                log_message(
                    f"Unable to convert column `{time_field}` to float type",
                    message_type="warning",
                )
        else:
            adata.obs["time_field"] = adata.obs[time_field]

        time_dict = dict(zip(adata.obs[time_field], adata.obs["time_field"]))
        if time_from not in time_dict.keys():
            log_message("`time_from` is incorrect", message_type="error")
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
        try:
            figures_dir = os.path.join(os.getcwd(), "figures")
            if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
                if not os.listdir(figures_dir):
                    os.rmdir(figures_dir)
                    log_message(
                        f"Removed empty figures directory: {figures_dir}",
                        message_type="info",
                        verbose=verbose,
                    )
        except Exception:
            pass

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


def CellTypistModels(on_the_fly=False, verbose=True):
    """
    Get a list of all available CellTypist models.

    Parameters
    ----------
    on_the_fly : bool
        If True, only show downloaded models.
        If False, show all available models (fetch list from server).
    verbose : bool
        Whether to show detailed information.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing model names and descriptions.
    """
    import os
    import platform

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"
    if is_apple_silicon:
        os.environ["PYTHONHASHSEED"] = "0"
        os.environ["PYTHONUNBUFFERED"] = "1"
        os.environ["MPLBACKEND"] = "Agg"
        os.environ["DISPLAY"] = ""

    try:
        import celltypist
        from celltypist import models
    except ImportError as e:
        log_message(
            f"celltypist import failed: {e}",
            message_type="error",
            verbose=verbose,
        )
        raise

    try:
        models_df = models.models_description(on_the_fly=on_the_fly)
        return models_df
    except Exception as e:
        log_message(
            f"Failed to get model list: {e}",
            message_type="error",
            verbose=verbose,
        )
        raise


def CellTypist(
    adata=None,
    h5ad=None,
    model="Immune_All_Low.pkl",
    mode="best match",
    p_thres=0.5,
    majority_voting=False,
    over_clustering=None,
    min_prop=0,
    use_GPU=False,
    insert_labels=True,
    insert_conf=True,
    insert_conf_by="predicted_labels",
    insert_prob=False,
    insert_decision=False,
    prefix="",
    verbose=True,
):
    """
    Run CellTypist cell type annotation.

    Parameters
    ----------
    adata : AnnData
        AnnData object (required if h5ad is not provided).
    h5ad : str
        Path to h5ad file (optional).
    model : str
        Model name or path. Default is 'Immune_All_Low.pkl'.
        Supports three formats:
        1. Model name (e.g., 'Immune_All_Low.pkl'): automatically searched in ~/.celltypist/data/models/
        2. Full path (contains '/'): use the provided path directly
        3. None: use default model via celltypist.models.get_default_model()
    mode : str
        Prediction mode: 'best match' or 'prob match'. Default is 'best match'.
    p_thres : float
        Probability threshold for 'prob match' mode. Default is 0.5.
    majority_voting : bool
        Whether to use majority voting. Default is False.
    over_clustering : str, list, or None
        Over-clustering result. Can be:
        - String: column name in adata.obs
        - List/array: over-clustering labels
        - None: use heuristic over-clustering
    min_prop : float
        Minimum proportion for majority voting. Default is 0.
    use_GPU : bool
        Whether to use GPU for over-clustering. Default is False.
    insert_labels : bool
        Whether to insert predicted labels into AnnData. Default is True.
    insert_conf : bool
        Whether to insert confidence scores. Default is True.
    insert_conf_by : str
        Which prediction type to base confidence on. Default is 'predicted_labels'.
    insert_prob : bool
        Whether to insert probability matrix. Default is False.
    insert_decision : bool
        Whether to insert decision matrix. Default is False.
    prefix : str
        Prefix for inserted columns. Default is empty string.
    verbose : bool
        Whether to show detailed information. Default is True.

    Returns
    -------
    AnnData
        AnnData object with CellTypist predictions inserted.
    """
    import os
    import platform
    import sys

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["KMP_WARNINGS"] = "0"
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"
    if is_apple_silicon:
        os.environ["PYTHONHASHSEED"] = "0"
        os.environ["PYTHONUNBUFFERED"] = "1"
        os.environ["SCANPY_SETTINGS"] = "scanpy_settings"
        os.environ["SCANPY_SETTINGS_VERBOSITY"] = "1"
        os.environ["MPLBACKEND"] = "Agg"
        os.environ["DISPLAY"] = ""
        os.environ["NUMBA_NUM_THREADS"] = "1"
        os.environ["NUMBA_DISABLE_JIT"] = "1"
        os.environ["NUMBA_THREADING_LAYER"] = "tbb"
        os.environ["NUMBA_DEFAULT_NUM_THREADS"] = "1"

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        log_message(
            f"matplotlib setup failed: {e}", message_type="warning", verbose=verbose
        )
        plt = None

    if is_apple_silicon:
        try:
            import numba

            numba.config.DISABLE_JIT = True
            numba.set_num_threads(1)
        except:
            pass

    try:
        import celltypist
        from celltypist import models
    except ImportError as e:
        log_message(
            f"celltypist import failed: {e}",
            message_type="error",
            verbose=verbose,
        )
        raise

    import warnings

    warnings.simplefilter("ignore", category=UserWarning)
    warnings.simplefilter("ignore", category=FutureWarning)
    warnings.simplefilter("ignore", category=DeprecationWarning)

    try:
        if adata is None and h5ad is None:
            log_message(
                "`adata` or `h5ad` must be provided",
                message_type="error",
                verbose=verbose,
            )
            return None

        if adata is None:
            import scanpy as sc

            adata = sc.read(h5ad)

        import numpy as np
        import scanpy as sc

        needs_preprocessing = False

        if adata.X is not None:
            sample_size = min(100, adata.n_obs)
            try:
                if hasattr(adata.X, "toarray"):
                    sample_data = adata.X[:sample_size].toarray()
                else:
                    sample_data = adata.X[:sample_size]

                min_val = float(np.min(sample_data))
                max_val = float(np.max(sample_data))

                if min_val < 0 or max_val > 9.22:
                    needs_preprocessing = True
                    log_message(
                        f"Data format invalid for CellTypist (min={min_val:.2f}, max={max_val:.2f}). "
                        f"Will attempt to preprocess...",
                        message_type="warning",
                        verbose=verbose,
                    )
            except Exception as e:
                log_message(
                    f"Could not check data format: {e}. Attempting preprocessing...",
                    message_type="warning",
                    verbose=verbose,
                )
                needs_preprocessing = True

        if needs_preprocessing:
            log_message(
                "Preprocessing data for CellTypist (normalize to 10000 counts per cell and log1p transform)...",
                message_type="info",
                verbose=verbose,
            )

            if adata.raw is not None:
                log_message(
                    "Using .raw data for preprocessing",
                    message_type="info",
                    verbose=verbose,
                )
                adata = adata.raw.to_adata()
            elif adata.X is not None:
                try:
                    if hasattr(adata.X, "toarray"):
                        sample_check = adata.X[: min(10, adata.n_obs)].toarray()
                    else:
                        sample_check = adata.X[: min(10, adata.n_obs)]

                    if np.all(sample_check >= 0) and np.max(sample_check) > 20:
                        adata.raw = adata.copy()
                        log_message(
                            "Detected raw counts in .X, storing in .raw",
                            message_type="info",
                            verbose=verbose,
                        )
                except:
                    pass

            sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
            sc.pp.log1p(adata)

        if model is None:
            model = models.get_default_model()
            log_message(
                f"Using default model: {model}",
                message_type="info",
                verbose=verbose,
            )
        elif "/" not in model:
            try:
                all_models = models.get_all_models()
                if model in all_models:
                    model = models.get_model_path(model)
                    log_message(
                        f"Using model: {model}",
                        message_type="info",
                        verbose=verbose,
                    )
                else:
                    log_message(
                        f"Model '{model}' not found in downloaded models. "
                        f"Will try to use as path or download if needed.",
                        message_type="warning",
                        verbose=verbose,
                    )
            except Exception as e:
                log_message(
                    f"Could not check model list: {e}. Using model as provided.",
                    message_type="warning",
                    verbose=verbose,
                )

        over_clustering_value = None
        if over_clustering is not None:
            if isinstance(over_clustering, str):
                if over_clustering in adata.obs.columns:
                    over_clustering_value = adata.obs[over_clustering].values
                else:
                    log_message(
                        f"'{over_clustering}' not found in adata.obs. "
                        f"Will use heuristic over-clustering.",
                        message_type="warning",
                        verbose=verbose,
                    )
            else:
                over_clustering_value = over_clustering

        log_message("Running CellTypist annotation...", verbose=verbose)

        predictions = celltypist.annotate(
            filename=adata,
            model=model,
            mode=mode,
            p_thres=p_thres,
            majority_voting=majority_voting,
            over_clustering=over_clustering_value,
            min_prop=min_prop,
            use_GPU=use_GPU,
        )

        adata = predictions.to_adata(
            insert_labels=insert_labels,
            insert_conf=insert_conf,
            insert_conf_by=insert_conf_by,
            insert_decision=insert_decision,
            insert_prob=insert_prob,
            prefix=prefix,
        )

        log_message(
            "CellTypist annotation completed successfully",
            message_type="success",
            verbose=verbose,
        )

        return adata

    except Exception as e:
        log_message(
            f"CellTypist annotation failed: {e}",
            message_type="error",
            verbose=verbose,
        )
        raise
