#!/usr/bin/env python3

import sys

last_line_was_latex = True

for line in sys.stdin:
    if line.startswith("//"):
        latex_line = line[2:]
        if not last_line_was_latex:
            print("\\end{verbatim}")
        print(latex_line if latex_line.isspace() else latex_line.lstrip(), end="")
        last_line_was_latex = True
    else:
        if last_line_was_latex:
            print("\\begin{verbatim}")
        print(line, end="")
        last_line_was_latex = False
