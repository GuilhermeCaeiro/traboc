const output = IOBuffer()

using REPL

const out_terminal = REPL.Terminals.TerminalBuffer(output)

const basic_repl = REPL.BasicREPL(out_terminal)

const basic_display = REPL.REPLDisplay(basic_repl)

Base.pushdisplay(basic_display)