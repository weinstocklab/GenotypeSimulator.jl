#!/usr/bin/env bash

JULIAC=$(julia +1.12 -e 'print(joinpath(Sys.BINDIR, "..", "share", "julia", "juliac", "juliac.jl"))')

julia +1.12 --project=. --experimental ${JULIAC} --verbose --output-exe simulator --experimental --trim=unsafe-warn scripts/main.jl