#!/bin/sh
# Locate the rust-lld bundled with the active toolchain and forward all arguments.
exec "$(rustc --print sysroot)/lib/rustlib/$(rustc -vV | sed -n 's/^host: //p')/bin/rust-lld" "$@"
