//! Rounding operations for Dilithium.

const params = @import("params.zig");

const Q = params.Q;
const D = params.D;
const GAMMA2 = params.GAMMA2;

/// Result of power2round operation.
pub const Power2RoundResult = struct {
    a0: i32,
    a1: i32,
};

/// For finite field element a, compute a0, a1 such that
/// a mod^+ Q = a1*2^D + a0 with -2^{D-1} < a0 <= 2^{D-1}.
/// Assumes a to be standard representative.
pub fn power2round(a: i32) Power2RoundResult {
    const a1 = (a + (1 << (D - 1)) - 1) >> D;
    const a0 = a - (a1 << D);
    return .{ .a0 = a0, .a1 = a1 };
}

/// Result of decompose operation.
pub const DecomposeResult = struct {
    a0: i32,
    a1: i32,
};

/// For finite field element a, compute high and low bits a0, a1 such that
/// a mod^+ Q = a1*ALPHA + a0 with -ALPHA/2 < a0 <= ALPHA/2 except
/// if a1 = (Q-1)/ALPHA where we set a1 = 0 and
/// -ALPHA/2 <= a0 = a mod^+ Q - Q < 0.
/// Assumes a to be standard representative.
pub fn decompose(a: i32) DecomposeResult {
    var a1: i32 = (a + 127) >> 7;

    if (GAMMA2 == (Q - 1) / 32) {
        a1 = (a1 *% 1025 + (1 << 21)) >> 22;
        a1 &= 15;
    } else if (GAMMA2 == (Q - 1) / 88) {
        a1 = (a1 *% 11275 + (1 << 23)) >> 24;
        a1 ^= ((43 - a1) >> 31) & a1;
    }

    var a0: i32 = a - a1 *% (2 * GAMMA2);
    a0 -= (((Q - 1) / 2 - a0) >> 31) & Q;

    return .{ .a0 = a0, .a1 = a1 };
}

/// Compute hint bit indicating whether the low bits of the
/// input element overflow into the high bits.
///
/// Returns true if overflow.
pub fn makeHint(a0: i32, a1: i32) bool {
    // Constant-time implementation to avoid data-dependent branches
    const c1 = @intFromBool(a0 > GAMMA2);
    const c2 = @intFromBool(a0 < -GAMMA2);
    const c3 = @intFromBool(a0 == -GAMMA2) & @intFromBool(a1 != 0);
    return (c1 | c2 | c3) != 0;
}

/// Correct high bits according to hint.
///
/// Returns corrected high bits.
pub fn useHint(a: i32, hint: bool) i32 {
    const result = decompose(a);
    const a0 = result.a0;
    const a1 = result.a1;

    if (!hint) return a1;

    if (GAMMA2 == (Q - 1) / 32) {
        return if (a0 > 0)
            (a1 + 1) & 15
        else
            (a1 - 1) & 15;
    } else if (GAMMA2 == (Q - 1) / 88) {
        return if (a0 > 0)
            (if (a1 == 43) @as(i32, 0) else a1 + 1)
        else
            (if (a1 == 0) @as(i32, 43) else a1 - 1);
    }

    unreachable;
}

// Optional: C-compatible API if needed for interop
pub fn power2roundPtr(a0: *i32, a: i32) i32 {
    const result = power2round(a);
    a0.* = result.a0;
    return result.a1;
}

pub fn decomposePtr(a0: *i32, a: i32) i32 {
    const result = decompose(a);
    a0.* = result.a0;
    return result.a1;
}
